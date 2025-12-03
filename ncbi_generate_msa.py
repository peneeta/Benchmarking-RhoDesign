import requests
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from Bio import AlignIO
import argparse
import logging
import subprocess
import os

"""
Generate MSAs using NCBI BLAST API for RNA sequences
Creates MSA in .a3m format which is required for RhoFold+
Now includes MAFFT alignment for proper A3M format
REQUIRES INSTALLING MAFFT AND BIOPYTHON!!
"""

# set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# dummy entrez email setting
Entrez.email = "your.email@example.com"

class NCBIBlastMSAGenerator:
    """Generate MSAs using NCBI BLAST API with proper alignment."""
    
    def __init__(self, output_dir, email, min_identity=70, max_hits=100, database="nt"):
        """
        Initialize the BLAST MSA generator.
        
        Args:
            output_dir: Directory to save output files
            email: Your email (required by NCBI)
            min_identity: Minimum sequence identity threshold (default: 70%)
            max_hits: Maximum number of hits to retrieve (default: 100)
            database: NCBI database to search (default: "nt" for nucleotide)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_identity = min_identity
        self.max_hits = max_hits
        self.database = database
        Entrez.email = email
        self.base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        
    def submit_blast_search(self, sequence):
        """
        Submit a BLAST search to NCBI.
        
        Args:
            sequence: RNA sequence string
            
        Returns:
            Request ID (RID) for the search
        """
        params = {
            "CMD": "Put",
            "PROGRAM": "blastn",
            "DATABASE": self.database,
            "QUERY": sequence,
            "HITLIST_SIZE": self.max_hits,
            "FORMAT_TYPE": "XML"
        }
        
        logger.info("Submitting BLAST search to NCBI...")
        
        try:
            response = requests.post(self.base_url, data=params, timeout=30)
            response.raise_for_status()
            
            # Parse the RID from response
            content = response.text
            for line in content.split('\n'):
                if 'RID =' in line:
                    rid = line.split('=')[1].strip()
                    logger.info(f"BLAST search submitted. RID: {rid}")
                    return rid
            
            logger.error("Could not find RID in response")
            return None
            
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to submit BLAST search: {e}")
            return None
    
    def check_blast_status(self, rid, max_wait=600, poll_interval=30):
        """
        Check the status of a BLAST search.
        
        Args:
            rid: Request ID from submission
            max_wait: Maximum time to wait in seconds (default: 600)
            poll_interval: Time between status checks in seconds (default: 30)
            
        Returns:
            True if completed, False otherwise
        """
        params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_OBJECT": "SearchInfo"
        }
        
        start_time = time.time()
        
        while time.time() - start_time < max_wait:
            try:
                response = requests.get(self.base_url, params=params, timeout=30)
                response.raise_for_status()
                
                content = response.text
                
                if "Status=READY" in content:
                    logger.info(f"BLAST search {rid} completed")
                    return True
                elif "Status=WAITING" in content:
                    logger.info(f"BLAST search {rid} still running...")
                    time.sleep(poll_interval)
                else:
                    logger.error(f"BLAST search {rid} failed")
                    return False
                    
            except requests.exceptions.RequestException as e:
                logger.error(f"Error checking BLAST status: {e}")
                return False
        
        logger.error(f"BLAST search {rid} timed out after {max_wait} seconds")
        return False
    
    def get_blast_results(self, rid):
        """
        Retrieve BLAST search results.
        
        Args:
            rid: Request ID from submission
            
        Returns:
            List of hit dictionaries with sequence info
        """
        params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "XML"
        }
        
        try:
            response = requests.get(self.base_url, params=params, timeout=60)
            response.raise_for_status()
            
            # Parse XML results
            root = ET.fromstring(response.content)
            
            hits = []
            for iteration in root.findall('.//Iteration'):
                for hit in iteration.findall('.//Hit'):
                    hit_id = hit.find('Hit_id').text
                    hit_def = hit.find('Hit_def').text
                    
                    # Get HSP (High-scoring Segment Pair) info
                    for hsp in hit.findall('.//Hsp'):
                        identity = float(hsp.find('Hsp_identity').text)
                        align_len = float(hsp.find('Hsp_align-len').text)
                        hit_seq = hsp.find('Hsp_hseq').text
                        
                        pct_identity = (identity / align_len) * 100
                        
                        if pct_identity >= self.min_identity:
                            hits.append({
                                'id': hit_id,
                                'description': hit_def,
                                'sequence': hit_seq.replace('-', ''),
                                'identity': pct_identity
                            })
                        
                        break  # Only take first HSP per hit
            
            logger.info(f"Retrieved {len(hits)} hits with identity >= {self.min_identity}%")
            return hits[:self.max_hits]
            
        except Exception as e:
            logger.error(f"Failed to retrieve BLAST results: {e}")
            return []
    
    def align_sequences_with_mafft(self, query_seq, query_id, hits, temp_dir):
        """
        Align sequences using MAFFT.
        
        Args:
            query_seq: Query sequence
            query_id: Query sequence ID
            hits: List of hit dictionaries
            temp_dir: Directory for temporary files
            
        Returns:
            List of aligned SeqRecord objects or None if alignment fails
        """
        if len(hits) == 0:
            logger.warning("No hits to align, returning query only")
            return [SeqRecord(Seq(str(query_seq)), id=query_id, description="")]
        
        # Create temporary FASTA file with all sequences
        temp_fasta = temp_dir / "temp_unaligned.fasta"
        temp_aligned = temp_dir / "temp_aligned.fasta"
        
        try:
            # Write unaligned sequences
            with open(temp_fasta, 'w') as f:
                f.write(f">{query_id}\n{str(query_seq)}\n")
                for idx, hit in enumerate(hits):
                    hit_id = hit['id'].replace(' ', '_')
                    f.write(f">{hit_id}\n{hit['sequence']}\n")
            
            # Run MAFFT
            logger.info("Running MAFFT alignment...")
            mafft_cmd = f"mafft --auto --quiet --thread -1 {temp_fasta} > {temp_aligned}"
            result = subprocess.run(
                mafft_cmd, 
                shell=True, 
                capture_output=True, 
                text=True
            )
            
            if result.returncode != 0:
                logger.error(f"MAFFT failed: {result.stderr}")
                return None
            
            # Read aligned sequences
            aligned_seqs = list(SeqIO.parse(temp_aligned, "fasta"))
            logger.info(f"Alignment complete. Length: {len(aligned_seqs[0].seq)}")
            
            return aligned_seqs
            
        except FileNotFoundError:
            logger.error("MAFFT not found. Please install: conda install -c bioconda mafft")
            return None
        except Exception as e:
            logger.error(f"Alignment failed: {e}")
            return None
        finally:
            # Clean up temp files
            if temp_fasta.exists():
                temp_fasta.unlink()
            if temp_aligned.exists():
                temp_aligned.unlink()
    
    def convert_to_a3m(self, aligned_seqs):
        """
        Convert aligned sequences to A3M format.
        A3M format removes positions where the query has gaps.
        
        Args:
            aligned_seqs: List of aligned SeqRecord objects (query must be first)
            
        Returns:
            List of A3M formatted sequences
        """
        if not aligned_seqs:
            return []
        
        query_seq = str(aligned_seqs[0].seq)
        a3m_seqs = []
        
        for record in aligned_seqs:
            seq = str(record.seq)
            a3m_seq = []
            
            # Keep only positions where query has a residue (not a gap)
            for q_char, s_char in zip(query_seq, seq):
                if q_char != '-':
                    a3m_seq.append(s_char.upper())
            
            a3m_seqs.append({
                'id': record.id,
                'sequence': ''.join(a3m_seq)
            })
        
        return a3m_seqs
    
    def create_msa_from_hits(self, query_seq, query_id, hits, output_file, temp_dir):
        """
        Create an aligned MSA in A3M format from BLAST hits.
        
        Args:
            query_seq: Query sequence
            query_id: Query sequence ID
            hits: List of hit dictionaries
            output_file: Path to output A3M file
            temp_dir: Directory for temporary files
            
        Returns:
            Number of sequences in MSA
        """
        # Align sequences with MAFFT
        aligned_seqs = self.align_sequences_with_mafft(
            query_seq, query_id, hits, temp_dir
        )
        
        if aligned_seqs is None:
            logger.warning("Alignment failed, creating unaligned MSA")
            # Fallback to unaligned format
            with open(output_file, 'w') as f:
                f.write(f">{query_id}\n{str(query_seq)}\n")
                for idx, hit in enumerate(hits):
                    hit_id = hit['id'].replace(' ', '_')
                    f.write(f">{hit_id}\n{hit['sequence']}\n")
            return len(hits) + 1
        
        # Convert to A3M format
        a3m_seqs = self.convert_to_a3m(aligned_seqs)
        
        # Write A3M file
        with open(output_file, 'w') as f:
            for i, seq_data in enumerate(a3m_seqs):
                f.write(f">{seq_data['id']}\n")
                f.write(f"{seq_data['sequence']}\n")
        
        num_seqs = len(a3m_seqs)
        logger.info(f"Created aligned A3M with {num_seqs} sequences: {output_file}")
        
        return num_seqs
    
    def process_sequence(self, seq_record, seq_index, pdb_id, pdb_output_dir):
        """
        Process a single sequence through BLAST search and alignment.
        
        Args:
            seq_record: BioPython SeqRecord object
            seq_index: Index of the sequence
            pdb_id: PDB ID for naming
            pdb_output_dir: Output directory for this PDB
            
        Returns:
            Dictionary with results
        """
        output_prefix = f"{pdb_id}_seq_{seq_index}"
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing {pdb_id} sequence {seq_index}: {seq_record.id}")
        logger.info(f"Sequence length: {len(seq_record.seq)}")
        logger.info(f"{'='*60}")
        
        result = {
            "pdb_id": pdb_id,
            "sequence_id": seq_record.id,
            "index": seq_index,
            "status": "failed",
            "msa_file": None,
            "num_sequences": 0,
            "error": None
        }
        
        # Create temp directory for this sequence
        temp_dir = pdb_output_dir / f"temp_{output_prefix}"
        temp_dir.mkdir(exist_ok=True)
        
        try:
            # Submit BLAST search
            rid = self.submit_blast_search(str(seq_record.seq))
            
            if not rid:
                result["error"] = "Failed to submit BLAST search"
                return result
            
            # Wait for completion
            if not self.check_blast_status(rid):
                result["error"] = "BLAST search failed or timed out"
                return result
            
            # Get results
            hits = self.get_blast_results(rid)
            
            # Create aligned MSA
            output_a3m = pdb_output_dir / f"{output_prefix}.a3m"
            num_seqs = self.create_msa_from_hits(
                seq_record.seq,
                seq_record.id,
                hits,
                output_a3m,
                temp_dir
            )
            
            result["msa_file"] = str(output_a3m)
            result["num_sequences"] = num_seqs
            
            if num_seqs == 1:
                result["status"] = "no_hits"
                result["error"] = "No homologous sequences found"
                logger.warning(f"No hits found for {output_prefix}, created single-sequence MSA")
            elif num_seqs < 10:
                result["status"] = "success"
                logger.warning(f"Only {num_seqs} sequences found for {output_prefix}.")
            else:
                result["status"] = "success"
            
        except Exception as e:
            result["error"] = str(e)
            logger.error(f"Failed to process {output_prefix}: {e}")
        finally:
            # Clean up temp directory
            if temp_dir.exists():
                import shutil
                shutil.rmtree(temp_dir, ignore_errors=True)
        
        return result
    
    def process_pdb_directory(self, pdb_dir):
        """
        Process all sequences from a PDB directory.
        
        Args:
            pdb_dir: Path to PDB directory containing FASTA file
            
        Returns:
            List of result dictionaries
        """
        pdb_dir = Path(pdb_dir)
        pdb_id = pdb_dir.name
        results = []
        
        # Find FASTA file
        fasta_files = list(pdb_dir.glob("*.fasta")) + list(pdb_dir.glob("*.fa"))
        
        if not fasta_files:
            logger.warning(f"No FASTA file found in {pdb_dir}")
            return results
        
        if len(fasta_files) > 1:
            logger.warning(f"Multiple FASTA files found in {pdb_dir}, using: {fasta_files[0]}")
        
        input_fasta = fasta_files[0]
        logger.info(f"\nProcessing PDB: {pdb_id}")
        logger.info(f"FASTA file: {input_fasta}")
        
        # Read sequences
        sequences = list(SeqIO.parse(input_fasta, "fasta"))
        logger.info(f"Found {len(sequences)} sequences")
        
        # Create output directory
        pdb_output_dir = self.output_dir / pdb_id
        pdb_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Process each sequence
        for idx, seq_record in enumerate(sequences):
            result = self.process_sequence(seq_record, idx, pdb_id, pdb_output_dir)
            results.append(result)
            
            # NCBI requests rate limiting (max 3 requests per second)
            time.sleep(1)
        
        return results
    
    def process_all_directories(self, base_dir):
        """
        Process all PDB directories.
        
        Args:
            base_dir: Base directory containing PDB subdirectories
            
        Returns:
            Dictionary mapping PDB IDs to results
        """
        base_dir = Path(base_dir)
        all_results = {}
        
        pdb_dirs = [d for d in base_dir.iterdir() if d.is_dir()]
        
        if not pdb_dirs:
            logger.error(f"No subdirectories found in {base_dir}")
            return all_results
        
        logger.info(f"Found {len(pdb_dirs)} PDB directories to process")
        
        for pdb_dir in sorted(pdb_dirs):
            results = self.process_pdb_directory(pdb_dir)
            all_results[pdb_dir.name] = results
            
            successful = sum(1 for r in results if r["status"] == "success")
            logger.info(f"\nCompleted {pdb_dir.name}: {successful}/{len(results)} with hits")
        
        return all_results


def main():
    parser = argparse.ArgumentParser(
        description="Generate aligned MSAs using NCBI BLAST + MAFFT"
    )
    parser.add_argument(
        "base_dir",
        help="Base directory containing PDB subdirectories with FASTA files"
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Your email address (required by NCBI)"
    )
    parser.add_argument(
        "--output-dir",
        default="blast_msa_output",
        help="Output directory (default: blast_msa_output)"
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=70.0,
        help="Minimum sequence identity threshold (default: 70)"
    )
    parser.add_argument(
        "--max-hits",
        type=int,
        default=100,
        help="Maximum number of hits per sequence (default: 100)"
    )
    parser.add_argument(
        "--database",
        default="nt",
        help="NCBI database to search (default: nt)"
    )
    
    args = parser.parse_args()
    
    # Check if MAFFT is available
    try:
        subprocess.run(["mafft", "--version"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        logger.error("MAFFT not found! Please install it:")
        logger.error("  conda install -c bioconda mafft")
        logger.error("  OR: sudo apt-get install mafft")
        logger.error("  OR: brew install mafft")
        return
    
    # Initialize generator
    generator = NCBIBlastMSAGenerator(
        output_dir=args.output_dir,
        email=args.email,
        min_identity=args.min_identity,
        max_hits=args.max_hits,
        database=args.database
    )
    
    # Process all directories
    all_results = generator.process_all_directories(args.base_dir)
    
    # Save summary
    summary_file = Path(args.output_dir) / "blast_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("NCBI BLAST + MAFFT MSA GENERATION SUMMARY\n")
        f.write("="*80 + "\n\n")
        
        total_sequences = 0
        total_with_hits = 0
        total_msa_seqs = 0
        
        for pdb_id, results in sorted(all_results.items()):
            f.write(f"\nPDB ID: {pdb_id}\n")
            f.write("-" * 70 + "\n")
            
            for result in results:
                total_sequences += 1
                f.write(f"  Sequence {result['index']}: {result['sequence_id']}\n")
                f.write(f"    Status: {result['status']}\n")
                
                if result['status'] == 'success':
                    total_with_hits += 1
                    total_msa_seqs += result['num_sequences']
                    f.write(f"    MSA size: {result['num_sequences']} sequences\n")
                    f.write(f"    Output: {result['msa_file']}\n")
                elif result['status'] == 'no_hits':
                    f.write(f"    Single sequence MSA created\n")
                    f.write(f"    Output: {result['msa_file']}\n")
                else:
                    f.write(f"    Error: {result['error']}\n")
                f.write("\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write(f"TOTAL RESULTS:\n")
        f.write(f"  Sequences processed: {total_sequences}\n")
        f.write(f"  Sequences with hits: {total_with_hits}\n")
        f.write(f"  Total MSA sequences: {total_msa_seqs}\n")
        if total_with_hits > 0:
            f.write(f"  Average MSA size: {total_msa_seqs/total_with_hits:.1f} sequences\n")
        f.write("="*80 + "\n")
    
    logger.info(f"\nAll processing complete!")
    logger.info(f"Summary saved to {summary_file}")


if __name__ == "__main__":
    main()
    

# HOW TO RUN THIS
# /opt/homebrew/Caskroom/mambaforge/base/bin/python /Users/peneeta/Desktop/grad/courses/fall25/gen_ai_biomed/project/Benchmarking-RhoDesign/ncbi_generate_msa.py ./generated_seqs \
# --email pawojcik@andrew.cmu.edu \
#     --output-dir msa_results \
#     --min-identity 70 \
#     --max-hits 100