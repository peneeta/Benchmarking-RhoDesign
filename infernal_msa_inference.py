#!/usr/bin/env python3
"""
Automated Infernal pipeline for generating MSAs from RNA sequences for RhoFold+
Ended up not getting any hits with this method... Pivoted to NCBI instead (no hits found)
"""

import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio import AlignIO
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class InfernalPipeline:
    """Pipeline for running Infernal to generate MSAs for RNA sequences."""
    
    def __init__(self, rfam_cm_path, output_dir, database_path=None):
        """
        Initialize the Infernal pipeline.
        
        Args:
            rfam_cm_path: Path to Rfam.cm covariance model file
            output_dir: Directory to save output files
            database_path: Optional path to sequence database for homolog searching
        """
        self.rfam_cm_path = rfam_cm_path
        self.output_dir = Path(output_dir)
        self.database_path = database_path
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def extract_individual_sequence(self, seq_record, seq_index):
        """
        Extract and save a single sequence to a temporary file.
        
        Args:
            seq_record: BioPython SeqRecord object
            seq_index: Index of the sequence
            
        Returns:
            Path to the temporary sequence file
        """
        temp_dir = self.output_dir / "temp"
        temp_dir.mkdir(exist_ok=True)
        
        seq_file = temp_dir / f"seq_{seq_index}_{seq_record.id}.fasta"
        SeqIO.write(seq_record, seq_file, "fasta")
        
        return seq_file
    
    def find_best_cm_family(self, seq_file, output_prefix):
        """
        Find the best matching CM family for a sequence using cmscan.
        
        Args:
            seq_file: Path to query sequence file
            output_prefix: Prefix for output files
            
        Returns:
            Tuple of (best_cm_name, best_cm_file_path) or (None, None) if no match
        """
        scan_out = self.output_dir / f"{output_prefix}_scan.txt"
        tbl_out = self.output_dir / f"{output_prefix}_scan.tbl"
        
        cmd = [
            "cmscan",
            "--cut_ga",
            "--rfam",
            "--nohmmonly",
            "-o", str(scan_out),
            "--tblout", str(tbl_out),
            str(self.rfam_cm_path),
            str(seq_file)
        ]
        
        logger.info(f"Running cmscan to find best CM family for {output_prefix}...")
        
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse the table output to get the best hit
            best_cm = None
            with open(tbl_out, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    # First column is the CM name
                    parts = line.split()
                    if parts:
                        best_cm = parts[1]  # Target name (CM family)
                        break
            
            if best_cm:
                logger.info(f"Best matching CM family: {best_cm}")
                
                # Extract this specific CM to a temporary file
                cm_file = self.output_dir / "temp" / f"{best_cm}.cm"
                extract_cmd = [
                    "cmfetch",
                    str(self.rfam_cm_path),
                    best_cm
                ]
                
                with open(cm_file, 'w') as out:
                    subprocess.run(extract_cmd, stdout=out, check=True)
                
                return best_cm, cm_file
            else:
                logger.warning(f"No CM family match found for {output_prefix}")
                return None, None
                
        except subprocess.CalledProcessError as e:
            logger.error(f"cmscan failed for {output_prefix}: {e.stderr}")
            return None, None
    
    def create_alignment(self, seq_file, output_prefix, cm_file=None):
        """
        Create multiple sequence alignment using cmalign.
        
        Args:
            seq_file: Path to sequence file
            output_prefix: Prefix for output files
            cm_file: Optional specific CM file to use (defaults to Rfam.cm)
            
        Returns:
            Path to alignment file
        """
        alignment_out = self.output_dir / f"{output_prefix}_alignment.sto"
        
        if cm_file is None:
            cm_file = self.rfam_cm_path
        
        cmd = [
            "cmalign",
            "-o", str(alignment_out),
            str(cm_file),
            str(seq_file)
        ]
        
        logger.info(f"Running cmalign for {output_prefix}...")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info(f"Alignment created: {alignment_out}")
            return alignment_out
        except subprocess.CalledProcessError as e:
            logger.error(f"cmalign failed for {output_prefix}: {e.stderr}")
            raise
    
    def convert_stockholm_to_a3m(self, sto_file, output_a3m):
        """
        Convert Stockholm format to A3M format for RhoFold.
        
        Args:
            sto_file: Path to Stockholm alignment file
            output_a3m: Path for output A3M file
            
        Returns:
            Path to A3M alignment file
        """
        try:
            # Read Stockholm alignment
            alignment = AlignIO.read(sto_file, "stockholm")
            
            # Write in A3M format (FASTA-like but with specific formatting)
            with open(output_a3m, 'w') as f:
                for record in alignment:
                    # A3M format: header line with '>' followed by sequence
                    # Lowercase letters indicate insertions relative to match states
                    f.write(f">{record.id}\n")
                    # Remove gaps represented as '.' or '-' for insertions
                    seq_str = str(record.seq).replace('.', '').replace('-', '')
                    f.write(f"{seq_str}\n")
            
            logger.info(f"Converted to A3M: {output_a3m}")
            return output_a3m
        except Exception as e:
            logger.error(f"Failed to convert {sto_file} to A3M: {e}")
            raise
    
    def process_sequence(self, seq_record, seq_index, pdb_id, output_dir):
        """
        Process a single sequence through the Infernal pipeline.
        
        Args:
            seq_record: BioPython SeqRecord object
            seq_index: Index of the sequence
            pdb_id: PDB ID for naming
            output_dir: Directory to save output files
            
        Returns:
            Dictionary with results
        """
        output_prefix = f"{pdb_id}_seq_{seq_index}"
        
        logger.info(f"Processing {pdb_id} sequence {seq_index}: {seq_record.id}")
        
        result = {
            "pdb_id": pdb_id,
            "sequence_id": seq_record.id,
            "index": seq_index,
            "status": "failed",
            "cm_family": None,
            "alignment_a3m": None,
            "error": None
        }
        
        try:
            # Extract sequence to temporary file
            seq_file = self.extract_individual_sequence(seq_record, seq_index)
            
            # Find best matching CM family
            cm_family, cm_file = self.find_best_cm_family(seq_file, output_prefix)
            
            if cm_file is None:
                result["error"] = "No matching CM family found in Rfam"
                logger.warning(f"Skipping {pdb_id}_seq_{seq_index}: no CM match")
                return result
            
            result["cm_family"] = cm_family
            
            # Create alignment using the specific CM
            alignment_sto = self.create_alignment(seq_file, output_prefix, cm_file=cm_file)
            
            # Convert to A3M format
            alignment_a3m = output_dir / f"{output_prefix}.a3m"
            self.convert_stockholm_to_a3m(alignment_sto, alignment_a3m)
            result["alignment_a3m"] = str(alignment_a3m)
            
            result["status"] = "success"
            logger.info(f"Successfully processed {pdb_id}_seq_{seq_index} with CM family {cm_family}")
            
        except Exception as e:
            result["error"] = str(e)
            logger.error(f"Failed to process {pdb_id}_seq_{seq_index}: {e}")
        
        return result
    
    def process_pdb_directory(self, pdb_dir):
        """
        Process all sequences from a PDB directory containing a FASTA file.
        
        Args:
            pdb_dir: Path to PDB directory
            
        Returns:
            List of result dictionaries
        """
        pdb_dir = Path(pdb_dir)
        pdb_id = pdb_dir.name
        results = []
        
        # Find FASTA file in the directory
        fasta_files = list(pdb_dir.glob("*.fasta")) + list(pdb_dir.glob("*.fa"))
        
        if not fasta_files:
            logger.warning(f"No FASTA file found in {pdb_dir}")
            return results
        
        if len(fasta_files) > 1:
            logger.warning(f"Multiple FASTA files found in {pdb_dir}, using first one: {fasta_files[0]}")
        
        input_fasta = fasta_files[0]
        logger.info(f"Processing {pdb_id}: {input_fasta}")
        
        # Read sequences
        sequences = list(SeqIO.parse(input_fasta, "fasta"))
        logger.info(f"Found {len(sequences)} sequences in {pdb_id}")
        
        # Create output directory for this PDB
        pdb_output_dir = self.output_dir / pdb_id
        pdb_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Process each sequence
        for idx, seq_record in enumerate(sequences):
            result = self.process_sequence(seq_record, idx, pdb_id, pdb_output_dir)
            results.append(result)
        
        return results
    
    def process_all_pdb_directories(self, base_dir):
        """
        Process all PDB directories in the base directory.
        
        Args:
            base_dir: Base directory containing PDB subdirectories
            
        Returns:
            Dictionary mapping PDB IDs to their results
        """
        base_dir = Path(base_dir)
        all_results = {}
        
        # Find all subdirectories (each should be a PDB ID)
        pdb_dirs = [d for d in base_dir.iterdir() if d.is_dir()]
        
        if not pdb_dirs:
            logger.error(f"No subdirectories found in {base_dir}")
            return all_results
        
        logger.info(f"Found {len(pdb_dirs)} PDB directories to process")
        
        for pdb_dir in sorted(pdb_dirs):
            logger.info(f"\n{'='*60}")
            logger.info(f"Processing PDB directory: {pdb_dir.name}")
            logger.info(f"{'='*60}")
            
            results = self.process_pdb_directory(pdb_dir)
            all_results[pdb_dir.name] = results
            
            successful = sum(1 for r in results if r["status"] == "success")
            logger.info(f"Completed {pdb_dir.name}: {successful}/{len(results)} sequences successful")
        
        return all_results


def main():
    parser = argparse.ArgumentParser(
        description="Generate A3M MSAs for RNA sequences using Infernal for RhoFold structure prediction"
    )
    parser.add_argument(
        "base_dir",
        help="Base directory containing PDB subdirectories, each with a FASTA file"
    )
    parser.add_argument(
        "--rfam-cm",
        required=True,
        help="Path to Rfam.cm covariance model file"
    )
    parser.add_argument(
        "--output-dir",
        default="infernal_output",
        help="Output directory for MSA files (default: infernal_output)"
    )
    
    args = parser.parse_args()
    
    # Initialize pipeline
    pipeline = InfernalPipeline(
        rfam_cm_path=args.rfam_cm,
        output_dir=args.output_dir
    )
    
    # Process all PDB directories
    all_results = pipeline.process_all_pdb_directories(args.base_dir)
    
    # Save comprehensive summary
    summary_file = Path(args.output_dir) / "processing_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("INFERNAL MSA GENERATION SUMMARY\n")
        f.write("="*80 + "\n\n")
        
        total_sequences = 0
        total_successful = 0
        
        for pdb_id, results in sorted(all_results.items()):
            f.write(f"\nPDB ID: {pdb_id}\n")
            f.write("-" * 60 + "\n")
            
            for result in results:
                total_sequences += 1
                f.write(f"  Sequence {result['index']}: {result['sequence_id']}\n")
                f.write(f"    Status: {result['status']}\n")
                
                if result['status'] == 'success':
                    total_successful += 1
                    f.write(f"    CM Family: {result['cm_family']}\n")
                    f.write(f"    Output: {result['alignment_a3m']}\n")
                else:
                    f.write(f"    Error: {result['error']}\n")
                f.write("\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write(f"TOTAL: {total_successful}/{total_sequences} sequences processed successfully\n")
        f.write("="*80 + "\n")
    
    logger.info(f"\n{'='*60}")
    logger.info(f"All processing complete!")
    logger.info(f"Summary saved to {summary_file}")
    logger.info(f"Total: {total_successful}/{total_sequences} sequences successful")
    logger.info(f"{'='*60}")


if __name__ == "__main__":
    main()
    
    
# TO RUN THIS
# /opt/homebrew/Caskroom/mambaforge/base/bin/python /Users/peneeta/Desktop/grad/courses/fall25/gen_ai_biomed/project/Benchmarking-RhoDesign/infernal_msa_inference.py ./generated_seqs --rfam-cm ./Rfam/Rfam.cm