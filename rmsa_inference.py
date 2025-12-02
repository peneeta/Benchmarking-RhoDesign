#!/usr/bin/env python3
"""
Batch processor for running rMSA on multiple RNA sequences organized by PDB ID.
Generates A3M MSA files for each sequence suitable for RhoFold structure prediction.
"""

import subprocess
import os
from pathlib import Path
from Bio import SeqIO
import argparse
import logging
import shutil

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class rMSABatchProcessor:
    """Process multiple RNA sequences through rMSA pipeline."""
    
    def __init__(self, rmsa_path, output_dir, cpu=4, fast=True, timeout=3600):
        """
        Initialize the rMSA batch processor.
        
        Args:
            rmsa_path: Path to rMSA.pl script
            output_dir: Directory to save output files
            cpu: Number of CPU threads to use
            fast: Use fast mode (default True)
            timeout: Timeout in seconds per sequence (default 3600)
        """
        self.rmsa_path = Path(rmsa_path)
        self.output_dir = Path(output_dir)
        self.cpu = cpu
        self.fast = 1 if fast else 0
        self.timeout = timeout
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Verify rMSA exists
        if not self.rmsa_path.exists():
            raise FileNotFoundError(f"rMSA.pl not found at {self.rmsa_path}")
    
    def extract_individual_sequence(self, seq_record, temp_dir, seq_index):
        """
        Extract a single sequence to a temporary FASTA file.
        
        Args:
            seq_record: BioPython SeqRecord object
            temp_dir: Temporary directory path
            seq_index: Index of the sequence
            
        Returns:
            Path to the temporary sequence file
        """
        temp_dir.mkdir(exist_ok=True, parents=True)
        seq_file = temp_dir / f"temp_seq_{seq_index}.fasta"
        SeqIO.write(seq_record, seq_file, "fasta")
        return seq_file
    
    def run_rmsa(self, input_fasta, output_prefix, temp_dir):
        """
        Run rMSA on a single sequence.
        
        Args:
            input_fasta: Path to input FASTA file
            output_prefix: Prefix for output files
            temp_dir: Temporary directory for intermediate files
            
        Returns:
            Path to output A3M file or None if failed
        """
        # rMSA outputs to the same directory as input with .afa extension
        # We'll run it in temp dir then move the output
        temp_input = temp_dir / f"{output_prefix}_input.fasta"
        shutil.copy(input_fasta, temp_input)
        
        cmd = [
            str(self.rmsa_path),
            str(temp_input),
            f"-cpu={self.cpu}",
            f"-fast={self.fast}",
            f"-timeout={self.timeout}",
            f"-tmpdir={temp_dir / 'rmsa_tmp'}"
        ]
        
        logger.info(f"Running rMSA for {output_prefix}...")
        logger.info(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout + 300  # Add buffer to timeout
            )
            
            # rMSA outputs <input>.afa
            output_afa = temp_input.with_suffix('.afa')
            
            if output_afa.exists():
                logger.info(f"rMSA completed successfully for {output_prefix}")
                return output_afa
            else:
                logger.error(f"rMSA did not produce output for {output_prefix}")
                logger.error(f"stdout: {result.stdout}")
                logger.error(f"stderr: {result.stderr}")
                return None
                
        except subprocess.TimeoutExpired:
            logger.error(f"rMSA timed out for {output_prefix}")
            return None
        except Exception as e:
            logger.error(f"rMSA failed for {output_prefix}: {e}")
            return None
    
    def convert_afa_to_a3m(self, afa_file, a3m_file):
        """
        Convert .afa (A2M) format to .a3m format.
        
        Args:
            afa_file: Path to input AFA file
            a3m_file: Path to output A3M file
            
        Note: rMSA outputs in A2M format which is essentially the same as A3M.
        We just rename and ensure proper format.
        """
        try:
            # Read and write to ensure proper formatting
            with open(afa_file, 'r') as f_in, open(a3m_file, 'w') as f_out:
                for line in f_in:
                    # A3M format: similar to A2M/AFA
                    # Lowercase = insertions, uppercase = match states
                    f_out.write(line)
            
            logger.info(f"Converted to A3M: {a3m_file}")
            return a3m_file
        except Exception as e:
            logger.error(f"Failed to convert {afa_file} to A3M: {e}")
            raise
    
    def process_sequence(self, seq_record, seq_index, pdb_id, pdb_output_dir):
        """
        Process a single sequence through rMSA.
        
        Args:
            seq_record: BioPython SeqRecord object
            seq_index: Index of the sequence
            pdb_id: PDB ID for naming
            pdb_output_dir: Output directory for this PDB
            
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
            "msa_file": None,
            "num_sequences": 0,
            "error": None
        }
        
        # Create temp directory for this sequence
        temp_dir = self.output_dir / "temp" / pdb_id / f"seq_{seq_index}"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # Extract sequence to temp file
            seq_file = self.extract_individual_sequence(seq_record, temp_dir, seq_index)
            
            # Run rMSA
            afa_file = self.run_rmsa(seq_file, output_prefix, temp_dir)
            
            if afa_file is None:
                result["error"] = "rMSA did not produce output"
                return result
            
            # Convert to A3M and save to final location
            final_a3m = pdb_output_dir / f"{output_prefix}.a3m"
            self.convert_afa_to_a3m(afa_file, final_a3m)
            
            # Count sequences in MSA
            num_seqs = sum(1 for _ in SeqIO.parse(final_a3m, "fasta"))
            
            result["msa_file"] = str(final_a3m)
            result["num_sequences"] = num_seqs
            result["status"] = "success"
            
            logger.info(f"Successfully processed {output_prefix}: {num_seqs} sequences in MSA")
            
            # Clean up temp directory
            shutil.rmtree(temp_dir, ignore_errors=True)
            
        except Exception as e:
            result["error"] = str(e)
            logger.error(f"Failed to process {output_prefix}: {e}")
        
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
    
    def process_all_directories(self, base_dir):
        """
        Process all PDB directories in the base directory.
        
        Args:
            base_dir: Base directory containing PDB subdirectories
            
        Returns:
            Dictionary mapping PDB IDs to their results
        """
        base_dir = Path(base_dir)
        all_results = {}
        
        # Find all subdirectories
        pdb_dirs = [d for d in base_dir.iterdir() if d.is_dir()]
        
        if not pdb_dirs:
            logger.error(f"No subdirectories found in {base_dir}")
            return all_results
        
        logger.info(f"Found {len(pdb_dirs)} PDB directories to process")
        
        for pdb_dir in sorted(pdb_dirs):
            logger.info(f"\n{'='*70}")
            logger.info(f"Processing PDB directory: {pdb_dir.name}")
            logger.info(f"{'='*70}")
            
            results = self.process_pdb_directory(pdb_dir)
            all_results[pdb_dir.name] = results
            
            successful = sum(1 for r in results if r["status"] == "success")
            total_seqs = sum(r["num_sequences"] for r in results if r["status"] == "success")
            
            logger.info(f"Completed {pdb_dir.name}: {successful}/{len(results)} sequences successful")
            logger.info(f"Total MSA sequences generated: {total_seqs}")
        
        return all_results


def main():
    parser = argparse.ArgumentParser(
        description="Generate MSAs for RNA sequences using rMSA for RhoFold structure prediction"
    )
    parser.add_argument(
        "base_dir",
        help="Base directory containing PDB subdirectories with FASTA files"
    )
    parser.add_argument(
        "--rmsa-path",
        default="./rMSA/rMSA.pl",
        help="Path to rMSA.pl script (default: ./rMSA/rMSA.pl)"
    )
    parser.add_argument(
        "--output-dir",
        default="rmsa_output",
        help="Output directory for MSA files (default: rmsa_output)"
    )
    parser.add_argument(
        "--cpu",
        type=int,
        default=4,
        help="Number of CPU threads (default: 4)"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=3600,
        help="Timeout per sequence in seconds (default: 3600)"
    )
    parser.add_argument(
        "--no-fast",
        action="store_true",
        help="Disable fast mode (more thorough but slower)"
    )
    
    args = parser.parse_args()
    
    # Initialize processor
    processor = rMSABatchProcessor(
        rmsa_path=args.rmsa_path,
        output_dir=args.output_dir,
        cpu=args.cpu,
        fast=not args.no_fast,
        timeout=args.timeout
    )
    
    # Process all directories
    all_results = processor.process_all_directories(args.base_dir)
    
    # Save comprehensive summary
    summary_file = Path(args.output_dir) / "rmsa_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("rMSA BATCH PROCESSING SUMMARY\n")
        f.write("="*80 + "\n\n")
        
        total_sequences = 0
        total_successful = 0
        total_msa_seqs = 0
        
        for pdb_id, results in sorted(all_results.items()):
            f.write(f"\nPDB ID: {pdb_id}\n")
            f.write("-" * 70 + "\n")
            
            for result in results:
                total_sequences += 1
                f.write(f"  Sequence {result['index']}: {result['sequence_id']}\n")
                f.write(f"    Status: {result['status']}\n")
                
                if result['status'] == 'success':
                    total_successful += 1
                    total_msa_seqs += result['num_sequences']
                    f.write(f"    MSA size: {result['num_sequences']} sequences\n")
                    f.write(f"    Output: {result['msa_file']}\n")
                else:
                    f.write(f"    Error: {result['error']}\n")
                f.write("\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write(f"TOTAL RESULTS:\n")
        f.write(f"  Processed: {total_successful}/{total_sequences} sequences successful\n")
        f.write(f"  Total MSA sequences generated: {total_msa_seqs}\n")
        f.write(f"  Average MSA size: {total_msa_seqs/total_successful:.1f} sequences\n" if total_successful > 0 else "")
        f.write("="*80 + "\n")
    
    logger.info(f"\n{'='*70}")
    logger.info(f"All processing complete!")
    logger.info(f"Summary saved to {summary_file}")
    logger.info(f"Results: {total_successful}/{total_sequences} sequences successful")
    logger.info(f"{'='*70}")


if __name__ == "__main__":
    main()