import os
from pathlib import Path

def SplitFASTA(base_directory):
    """
    Splits fasta files in subdirectories into individual sequence files.
    
    Args:
        base_directory: Path to the directory containing subdirectories with fasta files
    """
    base_path = Path(base_directory)
    
    # Iterate through each subdirectory
    for subdir in base_path.iterdir():
        if not subdir.is_dir():
            continue
            
        # Find fasta files in the subdirectory
        fasta_files = list(subdir.glob("*.fasta")) + list(subdir.glob("*.fa"))
        
        if not fasta_files:
            print(f"No fasta files found in {subdir.name}")
            continue
        
        # Process the first fasta file found
        fasta_file = fasta_files[0]
        print(f"Processing {fasta_file.name} in {subdir.name}")
        
        # Create output folder within the subdirectory
        output_folder = subdir / "split_sequences"
        output_folder.mkdir(exist_ok=True)
        
        # Read and split the fasta file
        with open(fasta_file, 'r') as f:
            sequences = []
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(''.join(current_seq))
                        current_seq = []
                    current_seq.append(line + '\n')
                else:
                    current_seq.append(line + '\n')
            
            # Add the last sequence
            if current_seq:
                sequences.append(''.join(current_seq))
        
        # Write individual sequence files
        for i, seq in enumerate(sequences):
            output_file = output_folder / f"{subdir.name}_seq_{i}.fasta"
            with open(output_file, 'w') as f:
                f.write(seq)
        
        print(f"Created {len(sequences)} files in {output_folder}")

# for generated_seqs (temp=1)
SplitFASTA("./generated_seqs")