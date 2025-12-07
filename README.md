# Benchmarking-RhoDesign

## Peneeta Wojcik

Course project for Generative AI For Biomedicine @ CMU

<hr>

## Description of Directories

- [aptamers](/aptamers/) - reference aptamers in FASTA and PDB formats
- [figs](/figs/) - figures generated for the report
- [generated_seqs](/generated_seqs/) - RhoDesign-generated aptamers at temp=1
- [generated_seqs_lowtemp](/generated_seqs_lowtemp/) - RhoDesign-generated aptamers at temp=1e-5
- [msa_results](/msa_results/) - MSAs for aptamers at temp=1
- [msa_results](/msa_results_lowtemp/) - MSAs for aptamers at temp=1e-5
- [rhofold_pred](/rhofold_pred/) - all successfully-predicted RhoFold+ output files for generated aptamers
- [superimposed](/superimposed/) - example superimposed structures to view in PyMOL

## Description of Files

- [analyze_seqs.ipynb](/analyze_seqs.ipynb) - Sequence analysis for generated structures
- [analyze_structures.ipynb](/analyze_structures.ipynb) - Structure analysis for generated sequences
- [ncbi_generate_msa.py](/ncbi_generate_msa.py) - Automated pipeline for performing BLAST queries and conducting alignment with MAFFT
- [split_fasta_files.py](/split_fasta_files.py) - Splits larger FASTA files into individual FASTA files for easier RhoFold+ input
- [structure_analysis.py](/structure_analysis.py) - Helper functions to calculate RMSD and TM scores and superimpose structures
