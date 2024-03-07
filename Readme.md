# SCA pipeline for *trans*-AT PKS coevolution analysis
This repo contains the scripts to reproduce the SCA of *trans*-AT PKSs, as published in 
Mabesoone, Leopold-Messer, Minas et al., *Evolution-guided engineering of trans-acyltransferase polyketide synthases*.

`SCA_pipeline.py` is the main script to run the SCA analysis. The SCA scripts as developed by Rama Ranganathan and coworkers are included in the `_sca` folder. [[source](https://reynoldsk.github.io/pySCA/index.html)]

The pipeline is designed to extract domain sequences from Genbank files that are annotated with antiSMASH. 
## Requirements:
- Python 3.6
- [Biopython](https://biopython.org/)
- [NumPy](https://numpy.org/)
- [SciPy](https://www.scipy.org/)
- [Pandas](https://pandas.pydata.org/)
- [Matplotlib](https://matplotlib.org/)

For the multiple sequence alignments:
- [Clustal Omega](http://www.clustal.org/omega/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [MUSCLE](https://www.drive5.com/muscle/)

In our experience, the SCA on sequence alignments of multidomain sequences require > 30 GB of RAM. Hence the pipeline is designed to run on a cluster with a high amount of RAM. The pipeline script is designed to work with  SLURM-based job schedulers.


## Commands
The main script `SCA_pipeline.py` can be run with the following command line arguments:
- `-h` or `--help`: Show the help message and exit.
- `-motif_file`: Path to the file containing the motifs to be extracted. Default is current working directory.
- `-genbank_dir`: Path to directory containing Genbank files to be extracted.
- `-algorithm`: Alignment algorithm. Default is MUSCLE.
- `-leading`: Number of amino acids upstream of the domains to be extracted. Default=0.
- `-trailing`: Number of trailing amino acids of the domains to be extracted. Default=0
- `-in_between`: Number of leading and trailing amino acids in between consecutive domains between extracted. Default=0
- `-minSequenceID`: Minimum sequence identity to be included in the SCA. Default=0.2
- `-maxSequenceID`: Maximum sequence identity to be included in the SCA. Default=0.8
- `-maxGapsFrequency`: Maximum frequency of gaps at positions in the MSA. Default=0.2
- `-minGapsFrequency`: Maximum frequency of gaps at positions in the MSA. Default=0.2
- `-extract`: Switch to turn on or off extraction from the database. Only perform SCA and processing of the SCA results. This requires the sequences to be extracted by a prior run.
- `-any_binding`: Any_binding option of the GenbankExtractor. This treats PCP, ACP and PP-binding domains as similar.
- `-only_analyse`: Switches off sequence extraction, alignment and SCA calculation. Only performs analysis of SCA results.
- `-hybrid_sequence`: Path to fasta file with the sequence of the artificial hybrid. This option ensures inclusion of a custom sequence that contains the query domain motif with a sequence weight of 1.0. This option is not used in the paper.

## Usage:
1. Clone the repository to a cluster running on Linux.
2. Run the script `SCA_pipeline.py` with the desired command line arguments.
3. The script will extract the sequences from the Genbank files, align them, perform the SCA analysis and organize and analyse the output of the SCA.

