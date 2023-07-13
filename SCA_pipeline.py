# -*- coding: utf-8 -*-
'''
This script contains the things necessary for the SCA pipeline and shuttles
the data from the GenBank extraction, to the multiple sequence alignment,
pySCA software and finally through the analysis package.

The aim of the SCA pipeline is to automate the analysis of various domain
motifs in trans-AT PKSs with statistical coupling analysis.

Mathijs Mabesoone, ETH Zurich, March 2020.
'''
import argparse
import sys
import os
from datetime import datetime
import subprocess

sys.path.append(os.getcwd()+'/_analysis/')
sys.path.append(os.getcwd()+'/_sca/')
import sequence_extractor

def parse_command_line_args(motif_file, leader, follower, intermediate):
    # This function checks if the command line argument is a file containing
    # the domain motifs and returns a list of the domain motifs

    motifs = [motif.replace('\n','').split(', ') for motif in open(motif_file, 'r').readlines() if motif != '']
    sys.stdout.write(f"INFO: Found motifs: {motifs}\n")

    basepairs = None
    if intermediate > 0:
        basepairs = [leader*3-1, follower*3-1, intermediate*3-1]
    else:
        basepairs = [leader*3-1, follower*3-1]

    return motifs, basepairs


def new_folder(folderName):
    folderName  = ''.join(x for x in folderName if x.isalnum() or x in ['_', '-', ',','\\', '/',':'])
    suffix = ''
    while os.path.isdir(folderName+str(suffix)):
        if suffix == '':
            suffix = 1
        else:
            suffix += 1

    folderName = folderName + str(suffix)
    os.mkdir(folderName)
    return folderName


def new_file(fileName):
    if len(fileName.split('.')) > 1:
        extension = '.' + fileName.split('.')[-1]
        fileName = ''.join(fileName.split('.')[0:-1])
    else:
        extension = ''
    fileName  = ''.join(x for x in fileName if x.isalnum() or x in ['_', '-', ',','/'])

    suffix = ''
    while os.path.isfile(fileName+str(suffix)+extension):
        if suffix == '':
            suffix = 1
        else:
            suffix += 1

    fileName = fileName + str(suffix) + extension

    open_file = open(fileName, 'w+')
    return fileName, open_file

def add_hybrid_sequence_to_AA_txt(motif_folder, hybrid_path):
    'Function to add the hybrid fasta to the extracted amino acid sequences to feed into the MSA.'
    AA_file = [os.path.join(motif_folder, file) for file in os.listdir(motif_folder) if '_AA.fasta' in file]
    if len(AA_file)>1:
        raise IndexError('AA_file too long')
    with open(AA_file[0], 'a') as destination:
        with open(hybrid_path, 'r') as source:
            sequence = source.read()
            destination.write(sequence)

def MSA_parser(motif, home_folder, motif_folder, algorithm):
    # This function creates the bash script for the MSA
    motif = ''.join(x for x in motif if x.isalnum() or x in ['_', '-', ',', ' ']).replace(' ','_')
    fname = 'shell_MSA_' + motif + '.sh'
    log_folder = new_folder(os.path.join(motif_folder, 'logfolder'))
    shell_folder = new_folder(os.path.join(motif_folder,'shellfolder'))
    fasta = os.path.join(motif_folder, motif+'_AA.fasta')
    if not os.path.isfile(fasta):
        fasta = [os.path.join(home_folder, file) for file in os.listdir(home_folder) if '.fasta' in file][0]
        sys.stdout.write(f"Analyzing {fasta}\n")
    with open(os.path.join(shell_folder, fname), 'w+') as f:
        f.write("#!/bin/sh\n")
        f.write(f"#$ -wd {motif_folder}           # run in motif directory\n")
        f.write("#$ -S /bin/bash   # interpreting shell for the job\n")
        f.write(f"#$ -N MSA_{motif}       # name of the job\n")
        f.write("#$ -V             # .bashrc is read in all nodes\n")
        f.write("#$ -pe smp 5      # number of threads to be reserved\n")
        #  If you get errors, try dynamically adapting memory??
        f.write("#$ -l h_vmem=32G  # memory required per thread\n")
        f.write(f"#$ -e {log_folder}/error_MSA-{motif}.log   # error file\n")
        f.write("#S -o output_MSA-' + motif + '    # output file\n\n")
        if algorithm.lower() == 'clustal':
            f.write("module load Clustal-Omega/1.2.4-foss-2018b\n")
            f.write(f"clustalo -i {fasta} --threads=5 --seqtype=Protein -o {os.path.join(motif_folder, motif)}.o123\n")
        elif algorithm.lower() == 'muscle':
            f.write('module load MUSCLE/3.8.31-foss-2018b\n')
            f.write(f'muscle -in {fasta} -out {os.path.join(motif_folder, "MSA_" + motif)}.o123 \n')
        elif 'mafft' in algorithm.lower():
            ep = float(algorithm.split('_')[1])
            f.write('module load MAFFT/7.305-foss-2018b-with-extensions\n')
            f.write(f"mafft --auto --thread 5 --ep {ep} {fasta} > {os.path.join(motif_folder, motif)}.o123 \n")
        f.close()

    return os.path.join(shell_folder, fname)


def SCA_parser(motif, home_folder, motif_folder, sca_parameters):
    # This function performs the SCA calculations using the published scripts
    # Note that this function is called from wihtin the motif folder
    motif = ''.join(x for x in motif if x.isalnum() or x in ['_', '-', ',', ' ']).replace(' ','_')
    if type(sca_parameters) != list:
        sca_parameters = [0.2, 0.2, 0.2, 0.8]
    elif len(sca_parameters)==0:
        sca_parameters = [0.2, 0.2, 0.2, 0.8]

    log_folder = os.path.join(motif_folder, 'logfolder')
    shell_folder = os.path.join(motif_folder, 'shellfolder')
    with open(os.path.join(shell_folder, f'shell_SCA_{motif}.sh'), 'w+') as f:
        f.write('#!/bin/sh\n')
        f.write(f'#$ -wd {motif_folder}          # run in sca directory\n')
        f.write('#$ -S /bin/bash   # interpreting shell for the job\n')
        f.write(f'#$ -N SCA_{motif}       # name of the job\n')
        f.write('#$ -V             # .bashrc is read in all nodes\n')
        f.write('#$ -pe smp 1      # number of threads to be reserved\n')
        #  If you get errors, try dynamically adapting memory??
        f.write('#$ -l h_vmem=64G  # memory required per thread\n')
        f.write(f'#$ -e {log_folder}/error_SCA-{motif}.log   # error file\n')
        f.write(f'#S -o output_SCA-{motif}    # output file\n\n')
        f.write('module load Python/3.6.6-foss-2018b\n')

        MSA_files = [os.path.join(motif_folder, file) for file in os.listdir(motif_folder) if '.o123' in file]
        filesize = 0
        for file in MSA_files:
            if os.stat(file).st_size > filesize:
                filesize = os.stat(file).st_size
                MSA_file = file

        f.write(f"python {home_folder}/_sca/scaProcessMSA.py -a {os.path.join(motif_folder, MSA_file)} -d {motif_folder} -p {' '.join([str(i) for i in sca_parameters])}\n")
        f.write(f"python {home_folder}/_sca/scaCore.py -i {os.path.join(motif_folder, motif + '.db')} -o {os.path.join(motif_folder, motif + '_sca.db')}\n")
        f.write(f"python {home_folder}/_sca/scaSectorID.py -i {os.path.join(motif_folder, motif + '_sca.db')} -o {os.path.join(motif_folder, motif + '_sectors.db')}\n")
        f.close()

    return os.path.join(shell_folder, f'shell_SCA_{motif}.sh')


def SCAanalysis_parser(motif, home_folder, motif_folder):
    # This function performs the SCA calculations using the published scripts
    motif = ''.join(x for x in motif if x.isalnum() or x in ['_', '-', ',', ' ']).replace(' ','_')
    shell_folder = os.path.join(motif_folder, 'shellfolder')
    fname = 'shell_SCAanalysis_' + motif + '.sh'
    with open(os.path.join(shell_folder, fname), 'w+') as f:
        f.write('#!/bin/sh\n')
        f.write(f"#$ -wd {motif_folder}          # run in sca directory\n")
        f.write('#$ -S /bin/bash   # interpreting shell for the job\n')
        f.write(f'#$ -N SCAanl_{motif}       # name of the job\n')
        f.write('#$ -V             # .bashrc is read in all nodes\n')
        f.write('#$ -pe smp 1      # number of threads to be reserved\n')
        #  If you get errors, try dynamically adapting memory??
        f.write('#$ -l h_vmem=128G  # memory required per thread\n')
        f.write(f'#$ -e error_SCAanalysis-{motif}.log   # error file\n')
        f.write(f'#S -o output_SCAanalysis-{motif}    # output file\n\n')
        f.write('module load Python/3.6.6-foss-2018b\n')
        f.write(f'python {home_folder}/_analysis/KSanalysis.py -motif_dir {motif_folder} -pysca_dir {home_folder}/_sca \n')
        f.close()

    return os.path.join(shell_folder, fname)


def main():
    # Need as arguments:
    # 1: motif file - 2: genbank directory - 3: leading basepairs - 4: trailing
    # basepairs

    parser = argparse.ArgumentParser(description='Pipeline script to extract domain sequences from a database, allign them and perform and analyze SCA.')
    parser.add_argument('-motif_file',help='Path to the file containing the motifs to be extracted. Default is current working directory.', type=str, default=os.getcwd())
    parser.add_argument('-genbank_dir', help='Path to directory containing Genbank files to be extracted.', type=str)
    parser.add_argument('-algorithm', help='Alignment algorithm. Default is MUSCLE', type=str, default='MUSCLE')
    parser.add_argument('-leading', help='Number of amino acids upstream of the domains to be extracted. Default=0.', type=int, default=0)
    parser.add_argument('-trailing', help='Number of trailing amino acids of the domains to be extracted. Default=0', type=int, default=0)
    parser.add_argument('-in_between', help='Number of leading and trailing amino acids in between consecutive domains between extracted. Default=0', type=int, default=0)
    parser.add_argument('-minSequenceID', help='Minimum sequence identity to be included in the SCA. Default=0.2', type=float, default=0.2)
    parser.add_argument('-maxSequenceID', help='Maximum sequence identity to be included in the SCA. Default=0.8', type=float, default=0.8)
    parser.add_argument('-maxGapsFrequency',help='Maximum frequency of gaps at positions in the MSA. Default=0.2', type=float, default=0.2)
    parser.add_argument('-minGapsFrequency',help='Maximum frequency of gaps at positions in the MSA. Default=0.2', type=float, default=0.2)
    parser.add_argument('-extract', help='Switch to turn on or off extraction from the database. Values: 0 or 1. Default=1', type=bool, default=True)
    parser.add_argument('-any_binding', help='Any_binding option of the GenbankExtractor.')
    parser.add_argument('-only_analyse', help='Switches off sequence extraction, alignment and SCA calculation. Only performs analysis of SCA results. Values: 1, 0. Default=0.', type=bool, default=False)
    parser.add_argument('-hybrid_sequence', help='Path to fasta file with the sequence of the artificial hybrid.', type=str)
    
    args = parser.parse_args()

    if args.only_analyse:
        args.extract=False
    try:
        os.mkdir('log_files')
    except FileExistsError:
        pass

    stdout_old = sys.stdout
    if not sys.platform == 'win32':
        sys.stdout = new_file('log_files' + '/logfile_' + datetime.now().strftime("%Y_%m_%d-%Hh%m")+'.txt')[1]

    # Load some properties
    with open(args.motif_file, 'r') as f:
        motifs = f.readlines()
    motifs = [motif.replace('\n', '').split(', ') for motif in motifs]

    # Save the home folder
    home_folder = os.path.dirname(args.motif_file)
    if not os.path.isfile(args.motif_file):
        home_folder = args.motif_file

    # Loop over the motifs
    for motif in motifs:
        query_motif = motif
        if type(motif)==list:
            motif = '_'.join(motif)
        if not args.only_analyse:
            sys.stdout.write(f"INFO: Started analyzing {motif}\n")
            motif_folder = new_folder(os.path.join(home_folder,os.path.splitext(motif)[0]))

            # ==================================================
            # Perform the extractions and save data
            sys.stdout.write(f"INFO: Created folder: {motif_folder}\n")
            stdout_old.write(f"INFO: Created folder: {motif_folder}\n")
            if args.extract:
                sys.stdout.write(f'sequence_extractor.extract_domain_sequences({query_motif}, {args.genbank_dir}, {motif_folder}, any_binding = {args.any_binding})')
                numSeqs = sequence_extractor.extract_domain_sequences(query_motif, args.genbank_dir, motif_folder, any_binding = args.any_binding)[0]
                #numSeqs = gbe.main_extractor(args.genbank_dir, args.any_binding, query_motif, fname=motif.replace(' ',''),
                #                         dest_folder=os.path.join(home_folder,motif.replace(' ','')))[0]
                print('extracted')
                #add_hybrid_sequence_to_AA_txt(os.path.join(motif_folder), args.hybrid_sequence)
                sys.stdout.write(f"Extraction of {motif} successful.\n")
                sys.stdout.write('Exctracted ' + str(numSeqs) + ' sequences.\n')

            print('INFO: Parsing MSA shell script')
            try:
                
                MSA_file = MSA_parser(motif, home_folder, motif_folder, args.algorithm)
            except Exception as e:
                sys.stdout.write(e)
                sys.stdout.write('ERROR: MSA parser error\n')
                stdout_old.write('ERROR: MSA parser error\n')

            # Move to the motif folder and perform the MSA
            f = open(os.path.join(motif_folder, f'MaxFreqGapsList-{args.maxGapsFrequency}'), 'w+')
            f.write(str(args.maxGapsFrequency))
            f.close()

            if not sys.platform == 'win32':
                sys.stdout.write('INFO: Starting MSA\n')
                stdout_old.write('INFO: Starting MSA\n')
                subprocess.run('qsub -sync y ' + MSA_file, shell=True)
                sys.stdout.write('INFO: Proceeding to SCA analysis\n')
                stdout_old.write('INFO: Proceeding to SCA analysis\n')

            # ==================================================
            # Do the SCA
            stdout_old.write('INFO: Parsing SCA shell script')
            try:
                SCA_file = SCA_parser(motif, home_folder, motif_folder, [args.minGapsFrequency, args.maxGapsFrequency, args.minSequenceID, args.maxSequenceID])
            except Exception as e:
                sys.stdout.write('ERROR: SCA parser error\n')
                sys.stdout.write(f"{e}\n")

            if not sys.platform == 'win32':
                sys.stdout.write('INFO: Starting SCA calculations\n')
                subprocess.run('qsub -sync y ' + SCA_file, shell=True)
            sys.stdout.write('INFO: SCA completed\n')
            sys.stdout.write('INFO: Proceeding to SCA processing\n')

        # ==================================================
        # Analyze the SCA
        try:
            motif = ''.join(x for x in motif if x.isalnum() or x in ['_', '-', ','])
            SCAanalysis_file = SCAanalysis_parser(motif, home_folder, motif_folder)
            subprocess.run('qsub -sync y ' + SCAanalysis_file, shell=True)
            sys.stdout.write('INFO: Finished processing ' + motif + '\n')
        except Exception as e:
            sys.stdout.write(e)
            sys.stdout.write('FAILED: Failed to perform the SCA processing for ' + motif + '\n')

        sys.stdout.write('MOVING ON TO NEXT MOTIF\n')
        sys.stdout.write('====================================\n')


if __name__ == '__main__':
    main()
