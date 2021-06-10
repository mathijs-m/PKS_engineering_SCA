# -*- coding: utf-8 -*-
'''
This script contains the things necessary for the SCA pipeline and shuttles
the data from the GenBank extraction, to the multiple sequence alignment,
pySCA software and finally through the analysis package.

The aim of the SCA pipeline is to automate the analysis of various domain
motifs in trans-AT PKSs with statistical coupling analysis.

Mathijs Mabesoone, ETH Zurich, March 2020.
'''
import GenBankExtractor as gbe
from sys import argv
import sys
from os.path import isdir
from os.path import isfile, join
from os import mkdir
from os import listdir
from os import getcwd
from os import chdir
from os import stat
from datetime import datetime
import subprocess
import traceback
sys.path.append(getcwd()+'/_analysis/')
sys.path.append(getcwd()+'/_sca/')


def parse_command_line_args(args):
    # This function checks if the command line argument is a file containing
    # the domain motifs and returns a list of the domain motifs
    if len(args) < 4:
        print("I didn't receive enough arguments. Exiting.\n", file=logfile)
        return

    motif_file = args[1]

    motif_file = open(args[1]).readlines()

    motifs = []
    for motif in motif_file:
        if motif[-1] == '\n':
            motif = motif[:-1]
        if motif != '':
            motifs.append(list(motif.split(', ')))
    print('INFO: Found motifs:' + str(motifs), file=logfile)

    genbank_dir = argv[2]
    if not isdir(genbank_dir):
        raise ValueError('ERROR: Supplied directory with GenBank files not found.'
              + ' Exiting.\n', file=logfile)
        return
    basepairs = None
    only_annotations = None
    if len(args) == 7:
        basepairs = [int(args[3])*3-1, int(args[4])*3-1, int(args[5])*3-1]
        only_annotations = int(args[6])
    if len(args) == 6:
        basepairs = [int(args[3])*3-1, int(args[4])*3-1]
        only_annotations = int(args[5])

    print( [motifs, genbank_dir, basepairs, only_annotations])

    return [motifs, genbank_dir, basepairs, only_annotations]


def new_folder(folderName):
    folderName  = ''.join(x for x in folderName if x.isalnum() or x in ['_', '-', ',','\\', '/',':'])
    suffix = ''
    while isdir(folderName+str(suffix)):
        if suffix == '':
            suffix = 1
        else:
            suffix += 1

    folderName = folderName + str(suffix)
    mkdir(folderName)
    return folderName


def new_file(fileName):
    if len(fileName.split('.')) > 1:
        extension = '.' + fileName.split('.')[-1]
        fileName = ''.join(fileName.split('.')[0:-1])
    else:
        extension = ''
    fileName  = ''.join(x for x in fileName if x.isalnum() or x in ['_', '-', ',','/'])

    suffix = ''
    while isfile(fileName+str(suffix)+extension):
        if suffix == '':
            suffix = 1
        else:
            suffix += 1

    fileName = fileName + str(suffix) + extension

    open_file = open(fileName, 'w+')
    return fileName, open_file


def MSA_parser_CLUSTAL(motif, folder):
    # This function creates the bash script for the MSA
    motif = ''.join(x for x in motif if x.isalnum() or x in ['_', '-', ','])
    jobname = 'MSA_' + motif
    fname = 'shell_MSA_' + motif + '.sh'
    log_folder = new_folder(join(getcwd(), folder, 'logfolder'))
    shell_folder = new_folder(join(getcwd(),folder,'shellfolder'))
    with open(join(getcwd(), shell_folder, fname), 'w+') as f:
        f.write('#!/bin/sh\n')
        f.write('#$ -cwd           # run in motif directory\n')
        f.write('#$ -S /bin/bash   # interpreting shell for the job\n')
        f.write('#$ -N ' + jobname + '       # name of the job\n')
        f.write('#$ -V             # .bashrc is read in all nodes\n')
        f.write('#$ -pe smp 5      # number of threads to be reserved\n')
        #  If you get errors, try dynamically adapting memory??
        f.write('#$ -l h_vmem=32G  # memory required per thread\n')
        f.write('#$ -e ' + log_folder + '/error_MSA-' + motif + '.log   # error file\n')
        f.write('#S -o output_MSA-' + motif + '    # output file\n\n')
        f.write('module load Clustal-Omega/1.2.4-foss-2018b\n')
        f.write('clustalo -i ' + motif + '_AA.txt --threads=5 --seqtype=Protein\n')
        f.close()

    return join(getcwd(), folder, shell_folder, fname)


def MSA_parser_MAFFT(motif, folder, ep):
    # This function creates the bash script for the MSA
    motif = ''.join(x for x in motif if x.isalnum() or x in ['_', '-', ','])
    jobname = 'MSA_' + motif
    fname = 'shell_MSA_' + motif + '.sh'
    log_folder = new_folder(join(getcwd(), folder, 'logfolder'))
    shell_folder = new_folder(join(getcwd(),folder,'shellfolder'))
    with open(join(getcwd(), shell_folder, fname), 'w+') as f:
        f.write('#!/bin/sh\n')
        f.write('#$ -cwd           # run in motif directory\n')
        f.write('#$ -S /bin/bash   # interpreting shell for the job\n')
        f.write('#$ -N ' + jobname + '       # name of the job\n')
        f.write('#$ -V             # .bashrc is read in all nodes\n')
        f.write('#$ -pe smp 5      # number of threads to be reserved\n')
        #  If you get errors, try dynamically adapting memory??
        f.write('#$ -l h_vmem=32G  # memory required per thread\n')
        f.write('#$ -e ' + log_folder + '/error_MSA-' + motif + '.log   # error file\n')
        f.write('#S -o output_MSA-' + motif + '    # output file\n\n')
        f.write('module load MAFFT/7.305-foss-2018b-with-extensions\n')
        f.write('mafft --auto --thread 5 --ep ' + str(ep) + ' ' +
                motif + '_AA.txt > MSA_' + motif + '.o123 \n')
        f.close()

    return join(getcwd(), folder, shell_folder, fname)


def MSA_parser_MUSCLE(motif, folder):
    # This function creates the bash script for the MSA
    motif = ''.join(x for x in motif if x.isalnum() or x in ['_', '-', ','])
    jobname = 'MSA_' + motif
    fname = 'shell_MSA_' + motif + '.sh'
    log_folder = new_folder(join(getcwd(), folder, 'logfolder'))
    shell_folder = new_folder(join(getcwd(),folder,'shellfolder'))
    with open(join(getcwd(), shell_folder, fname), 'w+') as f:
        f.write('#!/bin/sh\n')
        f.write('#$ -cwd           # run in motif directory\n')
        f.write('#$ -S /bin/bash   # interpreting shell for the job\n')
        f.write('#$ -N ' + jobname + '       # name of the job\n')
        f.write('#$ -V             # .bashrc is read in all nodes\n')
        f.write('#$ -pe smp 5      # number of threads to be reserved\n')
        #  If you get errors, try dynamically adapting memory??
        f.write('#$ -l h_vmem=32G  # memory required per thread\n')
        f.write('#$ -e ' + log_folder + '/error_MSA-' + motif + '.log   # error file\n')
        f.write('#S -o output_MSA-' + motif + '    # output file\n\n')
        f.write('module load MUSCLE/3.8.31-foss-2018b\n')
        f.write('muscle ' + motif + '_AA.txt > MSA_' + motif + '.o123 \n')
        f.close()

    return join(getcwd(), folder, shell_folder, fname)


def SCA_parser(motif, home_folder, argv):
    # This function performs the SCA calculations using the published scripts
    # Note that this function is called from wihtin the motif folder
    motif = ''.join(x for x in motif if x.isalnum() or x in ['_', '-', ','])
    jobname = 'SCA_' + motif
    fname = 'shell_SCA_' + motif + '.sh'
    if type(argv[0]) is list and len(argv[0]) == 4:
        pars = argv[0]
    else:
        pars = [0.2, 0.2, 0.2, 0.8]

    log_folder = join(getcwd(), 'logfolder')
    shell_folder = join(getcwd(), 'shellfolder')
    with open(join(shell_folder, fname), 'w+') as f:
        f.write('#!/bin/sh\n')
        f.write('#$ -cwd          # run in sca directory\n')
        f.write('#$ -S /bin/bash   # interpreting shell for the job\n')
        f.write('#$ -N ' + jobname + '       # name of the job\n')
        f.write('#$ -V             # .bashrc is read in all nodes\n')
        f.write('#$ -pe smp 1      # number of threads to be reserved\n')
        #  If you get errors, try dynamically adapting memory??
        f.write('#$ -l h_vmem=64G  # memory required per thread\n')
        f.write('#$ -e ' + log_folder + '/error_SCA-' + motif + '.log   # error file\n')
        f.write('#S -o output_SCA-' + motif + '    # output file\n\n')
        f.write('module load Python/3.7.0-foss-2018b\n')

        files = listdir(getcwd())
        MSA_files = []
        for file in files:
            if all([file.find('MSA_') == 0, file.find('.o') != -1]):
                MSA_files.append(file)
        filesize = 0
        for file in MSA_files:
            if stat(file).st_size > filesize:
                filesize = stat(file).st_size
                MSA_file = file

        f.write('python ' + home_folder + '/_sca/scaProcessMSA.py -a ' + join(getcwd(), MSA_file) + ' -d ' + join(getcwd()) + ' -p ' + ' '.join([str(i) for i in pars]) + ' \n')
        f.write('python ' + home_folder + '/_sca/scaCore.py -i ' + join(getcwd(), 'MSA_' + motif + '.db') + ' -o ' + join(getcwd(), motif + '_sca.db') + ' \n')
        f.write('python ' + home_folder + '/_sca/scaSectorID.py -i ' + join(getcwd(), motif + '_sca.db') + ' -o ' + join(getcwd(), motif + '_sectors.db') + ' \n')
        f.close()

    return join(shell_folder, fname)


def SCAanalysis_parser(motif, home_folder):
    # This function performs the SCA calculations using the published scripts
    motif = ''.join(x for x in motif if x.isalnum() or x in ['_', '-', ','])
    jobname = 'SCAanl_' + motif
    shell_folder = join(getcwd(), 'shellfolder')
    fname = 'shell_SCAanalysis_' + motif + '.sh'
    with open(join(shell_folder, fname), 'w+') as f:
        f.write('#!/bin/sh\n')
        f.write('#$ -cwd          # run in sca directory\n')
        f.write('#$ -S /bin/bash   # interpreting shell for the job\n')
        f.write('#$ -N ' + jobname + '       # name of the job\n')
        f.write('#$ -V             # .bashrc is read in all nodes\n')
        f.write('#$ -pe smp 1      # number of threads to be reserved\n')
        #  If you get errors, try dynamically adapting memory??
        f.write('#$ -l h_vmem=128G  # memory required per thread\n')
        f.write('#$ -e error_SCAanalysis-' + motif + '.log   # error file\n')
        f.write('#S -o output_SCAanalysis-' + motif + '    # output file\n\n')
        f.write('module load Python/3.7.0-foss-2018b\n')

        f.write('python ' + home_folder + '/_analysis/KSanalysis.py ' + motif + ' ' + home_folder + '/_sca \n')
        f.close()

    return join(shell_folder, fname)


def main(max_freq_gaps, extract):
    # Need as arguments:
    # 1: motif file - 2: genbank directory - 3: leading basepairs - 4: trailing
    # basepairs

    minSeqID = 0.2
    maxSeqID = 0.8

    if extract:
        print('INFO: Starting...')
        cur_time = datetime.now()
        cur_time = (str(cur_time.year) + '_' + str(cur_time.month) + '_' +
                    str(cur_time.day) + '-' + str(cur_time.hour) + 'h' +
                    str(cur_time.minute))

        global logfile
        log_folder = 'log_files'
        try:
            mkdir(log_folder)
        except FileExistsError:
            pass
        [logfilename, logfile] = new_file(log_folder + '/logfile_' +
                                          cur_time+'.txt')
        sys.stdout = logfile

    motif_file = argv[1]
    genbank_dir = argv[2]
    algorithm = argv[3]
    leader = argv[4]
    follower = argv[5]
    intermediate = argv[6]
    only_annotations = argv[7]
    parse_input = ['', motif_file, genbank_dir, leader, follower, intermediate,
                   only_annotations]

    [motifs, genbank_dir, basepairs, only_annotations] = parse_command_line_args(parse_input)

    db_props = gbe.load_database_properties_api(genbank_dir)

    organisms = db_props['organisms']
    BGCs = db_props['BGCs']
    specificities = db_props['specificities']

    home_folder = getcwd()

    for motif in motifs:
        query_motif = motif
        motif = '_'.join(motif)
        print('\nINFO: Started analyzing ' + motif + '\n')
        motif_folder = new_folder(motif)

        # ==================================================
        # Perform the extractions and save data
        print('\nINFO: Created folder: ' + motif_folder, file=logfile)

        numSeqs = gbe.main_extractor(genbank_dir, query_motif, basepairs, BGCs,
                                     organisms, specificities, db_props,
                                     motif, motif_folder, only_annotations)[0]
        print('Extraction of ' + str(motif) + ' successful.', file=logfile)
        print('Exctracted ' + str(numSeqs) + ' sequences.', file=logfile)

        print('INFO: Parsing MSA shell script')
        try:
            if algorithm.lower() == 'clustal':
                MSA_file = MSA_parser_CLUSTAL(motif, motif_folder)
            elif algorithm.lower() == 'muscle':
                MSA_file = MSA_parser_MUSCLE(motif, motif_folder)
            elif 'mafft' in algorithm.lower():
                ep = float(algorithm.split('_')[1])
                MSA_file = MSA_parser_MAFFT(motif, motif_folder, ep)

        except Exception as e:
            print(e)
            print('ERROR: MSA parser error')
            traceback.print_exc()
        # Move to the motif folder and perform the MSA
        chdir('/'.join([home_folder, motif_folder]))
        f = open('/'.join([home_folder, motif_folder,
                           'MaxFreqGapsList-' + str(max_freq_gaps)]), 'w+')
        f.write(str(max_freq_gaps))
        f.close()

        if not sys.platform == 'win32':
            print('INFO: Starting MSA')
            subprocess.run('qsub -sync y ' + MSA_file, shell=True)
            print('INFO: Proceeding to SCA analysis')

        # ==================================================
        # Now do the SCA
        print('INFO: Parsing SCA shell script')
        try:
            pars = [max_freq_gaps, max_freq_gaps, minSeqID, maxSeqID]
            SCA_file = SCA_parser(motif, home_folder, pars)
        except Exception:
            print('ERROR: SCA parser error')
            traceback.print_exc()

        if not sys.platform == 'win32':
            print('INFO: Starting SCA calculations')
            subprocess.run('qsub -sync y ' + SCA_file, shell=True)
        print('INFO: SCA completed')
        print('\nINFO: Proceeding to SCA processing')

        # ==================================================
        # Now analyze the SCA
        try:
            motif = ''.join(x for x in motif if x.isalnum() or x in ['_', '-', ','])
            SCAanalysis_file = SCAanalysis_parser(motif, home_folder)
            subprocess.run('qsub -sync y ' + SCAanalysis_file, shell=True)
            print('INFO: Finished processing ' + motif + '\n')
        except Exception:
            traceback.print_exc()
            print('FAILED: Failed to perform the SCA processing for ' + motif + '\n')

        chdir(home_folder)
        print('MOVING ON TO NEXT MOTIF\n')
        print('====================================')


if __name__ == '__main__':
    max_freq_gaps_list = [0.2]
    extract = 1
    for max_freq_gaps in max_freq_gaps_list:
        main(max_freq_gaps, extract)
        extract = 0
