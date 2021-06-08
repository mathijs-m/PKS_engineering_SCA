# -*- coding: utf-8 -*-
"""
Script to analyze the SCA of the KS alignment as done by Allison

Mathijs Mabesoone, 01.02.2021
"""

import os
import shutil
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import scipy.cluster.hierarchy as sch
from scipy.stats import scoreatpercentile
import pickle
import traceback
import sys
sys.path.append('../_analysis')
sys.path.append('../_sca')
import KStools
from pysca import sca


def newFolder(folderName):
    if not os.path.isdir(folderName):
        os.mkdir(folderName)
    return folderName


def construct_phylogeny_and_annotations(enzymelist):
    annot = dict()
    spec = dict()

    for num, entry in enumerate(enzymelist):
        s1 = entry.split('|')
        bgc_name = s1[0].strip()
        gene = s1[1].strip()
        organism = s1[3].strip()
        print(bgc_name)
        print(gene)
        print(organism)

        annot = ''
        if (bgc_name in spec.keys()):
            spec[bgc_name].append(gene)
        else:
            spec[bgc_name] = [gene]
    return [annot, spec]


def obtain_reordered_sector_data(sector_groups, IC_list):
    # This function obtains the data for the reclustered sectors and returns
    # these reorederded sectors

    sectors = list()

    for n, reordered_sector in enumerate(sector_groups):
        new_sector = KStools.Unit()
        all_residues = list()
        all_vector_values = list()
        for subsector in reordered_sector:
            all_residues = all_residues+IC_list[subsector].items
            all_vector_values = all_vector_values+list(IC_list[subsector].vect)
        inds_sorted_vector_values = np.argsort(all_vector_values)
        new_sector.items = [all_residues[i] for i in inds_sorted_vector_values]
        new_sector.col = 'r'
        sectors.append(new_sector)

    return sectors


def organize_input(input_file, motif_folder):
    # Transfer all calculations and input files to the 0_InputData folder
    if input_file[-3:] != '.db':
        input_file = ''.join([input_file, '.db'])

    folderName = newFolder('/'.join([motif_folder, '0_InputData']))
    for file in os.listdir(motif_folder):
        if file.split('.')[-1] in ['db', 'fasta'] or \
                file.split('.')[-1][0] == 'o' :
            src = '/'.join([motif_folder, file])
            dest = folderName
            shutil.move(src, dest)

    return folderName


def plot_seqsimilarity_matrix(listS, Dsca, Dseq, statistics_folder):
    # Cluster the sequence similarity matrix and plot it
    Z = sch.linkage(Dsca['simMat'], method='complete', metric='cityblock')
    R = sch.dendrogram(Z, no_plot=True)
    ind = R['leaves']

    plt.figure(figsize=[18, 8])
    plt.subplot(121)
    plt.hist(listS, int(round(Dseq['Npos']/2)))
    plt.xlabel('Pairwise sequence identities', fontsize=18)
    plt.ylabel('Number', fontsize=18)
    plt.subplot(122)
    plt.imshow(Dsca['simMat'][np.ix_(ind, ind)], vmin=0, vmax=1)
    plt.colorbar()
    plt.ylabel('Sequence nr.', fontsize=18)
    plt.xlabel('Sequence nr.', fontsize=18)
    plt.tight_layout()
    plt.savefig('/'.join([statistics_folder, 'SequenceSimilarity.png']),
                transparent=True, dpi=300)
    plt.close()


def analyze_BGC_statistics(Dseq, statistics_folder):
    # Analyse from which BGCs the KSs originate, plot that and save to txt

    [annot, spec] = construct_phylogeny_and_annotations(Dseq['hd'])
    print('BGCs and number of KSs')

    num = []
    xlab = []
    FileID = open('/'.join([statistics_folder, 'BGCstatistics.txt']), 'w+')
    FileID.write('BGC\t Count\t Sequences\n')

    for k in spec.keys():
        print(k + ': ' + str(len(spec[k])) + ' motifs.')
        num.append(len(spec[k]))
        xlab.append(k)
        string = k + '\t' + str(len(spec[k])) + '\t' + str(spec[k]) + '\n'
        FileID.write(string)
    FileID.close()

    # Plot a histogram
    plt.figure(figsize=[8, 14])
    plt.barh(range(0, len(xlab)), num, tick_label=xlab, align='center',
             color=color)
    plt.xlabel('Numer of KSs', fontsize=18)
    plt.tick_params(labelsize=14)
    plt.xticks(np.arange(max(num), step=2))
    plt.ylim([-0.5, len(xlab)])
    ax = plt.gca()
    ax.yaxis.grid(True)
    ax.set_axisbelow(True)
    plt.rcParams['grid.alpha'] = 0.0
    plt.grid()
    plt.tight_layout()
    plt.savefig('/'.join([statistics_folder, 'BGCstatistics.png']),
                transparent=True, dpi=300)
    plt.close()


def plot_positional_entropy(Dsca, Dseq, first_order_folder):
    plt.figure(figsize=[9, 5])
    fig, axs = plt.subplots(1, 1, figsize=(9, 4))
    xvals = [i+1 for i in range(len(Dsca['Di']))]
    tickDist = 50
    xticks = range(0, len(Dseq['ats']), tickDist)
    plt.bar(xvals, Dsca['Di'], color=color, width=1)
    plt.tick_params(labelsize=12)
    axs.set_xticks(xticks)
    labels = [Dseq['ats'][k] for k in xticks]
    axs.set_xticklabels(labels)
    plt.xlabel('Amino acid position', fontsize=18)
    plt.ylabel('$\it{D_i}$', fontsize=18)
    plt.tight_layout()
    plt.xlim([-10, len(Dseq['ats'])+10])
    plt.tight_layout()
    plt.savefig('/'.join([first_order_folder, 'PositionalEntropy.png']),
                transparent=True, dpi=300)
    plt.close()


def plot_SCA_matrix(Dsca, first_order_folder, consSeq):
    # Plot the SCA matrix
    plt.figure(figsize=[9, 8])
    plt.rcParams['figure.figsize'] = 9, 8
    plt.imshow(Dsca['Csca'], vmin=0, vmax=1.4, interpolation='none',
               aspect='equal')
    plt.xlabel('Amino Acid Position', fontsize=14)
    plt.ylabel('Amino Acid Position', fontsize=14)
    cbar = plt.colorbar(shrink=0.5)
    cbar.set_label('Positional entropy', fontsize=14, labelpad=-55)
    plt.tight_layout()
    plt.savefig('/'.join([first_order_folder, 'SCAmatrix.png']),
                transparent=True, dpi=300)
    plt.close()

    consSeq = ''.join([i + '\t' for i in consSeq['aa']])
    np.savetxt('/'.join([first_order_folder, 'SCA_matrix.txt']), Dsca['Csca'],
               delimiter='\t', header=consSeq)


def find_significant_eigenmodes(Dsca, Dseq, Dsect):
    # Identify the eigenmodes that are different from a random background
    hist0, bins = np.histogram(Dsca['Lrand'].flatten(), bins=Dseq['Npos'],
                               range=(0, Dsect['Lsca'].max()))
    hist1, bins = np.histogram(Dsect['Lsca'], bins=Dseq['Npos'],
                               range=(0, Dsect['Lsca'].max()))
    return hist0, hist1, bins


def plot_significant_eigenmodes(Dsca, Dseq, Dsect, first_order_folder):
    # First analyze...
    hist0, hist1, bins = find_significant_eigenmodes(Dsca, Dseq, Dsect)

    # ... and then plot the significant eigenmodes
    fig, ax = plt.subplots(figsize=[9, 8])
    fig.clf()
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    plt.bar(bins[:-1], hist1, np.diff(bins), color=color)
    plt.plot(bins[:-1], hist0/Dsca['Ntrials'], color=[0.03, 0.05, 0.5],
             linewidth=2)
    plt.tick_params(labelsize=11)
    plt.ylabel('Numbers', fontsize=18)
    ax2 = fig.add_subplot(gs[1, 0])
    plt.bar(bins[:-1], hist1, np.diff(bins), color=color)
    plt.plot(bins[:-1], hist0/Dsca['Ntrials'], color=[0.03, 0.05, 0.5],
             linewidth=2)
    plt.ylim([-0.1, 2])
    plt.ylabel('Numbers', fontsize=18)
    plt.xlabel('Eigenvalues', fontsize=18)
    plt.ylabel('Numbers', fontsize=18)
    print('Number of eigenmodes to keep is %i' % (Dsect['kpos']))
    plt.savefig('/'.join([first_order_folder, 'SignificantEigenmodes.png']),
                transparent=True, dpi=300)
    plt.tight_layout()
    plt.close()


def calculate_eigenmode_residues(ICA_vectors, no_eigenmodes, Dsect,
                                 first_order_folder):
    # For every significant eigenmode, calculate which residues contribute
    # most.

    plt.figure(figsize=[8, 18])
    for mode in range(no_eigenmodes):
        iqr = (scoreatpercentile(ICA_vectors[:, mode], 75) -
               scoreatpercentile(ICA_vectors[:, mode], 25))
        binwidth = 2*iqr*(len(ICA_vectors)**(-0.33))
        nbins = int(round((max(ICA_vectors[:, mode]) -
                           min(ICA_vectors[:, mode])) / binwidth))
        plt.subplot(no_eigenmodes, 1, mode+1)
        h_params = plt.hist(ICA_vectors[:, mode], nbins)
        x_dist = np.linspace(min(h_params[1]), max(h_params[1]), num=100)
        plt.plot(x_dist, Dsect['scaled_pd'][mode], color=[0.03, 0.05, 0.5],
                 linewidth=2)
        plt.plot([Dsect['cutoff'][mode], Dsect['cutoff'][mode]], [0, 60],
                 'k--', linewidth=1, color=color)
        plt.xlabel(r'$V^p_{%i}$' % (mode+1), fontsize=14)
        plt.ylabel('Number', fontsize=14)
    plt.tight_layout()
    plt.savefig('/'.join([first_order_folder, 'EigenvalueFitsICs.png']),
                transparent=True, dpi=300)
    plt.close()


def plot_sector_organized_SCA_matrix(reorderedSCA, ic_sizes, sec_order_folder):
    # Replot the SCA matrix ordered by segments
    plt.figure(figsize=[20, 20])
    plt.imshow(reorderedSCA, vmin=0, vmax=2, interpolation='none',
               aspect='equal', extent=[0, sum(ic_sizes),
                                       0, sum(ic_sizes)])
    tickPos = []
    tickLabel = []
    for i, j in enumerate(ic_sizes):
        if i == 0:
            tickPos.append(ic_sizes[i]/2)
            tickLabel.append(i+1)
        else:
            tickPos.append(float(ic_sizes[i-1]/2 + tickPos[i-1] +
                                 float(ic_sizes[i]/2)))
            tickLabel.append(i+1)

    line_index = 0
    for i in range(len(ic_sizes)):
        plt.plot([line_index+ic_sizes[i], line_index+ic_sizes[i]],
                 [0, sum(ic_sizes)], 'w', linewidth=1)
        plt.plot([0, sum(ic_sizes)],
                 [sum(ic_sizes) - line_index, sum(ic_sizes)-line_index],
                 'w', linewidth=1)
        line_index += ic_sizes[i]
    plt.tight_layout()
    plt.savefig('/'.join([sec_order_folder, 'SCAmatrix_ICordered.png']),
                transparent=True, dpi=300)
    plt.close()


def plot_reorganized_SCA_matrix(reordered_sectors, Dsca,
                                icsizes, sec_order_folder):
    # plot the re-ordered matrix
    sortpos = list()
    for sector in reordered_sectors:
        sortpos.extend(sector.items)

    reordered_Csca = Dsca['Csca'][np.ix_(sortpos, sortpos)]

    plt.figure(figsize=[10, 10])
    line_index = 0
    plt.imshow(reordered_Csca, vmin=0, vmax=2.2, interpolation='none',
               aspect='equal', extent=[0, len(sortpos), 0, len(sortpos)])

    # Plot lines at sector boundaries
    for sector in reordered_sectors:
        plt.plot([line_index+len(sector.items), line_index+len(sector.items)],
                 [0, len(sortpos)], 'w', linewidth=1)
        plt.plot([0, sum(icsizes)], [len(sortpos)-line_index,
                                     len(sortpos)-line_index], 'w',
                 linewidth=1)
        line_index += len(sector.items)
    plt.tight_layout()
    try:
        plt.savefig('/'.join([sec_order_folder, 'SCAmatrix_Reordered.png']),
                    transparent=True, dpi=300)
    except Exception:
        traceback.print_exc()
        print('ERROR: Saving the reoreded matrix at 300 dpi fialed!\n\n')
    try:
        plt.savefig('/'.join([sec_order_folder, 'SCAmatrix_Reordered.png']))
    except Exception:
        traceback.print_exc()
        print('ERROR: Saving the reoreded matrix in png failed!\n\n')
    try:
        plt.savefig('/'.join([sec_order_folder, 'SCAmatrix_Reordered.jpg']))
    except Exception:
        traceback.print_exc()
        print('ERROR: Saving the reoreded matrix in jpg failed!\n\n')
    plt.close()


def plot_positional_coupling_and_conservation(allignedSectors, Dseq, consSeq,
                                              sec_order_folder, reordered):
    # Plot the positional correlation strength
    max_coupling = 0
    for rank in allignedSectors:
        max_couplingSector = max(allignedSectors[rank].integerCoupling())
        if max_couplingSector > max_coupling:
            max_coupling = max_couplingSector

    if 'color' not in locals():
        color = [0.1, 0.6, 0.3]

    try:
        plt.figure(figsize=[2+0.3*len(Dseq['alg'][0]), 4*len(allignedSectors)])
        for rank in allignedSectors:
            ax1 = plt.subplot(len(allignedSectors)*2, 1, rank)
            plt.ylabel(' '.join(['Sector', str(rank)]))
            ax2 = ax1.twinx()
            ax1.bar(range(len(allignedSectors[rank].pos)),
                    allignedSectors[rank].coupling, width=0.9, color=color,
                    label='Coupling')
            ax2.bar(range(len(allignedSectors[rank].pos)),
                    [x for x in allignedSectors[rank].integerCons()],
                    width=0.9, color=[0.4, 0.5, 0.8], label='Conservation')
            plt.xticks(range(0, len(consSeq['aa']), 10), fontsize=12)
            plt.ylabel(''.join(['Sector', str(rank)]))
            ax1.set_xlim([-1, len(consSeq['aa'])])
            ax1.set_ylim([-max_coupling*1.1, 1.1*max_coupling])
            ax2.set_ylim([1.1, -1.1])
            plt.legend(frameon=False)
            if rank == 1:
                ax2 = ax1.twiny()
                ax2.set_xlim(ax1.get_xlim())
                ax2.set_xticks(range(0, len(consSeq['aa'])))
                ax2.set_xticklabels(consSeq['aa'])
        plt.tight_layout()
        if reordered:
            plt.savefig('/'.join([sec_order_folder,
                                  'PositionalCoupling_Reordered.png']),
                        transparent=True)
            plt.savefig('/'.join([sec_order_folder,
                                  'PositionalCoupling_Reordered.jpg']))
        else:
            plt.savefig('/'.join([sec_order_folder,
                                  'PositionalCoupling_Unordered.png']),
                        transparent=True)
            plt.savefig('/'.join([sec_order_folder,
                                  'PositionalCoupling_Unordered.jpg']))
        plt.close()
    except Exception:
        traceback.print_exc()
        print('ERROR: I did not print the positional coupling and conservation!\n\n' )


def plot_coupling_in_original_sequence(allignedSectors, Dsca,
                                       sec_order_folder, consSeq, reordered):
    # Map the couplings of every residue to the original sequence alignment
    consSeq = ''.join([i + '\t' for i in consSeq['aa']])
    try:
        for rank in allignedSectors:
            EntropyMat = np.empty(np.shape(Dsca['Csca']))
            EntropyMat[:] = np.nan
            for rownr, row in enumerate(EntropyMat):
                for col in range(len(row)):
                    if rownr in allignedSectors[rank].pos\
                            and col in allignedSectors[rank].pos:
                        # Csca is symmetric
                        EntropyMat[rownr, col] = Dsca['Csca'][rownr, col]
                        EntropyMat[col, rownr] = Dsca['Csca'][col, rownr]

            # Save to txt
            if reordered:
                fID = open('/'.join([sec_order_folder, 'SectorMatrix_sector' +
                                     str(rank)+'_Reordered.txt']), 'w+')
            else:
                fID = open('/'.join([sec_order_folder, 'SectorMatrix_sector' +
                                     str(rank)+'_Unordered.txt']), 'w+')

            fID.write('\t' + consSeq + '\n')
            for i, row in enumerate(EntropyMat):
                string = consSeq[2*i] + '\t'
                for j in row:
                    if np.isnan(j):
                        string += '0\t'
                    else:
                        string += '%.7f\t' % j
                string += '\n'
                fID.write(string)

            plt.figure(figsize=[10, 10])
            plot = plt.imshow(EntropyMat, vmin=0, vmax=1.4,
                              interpolation='none', aspect='equal')
            plot.axes.grid(which='major', linestyle='-', color='k', alpha=0.2)
            plt.tight_layout()
            if reordered:
                plt.savefig('/'.join([sec_order_folder, 'SectorMatrix_sector' +
                                      str(rank)+'_Reordered.png']),
                            transparent=True, dpi=600)
            else:
                plt.savefig('/'.join([sec_order_folder, 'SectorMatrix_sector' +
                                      str(rank)+'_Unordered.png']),
                            transparent=True, dpi=600)
            plt.close()
    except Exception:
        traceback.print_exc()
        print('ERROR: I did not print the coupling in the original matrix in png.')
        print('Failed at sector ' + str(rank) + ' of ' +
              str(len(allignedSectors)) + '\n\n')
    try:
        for rank in allignedSectors:
            EntropyMat = np.empty(np.shape(Dsca['Csca']))
            EntropyMat[:] = np.nan
            for rownr, row in enumerate(EntropyMat):
                for col in range(len(row)):
                    if rownr in allignedSectors[rank].pos\
                            and col in allignedSectors[rank].pos:
                        # Csca is symmetric
                        EntropyMat[rownr, col] = Dsca['Csca'][rownr, col]
                        EntropyMat[col, rownr] = Dsca['Csca'][col, rownr]
            # Save to txt
            np.savetxt('/'.join([sec_order_folder, 'SectorMatrix_sector' +
                                 str(rank)+'.txt']), EntropyMat,
                       delimiter='\t', header=consSeq)
            plt.figure(figsize=[10, 10])
            plot = plt.imshow(EntropyMat, vmin=0, vmax=1.4,
                              interpolation='none', aspect='equal')
            plot.axes.grid(which='major', linestyle='-', color='k', alpha=0.2)
            plt.tight_layout()
            if reordered:
                plt.savefig('/'.join([sec_order_folder, 'SectorMatrix_sector' +
                                  str(rank)+'_Reordered.jpg']))
            else:
                plt.savefig('/'.join([sec_order_folder, 'SectorMatrix_sector' +
                                  str(rank)+'_Unordered.jpg']))
            plt.close()
    except Exception:
        traceback.print_exc()
        print('ERROR: I did not print the coupling in the original matrix in JPG.')
        print('Failed at sector ' + str(rank) + ' of ' +
              str(len(allignedSectors)) + '\n\n')


def writePymol(pdb, sectors, outfilename, color, offset, chain="A", quit=1):
    """
    ****TAKEN FROM THE ORIGINAL RAMAGANDRAN SOFTWARE****
    Write basic a pymol script for displaying sectors and exporting an image.

    **Example**::

      writePymol(pdb, sectors, ics, ats, outfilename, chain='A',inpath=settings.path2structures, quit=1)
  
    """

    f = open(outfilename, "w")
    f.write("delete all\n")
    f.write("load %s.pdb, main\n" % pdb)
    f.write("hide all\n")
    f.write("bg_color white\n")
    f.write("show cartoon\n")
    f.write("color hydrogen\n\n")
    b, g, r = color
    # f.write("set_color color1, [%.3f,%.3f,%.3f]\n" % (b, g, r))
    # PyMol is 1-indexed
    for sector in sectors:
        f.write("create sector" + str(sector) + ", (resi %s) \n"
                % tuple([", resi ".join([str(pos+1-offset) 
                                         for pos in sectors[sector]
                                         if pos >= offset])]))
    # f.write("color color1, sector1\n")
    f.write("spectrum cyan_yellow_red, sector1\n")
    f.write("show spheres, sector1\n")
    f.write("show surface, sector1\n\n")
    f.write("zoom\n")
    f.write("set transparency, 0.4\n")
    f.close()
    f = open(outfilename.replace('.txt', '_save.txt'), 'w')
    f.write("ray\n")
    path_list = outfilename.split(os.sep)
    fn = path_list[-1]
    f.write("png %s\n" % fn.replace(".pml", ""))
    f.write("set transparency, 0.15\n")
    f.write('set surface_quality, 2\n')
    f.write('set cartoon_quality, 25\n')
    f.write('set geometry_export_mode, 1\n')
    f.write("save " + ''.join(fn.split('.')[0:-1]) + '.dae\n')
    if quit == 1:
        f.write("quit")
    f.close()
    

def map_to_original_sequence(motif, motif_folder, Dseq, model):
    # Map the SCA ats mapping back onto the original, prefiltered sequence and
    # reextract the original sequence

    # It's not the most beautiful way of doing it, but it works..

    MSA_file = [i for i in os.listdir(os.path.join(motif_folder,
                                                   '0_InputData'))
                if 'MSA_' + motif + '.o' in i][0]
    headers_full, sequences_full = sca.readAlg(motif_folder +
                                               '/0_InputData/' + MSA_file)

    for i in range(len(headers_full)):
        if Dseq['hd'][model] == headers_full[i]:
            index = i
            break
    pos_in_sca = 0
    full_seq_model = sequences_full[index]
    sca_seq = Dseq['alg'][model]
    full_seq = ''
    fil_in_orig = 0
    insert_count = 0
    # The new_ats maps the (map of the filtered sequence used for SCA) to the
    # positions in the original ClustalW sequence alingment
    filt_to_orig = dict()
    orig_to_filt = dict()
    for pos_in_seq, let in enumerate(full_seq_model):
        if let == '-':
            insert_count += 1
        else:
            full_seq += let
            x = str(fil_in_orig) + ' - ' + let
            while sca_seq[pos_in_sca] == '-':
                if pos_in_sca >= len(sca_seq)-1:
                    break
                pos_in_sca += 1
            if let == sca_seq[pos_in_sca]:
                x += ' : ' + sca_seq[pos_in_sca] + ' - ' + str(pos_in_sca)
                filt_to_orig[pos_in_sca] = fil_in_orig
                orig_to_filt[fil_in_orig] = pos_in_sca
                pos_in_sca += 1
            else:
                orig_to_filt[fil_in_orig] = None
            #if pos_in_seq <400:
            #    print(x)
            if pos_in_sca > len(sca_seq)-1:
                break
            fil_in_orig += 1

    return [full_seq, orig_to_filt, filt_to_orig]


def map_to_original_msa(motif, motif_folder, Dseq):
    # Map the SCA ats to the original sequence alignment
    MSA_file = [i for i in os.listdir(os.path.join(motif_folder,
                                                   '0_InputData'))
                if 'MSA_' + motif + '.o' in i][0]
    headers_full, sequences_full = sca.readAlg(motif_folder +
                                               '/0_InputData/' + MSA_file)

    sca_to_msa = dict()
    msa_to_sca = dict()
    mapping = []

    for pos_in_sca, let in enumerate(Dseq['alg'][0]):
        sca_to_msa[pos_in_sca] = None
    for pos_in_msa, let in enumerate(sequences_full[1]):
        msa_to_sca[pos_in_msa] = None
    seqs_in_sca = []
    for i, seq in enumerate(headers_full):
        index = None
        for j in range(len(Dseq['hd'])):
            if Dseq['hd'][j] == headers_full[i]:
                seqs_in_sca.append(i)
                break

    for sca_pos, msa_pos in enumerate(seqs_in_sca):
        # if sca_pos>3:
            # break
        seq_in_msa = sequences_full[msa_pos]
        try:
            sca_seq = Dseq['alg'][sca_pos]
        except Exception:
            continue
        insert_count = 0
        pos_in_sca = 0
        for pos_in_seq, let in enumerate(seq_in_msa):
            if let == '-':
                insert_count += 1
            else:
                while sca_seq[pos_in_sca] == '-':
                    if pos_in_sca >= len(sca_seq)-1:
                        break
                    pos_in_sca += 1
                if let == sca_seq[pos_in_sca]:
                    if sca_to_msa[pos_in_sca] is None:
                        sca_to_msa[pos_in_sca] = pos_in_seq
                    pos_in_sca += 1
            if pos_in_sca >= len(sca_seq) - 1:
                break
        mapping = [sca_to_msa[i] for i in sca_to_msa.keys()]
        if None not in mapping:
            break
    for sca_pos in sca_to_msa.keys():
        msa_to_sca[sca_to_msa[sca_pos]] = sca_pos

    return sca_to_msa, msa_to_sca


def map_sectors_to_lost_sequence(sca_to_msa, allignedSectors, 
                                 motif_folder, motif):
    # Create a dictionary that maps the msa position to the position in the
    # actual original sequence  
    MSA_file = [i for i in os.listdir(os.path.join(motif_folder,
                                                   '0_InputData'))
                if 'MSA_' + motif + '.o' in i][0]
    headers_full, sequences_full = sca.readAlg(motif_folder +
                                               '/0_InputData/' + MSA_file)

    for index, header in enumerate(headers_full):
        if 'Oocydin_Serratia | gene: ctg1_16 | start, end in domain sequence: [4, 5] | . '\
            in header:
                print(header)
                break
    msa_seq = sequences_full[index]

    msa_to_orig = dict()
    orig_to_msa = dict()
    pos_in_orig = 0
    orig = ''
    for i, let in enumerate(msa_seq):
        if let == '-':
            continue
        else:
            msa_to_orig[i] = pos_in_orig
            orig_to_msa[pos_in_orig] = i
            pos_in_orig += 1
            orig += let

    # Now map the sector on the original sequences
    sectorpos_in_orig = dict()
    sectorseqs = dict()
    for rank in allignedSectors:
        positions_in_orig = []
        sectorseq = ''
        for i in allignedSectors[rank].pos:
            if i == '-':
                sectorseq += i
                continue
            try:
                pos_in_msa = sca_to_msa[i]
            except KeyError:
                print('KeyError for sca_to_msa[' + str(i) + '] in rank ' + str(rank))
                continue
            try:
                pos_in_orig = msa_to_orig[pos_in_msa]
                positions_in_orig.append(pos_in_orig)
                sectorseq += orig[pos_in_orig]
            except KeyError:
                print('KeyError for msa_to_orig[' + str(i) + '] in rank ' + str(rank))
                pass
        sectorpos_in_orig[rank] = positions_in_orig
        sectorseqs[rank] = sectorseq

    return sectorpos_in_orig, sectorseqs


def SpyderInput():
    motifs = ['KS_trans-ATdocking','ACP_KS_trans-ATdocking', 'ACP_KS', 'KR_ACP_KS','KS',
               'KS_trans-ATdocking_KR',
              'AT','Condensation']
    #motifs = ['KS_trans-ATdocking', 'KS_trans-ATdocking_KR', 'AT','Condensation']
    cis_motifs = ['ACP_KS_AT', 'KS', 'KS_AT', 'KS_AT_KR']
    for motif in motifs:
        try:
            print(motif)
            input_folder = 'Z:/Scripts/S-005_SCA_Pipeline/Results/Results_2021-04-15_FullSequence/' + motif
            motif_folder = input_folder
            input_folder += '/0_InputData'
            # motif = 'KS_AT_KR'
            input_file = motif + '_sectors.db'
            import sys
            import os
            from importlib import reload
            sys.path.append('../_sca')
            sys.path.append('./_sca')
            sys.path.append('./_analysis')
            from pysca import sca
            #from KSanalysis import *
            import KStools
            from KSanalysis import *
            #reload(sys.modules['KSanalysis'])
            import pickle
            #motif_folder = os.getcwd()

            db = pickle.load(open('/'.join([input_folder, input_file]), 'rb'))
            Dseq = db['sequence']  # the results of scaProcessMSA
            Dsca = db['sca']       # the results of scaCore
            Dsect = db['sector']   # the results of scaSectorID

            consSeq = KStools.constructConsensus(Dseq['alg'])
            sca_to_msa, msa_to_sca = map_to_original_msa(motif, motif_folder, Dseq)

            sector_organizedSCA = Dsca['Csca'][np.ix_(Dsect['sortedpos'],
                                                      Dsect['sortedpos'])]
            ic_sizes = Dsect['icsize']
            sec_order_folder = newFolder('/'.join([motif_folder,'3_ReorderedSCAMatrix']))
            threshold = 1
            biggest_smaller_than_half = False
            while not biggest_smaller_than_half:
                [couplingMat, sec_groups, biggest_smaller_than_half] = \
                    KStools.find_correlated_sectors(sector_organizedSCA, ic_sizes,
                                                    threshold)
                threshold = threshold*1.01

            reordered_sectors = obtain_reordered_sector_data(sec_groups, Dsect['ics'])

            allignedSectors = dict()
            for rank, sector in enumerate(reordered_sectors):
                allignedSectors[rank+1] = KStools.allignedSector(rank+1, sector, consSeq, Dsca['Csca'], sca_to_msa)

            KStools.saveAlignment(allignedSectors, consSeq, sec_order_folder,
                              sca_to_msa, msa_to_sca)

            plot_positional_coupling_and_conservation(allignedSectors, Dseq, consSeq,
                                                      sec_order_folder)
            plot_coupling_in_original_sequence(allignedSectors, Dsca, sec_order_folder,
                                               consSeq)
        except Exception as e:
            print('Exception for ' + motif)
            print(e)

def main():
    motif = sys.argv[1]
    input_file = motif + '_sectors.db'
    sys.path.append(sys.argv[1])
    import pysca

    os.chdir('..')
    os.chdir(motif)
    # Main function to analyse the SCA results
    motif_folder = os.getcwd()
    global color
    color = [0.1, 0.6, 0.3]

    input_folder = organize_input(input_file, motif_folder)

    # Load the results from the SCA core calculations and show some statistics
    db = pickle.load(open('/'.join([input_folder, input_file]), 'rb'))
    Dseq = db['sequence']  # the results of scaProcessMSA
    Dsca = db['sca']       # the results of scaCore
    Dsect = db['sector']   # the results of scaSectorID

    consSeq = KStools.constructConsensus(Dseq['alg'])
    sca_to_msa, msa_to_sca = map_to_original_msa(motif, motif_folder, Dseq)

    print("After processing, the alignment size is %i sequences and %i positions" %
          (Dseq['Nseq'], Dseq['Npos']))
    print("With sequence weights, there are %i effective sequences" %
          (Dseq['effseqs']))

    # Create sequence similarity matrix
    statistics_folder = newFolder('/'.join([motif_folder,'1_GeneralStatistics']))
    listS = [Dsca['simMat'][i, j] for i in range(Dsca['simMat'].shape[0])
             for j in range(i+1, Dsca['simMat'].shape[1])]
    plot_seqsimilarity_matrix(listS, Dsca, Dseq, statistics_folder)

    # Analyze BGC statistics and phylogenies
    analyze_BGC_statistics(Dseq, statistics_folder)

    # First order statistics
    first_order_folder = newFolder('/'.join([motif_folder,'2_FirstOrderStatistics']))
    plot_positional_entropy(Dsca, Dseq, first_order_folder)
    plot_SCA_matrix(Dsca, first_order_folder, consSeq)
    plot_significant_eigenmodes(Dsca, Dseq, Dsect, first_order_folder)
    calculate_eigenmode_residues(Dsect['Vpica'], Dsect['kpos'], Dsect,
                                 first_order_folder)

    # Second order statistics
    sec_order_folder = newFolder('/'.join([motif_folder,'3_ReorderedSCAMatrix']))
    sector_organizedSCA = Dsca['Csca'][np.ix_(Dsect['sortedpos'],
                                              Dsect['sortedpos'])]
    ic_sizes = Dsect['icsize']
    plot_sector_organized_SCA_matrix(sector_organizedSCA, ic_sizes,
                                     sec_order_folder)

    # FIRST: SAVE THE DATA BEFORE REORDERING THE SECTORS
    # Align the sequences to the consensus sequence
    reordered = False

    allignedSectors = dict()
    for rank, sector in enumerate(Dsect['ics']):
        allignedSectors[rank+1] = KStools.allignedSector(rank+1, sector,
                                                         consSeq, Dsca['Csca'],
                                                         sca_to_msa)

    KStools.saveAlignment(allignedSectors, consSeq, Dseq['atsConsFiltering'],
                          sec_order_folder, sca_to_msa, msa_to_sca, reordered)

    # Finally, plot the positional coupling and conservation in 2 and 3D
    plot_positional_coupling_and_conservation(allignedSectors, Dseq,
                                              consSeq, sec_order_folder,
                                              reordered)
    plot_coupling_in_original_sequence(allignedSectors, Dsca, sec_order_folder,
                                       consSeq, reordered)

    # Reorder the ICs to cluster correlated ICs
    threshold = 1
    biggest_smaller_than_half = False
    while not biggest_smaller_than_half:
        [couplingMat, sec_groups, biggest_smaller_than_half] = \
            KStools.find_correlated_sectors(sector_organizedSCA, ic_sizes,
                                            threshold)
        threshold = threshold*1.01
    fileID = open(os.path.join(sec_order_folder, 'ReorderingThreshold.txt'),
                  'w+')
    fileID.write('The reordering threshold was: ' + str(threshold))
    fileID.close()
    reordered = True

    reordered_sectors = obtain_reordered_sector_data(sec_groups, Dsect['ics'])

    plot_reorganized_SCA_matrix(reordered_sectors, Dsca, ic_sizes,
                                sec_order_folder)

    # Align the sequences to the consensus sequence
    allignedSectors = dict()
    for rank, sector in enumerate(reordered_sectors):
        allignedSectors[rank+1] = KStools.allignedSector(rank+1, sector,
                                                         consSeq, Dsca['Csca'],
                                                         sca_to_msa)

    KStools.saveAlignment(allignedSectors, consSeq, sec_order_folder,
                          sca_to_msa, msa_to_sca, reordered)

    # Finally, plot the positional coupling and conservation in 2 and 3D
    plot_positional_coupling_and_conservation(allignedSectors, Dseq, consSeq,
                                              sec_order_folder, reordered)
    plot_coupling_in_original_sequence(allignedSectors, Dsca, sec_order_folder,
                                       consSeq, reordered)

    # Make the pymol scripts
    models = []
    for i, header in enumerate(Dseq['hd']):
        if 'oocydin ' in header.lower() or 'bacillaene' in header.lower():
            models.append(i)
    if 'cis' in motif_folder:
        for i, header in enumerate(Dseq['hd']):
            if 'erythr' in header.lower():
                print(header)
                models.append(i)

    pymol_folder = newFolder('/'.join([motif_folder, '4_PyMolFiles']))

    for model in models:
        # Find index in parsed alignment
        try:
            [full_seq,  orig_to_filt, filt_to_orig] = map_to_original_sequence(
                motif, motif_folder, Dseq, model)
        except IndexError as e:
            print('ERROR: IndexError during processing of ' + Dseq['hd'][model])
            print(e)
            print('-------------')
            continue
        gene = Dseq['hd'][model].split('|')[1].split(':')[1].strip()
        if gene is '':
            gene = motif
        gene = Dseq['hd'][model].split('|')[0].strip() + '_'  + \
            gene
        start = Dseq['hd'][model].split('|')[2].split(':')[1].strip()
        start = start.split(',')[0][1:]

        print('INFO: Making Pymol files of:')
        print('           ' + Dseq['hd'][model])
        pdb = gene + '_' + start
        try:
            fID = open(pdb + '_seqPos.txt', 'w+')
            fID.write('Pos in seq\tAA\tPosInFilteredSeq\tAA\n')
            for i, let in enumerate(full_seq):
                x = str(i) + '\t' + let
                if orig_to_filt[i] != None:
                    x += '\t' + str(orig_to_filt[i]) + '\t'  + Dseq['alg'][model][orig_to_filt[i]]
                x += '\n'
                fID.write(x)
            fID.close()
            color_i = [15, 161, 213]
            color_f = [239, 157, 1]
            numSect = len(allignedSectors)
            fname = '/'.join([pymol_folder, pdb + '.txt'])
            sectors = dict()
            offset = 0
            for rank in allignedSectors:
                pos_in_sectors = [i for i in allignedSectors[rank].pos if i!= '-']
                sectors[rank] = [filt_to_orig[i] for i in pos_in_sectors 
                             if i in filt_to_orig.keys()]
                writePymol(pdb, sectors, fname, color, offset, 'A', True)
        except KeyError as e:
            print('ERROR: KeyError during PyMol writing of ' + pdb)
            print(e)
            print('---------')
    return
if __name__ == '__main__':
    main()