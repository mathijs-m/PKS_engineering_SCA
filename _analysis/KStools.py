# -*- codposng: utf-8 -*-
"""
Collection of functions to analyze the KS analyses obtained from the SCA analysis

Mathijs Mabesoone, Piel group, February 2021
"""
import pandas as pd
import os
import numpy as np
import statistics


class Unit:
    """
    A class for units (sectors, sequence families, etc.)

    **Attributes**

      :name:  string describing the unit (ex: 'firmicutes')
      :items: set of member items (ex: indices for all firmicutes
              sequence in an alignment)
      :col:   color code associated to the unit (for plotting)
      :vect:  an additional vector describing the member items (ex: a list
              of sequence weights)
    """

    def __init__(self):
        self.name = ""
        self.items = set()
        self.col = 0
        self.vect = 0


class Annot:
    """
    A class for annotating sequences

    **Attributes**

      :desc:    description (often the sequence header)
      :species: species string
      :taxo:    taxonomy string
      :seq:     sequence
    """

    def __init__(self, descr, species, taxo, seq=""):
        self.descr = descr
        self.species = species
        self.taxo = taxo
        self.seq = seq


class allignedSector:
    def __init__(self, rank, sector, consSeq, Csca, sca_to_msa):
        [aa, cons, pos] = sectorAlignment(sector, consSeq)
        self.rank = rank
        self.aa = aa
        self.cons = cons
        self.pos = pos
        self.coupling = couplingStrength(pos, Csca)
        pos_in_msa = []
        for i in pos:
            try:
                pos_in_msa.append(sca_to_msa[i])
            except KeyError:
                pos_in_msa.append('-')
        self.pos_in_msa = pos_in_msa

    def integerCons(self):
        # Convert the dashes in conservation array to zeros
        intCons = []
        for i in range(len(self.cons)):
            if isinstance(self.cons[i], str):
                intCons.append(0)
            else:
                intCons.append(self.cons[i])
        return intCons

    def integerCoupling(self):
        intCoupling = []
        for i in range(len(self.coupling)):
            if isinstance(self.coupling[i], str):
                intCoupling.append(0)
            else:
                intCoupling.append(self.coupling[i])
        return intCoupling


def consensusAA(d, AAcount):
    """
    Find the maximum value in a dictionary
    a) create a list of the dict's keys and values;
    b) return the key with the max value"""
    v = list(d.values())
    k = list(d.keys())
    consAA = k[v.index(max(v))]
    conserv = d[consAA]/AAcount
    return [consAA, conserv]


def constructConsensus(MSA):
    # Construct a consensus sequence of the MSA based on the most encountered
    # amino acid

    # Returns a dictionary containing:
    #           1) 'aa', the consensus amino acid sequence
    #           2) 'cons', the fraction of sequences with that amino acid at
    #               the position

    consensus = dict()
    for pos in range(len(MSA[0])):
        AAcount = dict()
        for entry in range(len(MSA)):
            AA = MSA[entry][pos]
            if AA in AAcount:
                AAcount[AA] += 1
            else:
                AAcount[AA] = 1

        cons = consensusAA(AAcount, len(MSA))
        if 'aa' in consensus:
            consensus['aa'].append(cons[0])
            consensus['cons'].append(cons[1])
        else:
            consensus['aa'] = [cons[0]]
            consensus['cons'] = [cons[1]]
    consensus['aa'] = '' .join(consensus['aa'])
    return consensus


def sectorAlignment(sector, consSeq):
    # This function saves the sectors in a fasta-like file and also saves a
    # tab-delimited file with positional conservation

    # Returns a dictionary containing the sector sequences aligned against the
    # consensus sequence and the conservation of these positions

    aminoAcids = []
    posCons = []
    posInSector = []

    for pos in range(len(consSeq['aa'])):
        if pos in sector.items:
            aminoAcids.append(consSeq['aa'][pos])
            posCons.append(consSeq['cons'][pos])
            posInSector.append(pos+1)
        else:
            aminoAcids.append('-')
            posCons.append('-')
            posInSector.append('-')
    aminoAcids = ''.join(aminoAcids)
    return aminoAcids, posCons, posInSector


def couplingStrength(positions, Csca):
    # Calculates the sum of all coupling strengths in the Csca matrix for 
    # residues in the positions array

    # Returns couplingArr, an array containing all the coupling strengths

    coupling = []
    for i, entry in enumerate(positions):
        if str(entry).isnumeric():
            coupling.append(sum(Csca[i, :]))
        else:
            coupling.append(0)

    return coupling


def saveAlignment(allignedSectors, consSeq, folderName, sca_to_msa,
                  msa_to_sca, reordered):
    # This function saves all alligned sectors contained in the allignedSector
    # dictionary
    # Also saves a txt file with positional conservation in the MSA for every
    # position in the sector

    if not os.path.isdir(folderName):
        os.mkdir(folderName)

    if reordered:
        alignment_fname = '/'.join([folderName + '/SectorAlignment_Reordered.txt'])
        colSeq_fname = '/'.join([folderName + '/ColoredConsSeq_Reordered.rtf'])
        sectorCons_fname = '/'.join([folderName + '/SectorConservation_Reordered.txt'])
        sectorCoupling_fname = '/'.join([folderName + '/SectorCoupling_Reordered.txt'])
        sectorCoupling_full_fname = '/'.join([folderName + '/SectorCoupling_fullMSA_Reordered.txt'])
    if not reordered:
        alignment_fname = '/'.join([folderName + '/SectorAlignment_Unordered.txt'])
        colSeq_fname = '/'.join([folderName + '/ColoredConsSeq_Unordered.rtf'])
        sectorCons_fname = '/'.join([folderName + '/SectorConservation_Unordered.txt'])
        sectorCoupling_fname = '/'.join([folderName + '/SectorCoupling_Unordered.txt'])
        sectorCoupling_full_fname = '/'.join([folderName + '/SectorCoupling_fullMSA_Unordered.txt'])

    seqFile = open(alignment_fname, 'w+')
    seqFile.write('>SectorAlignment\n' + ''.join(consSeq['aa']) + '\n')
    consData = [[i for i in range(1, len(consSeq['aa'])+1)]]
    couplingData = [[i for i in range(1, len(consSeq['aa'])+1)]]
    consData.append([sca_to_msa[i] for i in range(len(consSeq['aa']))])
    consData.append([i for i in consSeq['aa']])
    couplingData.append([i for i in consSeq['aa']])

    colorTable = '{\colortbl;\\red0\green0\\blue0;\\red10\green40\\blue215;'+\
        '\\red120\green200\\blue205;\\red30\green225\\blue30;' + \
                '\\red255\green50\\blue205;\\red205\green0\\blue60;' + \
                '\\red200\green255\\blue0;\\red255\green255\\blue85;' + \
                '\\red0\green30\\blue128;\\red0\green128\\blue128;' + \
                '\\red0\green128\\blue0;\\red128\green0\\blue128;' + \
                '\\red128\green0\\blue0;\\red128\green128\\blue0;' + \
                '\\red128\green128\\blue128;\\red192\green192\\blue192;}'

    coloredSeq = dict()
    for i, let in enumerate(consSeq['aa']):
        coloredSeq[i+1] = {'let': let, 'col': '\cf1 '}
    try:
        for rank in allignedSectors:
            seqFile.write('>Sector ' + str(allignedSectors[rank].rank) + '\n' +
                          allignedSectors[rank].aa + '\n')
            consData.append(allignedSectors[rank].cons)
            couplingData.append(allignedSectors[rank].coupling)
            # Also write the colors for the sequence
            for pos in allignedSectors[rank].pos:
                if pos != '-':
                    coloredSeq[int(pos)]['col'] = '\cf' + str(rank+1) + ' '
    
        # Write the rtf file first
        print('writing at ' + colSeq_fname)
        fID = open(colSeq_fname, 'w+') 
        fID.write('{\\rtf1 \n' + colorTable + '\n')
        for pos in coloredSeq.keys():
            fID.write(coloredSeq[pos]['col'] + coloredSeq[pos]['let'])
        fID.write('\line \line')
        for rank in allignedSectors:
            fID.write('\cf' + str(rank+1) + ' Sector ' + str(rank) + '\line')
    
        fID.write('}')
        fID.close()
    except Exception:
        pass
    colnames = dict()
    colnames[0] = 'Position'
    colnames[2] = 'Consensus AA'
    for i in range(len(allignedSectors)):
        colnames[i+2] = str(i + 1)

    consOverview = pd.DataFrame.from_dict(data=consData).transpose().rename(
        columns=colnames)
    consOverview.to_csv(sectorCons_fname, sep='\t', index=0)

    couplingOverview = pd.DataFrame.from_dict(data=couplingData).transpose().rename(columns=colnames)
    couplingOverview.to_csv(sectorCoupling_fname, sep='\t', index=0)
    seqFile.close()


    colnames = dict()
    colnames[0] = 'Position'
    colnames[1] = 'PositionInMSA'
    colnames[2] = 'Consensus AA'
    for i in range(len(allignedSectors)):
        colnames[i+3] = str(i + 1)

    seqFile = open(sectorCoupling_full_fname, 'w+')
    seqFile.write('PosInMSA\tPosInSCA\tConsensusAA\t' +
                  '\t'.join([str(allignedSectors[i].rank) for i in allignedSectors]) + '\n')
    msa_positions = [i for i in msa_to_sca.keys()]
    for pos in msa_positions:
        string = str(pos)
        sca_pos = msa_to_sca[pos]
        if sca_pos != None:
            string += '\t' + str(msa_to_sca[pos]) + '\t' + consSeq['aa'][sca_pos]
            for rank in allignedSectors:
                if sca_pos in allignedSectors[rank].pos:
                    string += '\t' + str(allignedSectors[rank].coupling[sca_pos])
                else:
                    string += '\t-'
        string += '\n'
        seqFile.write(string)
    seqFile.close()

def find_correlated_sectors(sectMat, icsize, th):
    # This function finds the correlations between sectors in the SCA matrix
    # The selection is based on the average value in a non-diagonal domain 
    # being larger than the median of all non-diagonal domains
    off_diag = []
    ind_i = 0

    for i in range(len(icsize)):
        ind_j = 0
        for j in range(i):
            i_range = slice(ind_i, ind_i+icsize[i])
            j_range = slice(ind_j, ind_j+icsize[j])
            off_diag.extend(sectMat[i_range, j_range].flatten())
            ind_j += icsize[j]
        ind_i += icsize[i]

    med = np.median(off_diag)
    coupled_mat = np.zeros([len(icsize), len(icsize)])
    ind_i = 0

    avgMat = np.zeros(np.shape(sectMat))
    for i in range(len(icsize)):
        ind_j = 0
        for j in range(i+1):
            i_range = slice(ind_i, ind_i+icsize[i])
            j_range = slice(ind_j, ind_j+icsize[j])
            avg_interSect_cor = np.mean(sectMat[i_range, j_range])
            avgMat[i_range, j_range] = avg_interSect_cor > th*med
            avgMat[j_range, i_range] = avg_interSect_cor > th*med
            coupled_mat[i, j] = avg_interSect_cor > th*med
            coupled_mat[j, i] = avg_interSect_cor > th*med

            ind_j += icsize[j]

        ind_i += icsize[i]

    # Find the correlations between sectors and don't forget weakly internally
    # coupled sectors as individial sectors
    groups = []
    for i in range(len(icsize)):
        coupled = np.nonzero(coupled_mat[i])[0].tolist()
        if any(x in group for x in coupled for group in groups):
            for group in groups:
                if any(x in group for x in coupled):
                    group.extend(list(set(coupled)-set(group)))
        else:
            groups.append(coupled)
        if not coupled:
            groups.append([i])
    groups = [sorted(x) for x in groups if x]

    def filterDuplicates(thelist):
        seen = []
        for x in thelist:
            if x in seen:
                thelist.remove(x)
            seen.append(x)
        return thelist

    groups = filterDuplicates(groups)

    size = []
    for group in groups:
        size.append(sum([icsize[el] for el in group]))
    totsize = sum(icsize)
    biggest_smaller_than_half = max(size) < 0.5*totsize

    return [coupled_mat, groups, biggest_smaller_than_half]


def obtain_positional_frequencies_for_residue(alg, refseqno, refpos):
    # Obtian the amino acid frequencies at all positions for sequences in 
    # alignment alg that have the same amino acid at position refpos as
    # the reference sequence refseq in alignment at position refseqno.

    refseq = alg[refseqno]
    totseq = len(alg)
    refAA = refseq[refpos]
    code = "ACDEFGHIKLMNPQRSTVWY-"
    freq0 = [0.073,
             0.025,
             0.050,
             0.061,
             0.042,
             0.072,
             0.023,
             0.053,
             0.064,
             0.089, 
             0.023,
             0.043,
             0.052,
             0.040,
             0.052,
             0.073,
             0.056,
             0.063,
             0.013, 
             0.033,
             0]  # The last one is the frequency of gaps!

    freqs = []
    freqs_background_cor = []
    for pos in range(len(refseq)):
        AAfreq = dict(zip(code,[0 for i in range(len(code))]))
        AAfreq_bg = dict(zip(code, [-i for i in freq0]))
        for seqno in range(totseq):
            AAfreq[alg[seqno][pos]] += 1/totseq
            AAfreq_bg[alg[seqno][pos]] += 1/totseq
        freqs.append(AAfreq)
        freqs_background_cor.append(AAfreq_bg)

    return freqs, freqs_background_cor