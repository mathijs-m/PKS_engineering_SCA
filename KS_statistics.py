# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 10:14:55 2021

@author: mmabesoone
"""
from Bio import SeqIO
import os
from collections import defaultdict

def count_motif_in_domain_sequence(domain_sequences, motifs):
    motif_hits = list()
    for motif in motifs:
        i = 0
        hits = 0
        for cds in domain_sequences:
            for domain in cds:
                if i == len(motif):
                    hits += 1
                    i = 0
                if motif[i] == domain:
                    i += 1
                else:
                    i = 0
        motif_hits.append(hits)
    return motif_hits


def extract_domain_sequence(grouped_asdomains):
    domain_sequence = list()
    for cds in grouped_asdomains:
        cds_sequence = list()
        for asdomain in grouped_asdomains[cds]:
            cds_sequence.append(asdomain.qualifiers['aSDomain'][0])
        if cds_sequence != []:
            domain_sequence.append(cds_sequence)
    return domain_sequence


def obtain_asdomains(file):
    record = SeqIO.read(file, 'genbank')
    grouped_asdomains = defaultdict(list)
    asdomains = [feature for feature in record.features if feature.type=="aSDomain"]
    cdss = [feature for feature in record.features if feature.type=="CDS"]
    pairs = [(cds, asdomain) for asdomain in asdomains for cds in cdss if asdomain.location.start in cds]
    for k,v in pairs:
        grouped_asdomains[k].append(v)

    return grouped_asdomains


def main():
    folder = 'Z:\Scripts\_antiSMASH_genbanks'
    files = [os.path.join(folder, file) for file in os.listdir(folder) if '.gb' in file]
    motifs = [['PKS_KS'], ['PKS_KS', 'Trans-AT docking'], ['PKS_KS', 'Trans-AT docking', 'PKS_KR'], ['PP-binding', 'PKS_KS'], ['PKS_KR','PP-binding', 'PKS_KS']]
    hits = dict()
    for file in files:
        domain_sequences = extract_domain_sequence(obtain_asdomains(file))
        hits[os.path.splitext(os.path.basename(file))[0]] = count_motif_in_domain_sequence(domain_sequences, motifs)
    fID = open(os.path.join(os.getcwd(),"DomainStatistics.txt"), 'w')
    fID.write('File\t' + '\t'.join(motifs) + '\n')
    for file in hits:
        fID.write(file + '\t' + str('\t'.join(hits[file])) + '\n')
    fID.close()

if __name__ == '__main__':
    main()