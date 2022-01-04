# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 10:14:55 2021

@author: mmabesoone
"""
from Bio import SeqIO
import os
import sys
from collections import defaultdict
import multiprocessing as mp

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
    try:
        record = SeqIO.read(file, 'genbank')
    except Exception as e:
        sys.stdout.write('\tERROR for {file}: {e}\n')
        return dict()
    grouped_asdomains = defaultdict(list)
    asdomains = [feature for feature in record.features if feature.type=="aSDomain"]
    cdss = [feature for feature in record.features if feature.type=="CDS"]
    pairs = [(cds, asdomain) for asdomain in asdomains for cds in cdss if asdomain.location.start in cds]
    for k,v in pairs:
        grouped_asdomains[k].append(v)

    return grouped_asdomains


def confirm_in_house_file(file):
    if 'Run date     :: 2021' in open(file, 'r').read():
        return True
    else:
        return False

def map_inhouse_to_aSdb(inhouse_file, aSdb_files):
    mapping_to_aSdb = defaultdict(list)
    sys.stdout.write(f"Parsing {os.path.splitext(os.path.basename(inhouse_file))[0]}\n")
    inhouse_seq = SeqIO.read(inhouse_file, 'genbank').seq
    for aSdb_file in aSdb_files:
        aSdb_seq = SeqIO.read(aSdb_file, 'genbank').seq
        if inhouse_seq in aSdb_seq or aSdb_seq in inhouse_seq:
            mapping_to_aSdb[inhouse_file].append(aSdb_file)
            sys.stdout.write(f"\tFound for {os.path.splitext(os.path.basename(inhouse_file))[0]}\n")
    return mapping_to_aSdb


def main(export=False):
    pool = mp.Pool(33)
    folder = 'Z:\Scripts\_antiSMASH_genbanks'
    folder = 'Z:\Scripts\S-004_GenBankDownloader\GenBanks_cisAT_KS_type1'
    files = [os.path.join(folder, file) for file in os.listdir(folder) if '.gb' in file][0:200]

    inhouse_files = [file for file in files if 'region' not in file]
    aSdb_files = [file for file in files if file not in inhouse_files]

    mapping_to_aSdb = pool.starmap(map_inhouse_to_aSdb, [[inhouse_file, aSdb_files] for inhouse_file in inhouse_files])
    sys.stdout.write(str(mapping_to_aSdb))
    fID = open(os.path.join(os.getcwd(),'SequenceMappingToaSDatabase.txt'), 'w')
    fID.write('In-house file\tantiSMASH database file\n')
    for inhouse_file in mapping_to_aSdb:
            fID.write(inhouse_file + '\t' + '\t'.join(mapping_to_aSdb[inhouse_file]))
    fID.close()
    files = aSdb_files
    motifs = [['PKS_KS'], ['PKS_KS', 'Trans-AT_docking'], ['PKS_KS', 'Trans-AT_docking', 'PKS_KR'], ['PP-binding', 'PKS_KS'], ['PKS_KR','PP-binding', 'PKS_KS']]
    motifs = [['ACP','PKS_KS', 'PKS_AT', 'ACP']]
    hits = dict()
    for file in files:
        domain_sequences = extract_domain_sequence(obtain_asdomains(file))
        hits[os.path.splitext(os.path.basename(file))[0]] = count_motif_in_domain_sequence(domain_sequences, motifs)
    if export:
        fID = open(os.path.join(os.getcwd(),"DomainStatistics.txt"), 'w')
        fID.write('File\t' + '\t'.join(['-'.join([domain for domain in motif]) for motif in motifs]) + '\n')
        for file in hits:
            fID.write(file + '\t' + '\t'.join([str(hit) for hit in hits[file]])+ '\n')
        fID.close()
    return None

if __name__ == '__main__':
    main()