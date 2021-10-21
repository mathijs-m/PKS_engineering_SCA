# -*- coding: utf-8 -*-
"""
A script to extract concatenated AT-KS sequences

@author: mmabesoone
"""
from Bio import SeqIO
import os
import sys

def obtain_cds_and_asdomains(seqrecord):
    # Obtain all cdss and asdomains in a list of tuples from the seqrecord
    cdss = [feature for feature in seqrecord.features if feature.type == 'CDS']
    asdomains = [feature for feature in seqrecord.features if feature.type == 'aSDomain']
    cds_asdomains = [(cds, asdomain) for asdomain in asdomains for cds in cdss if asdomain.location.start.real in cds]
    return cds_asdomains

def find_at(cds_asdomains):
    ats = [cds_asdomain[1] for cds_asdomain in cds_asdomains if cds_asdomain[1].qualifiers['aSDomain'][0]=='PKS_AT']
    return ats

def find_kss(cds_asdomains):
    kss = [cds_asdomain[1] for cds_asdomain in cds_asdomains if cds_asdomain[1].qualifiers['aSDomain'][0]=='PKS_KS']
    return kss

def find_ks_adapter_domains(cds_asdomains):
    kss_adapters = []
    for i in range(len(cds_asdomains)-1):
        if cds_asdomains[i][1].qualifiers['aSDomain'] == ['PKS_KS'] and cds_asdomains[i+1][1].qualifiers['aSDomain'] == ['Trans-AT_docking']:
            kss_adapters.append((cds_asdomains[i][1], cds_asdomains[i+1][1]))
        elif cds_asdomains[i][1].qualifiers['aSDomain'] == ['PKS_KS'] and cds_asdomains[i-1][1].qualifiers['aSDomain'] == ['Trans-AT_docking']:
            kss_adapters.append((cds_asdomains[i][1], cds_asdomains[i-1][1]))
    return kss_adapters

def name_and_concatenate_sequences(file, ats, kss):
    concat_sequences = []
    for at in ats:
        for ks in kss:
            if type(ks) == tuple:
                header = '>' + os.path.basename(os.path.splitext(file)[0]) + '_' + at.qualifiers['locus_tag'][0] + '_' + '_'.join([domain.qualifiers['locus_tag'][0] for domain in ks]) + '\n'
                concat_sequences.append(header + at.qualifiers['translation'][0]+''.join(([domain.qualifiers['translation'][0] for domain in ks])))
            else:
                header = '>' + os.path.basename(os.path.splitext(file)[0]) + '_' + at.qualifiers['locus_tag'][0] + '_' + ks.qualifiers['locus_tag'][0] + '\n'
                concat_sequences.append(header + at.qualifiers['translation'][0]+ks.qualifiers['translation'][0])
    return concat_sequences

def save_sequences(folder, fname, concat_sequences):
    destination = os.path.join(folder, fname)
    with open(destination,'w') as fID:
        fID.write('\n'.join([seq for seq in concat_sequences]))

def main():
    genbank_folder = 'Z:\Scripts\_antiSMASH_genbanks'
    dest_folder = 'Z:\Scripts\S-005_SCA_Pipeline\AT_KS_sequences'
    genbanks = [os.path.join(genbank_folder, file) for file in os.listdir(genbank_folder) if '.gb' in file]
    ks_sequences = []
    ks_adapter_sequences = []
    for genbank in genbanks:
        seqrecord = SeqIO.read(genbank, 'genbank')
        cds_asdomains = obtain_cds_and_asdomains(seqrecord)
        ats = find_at(cds_asdomains)
        if ats == []:
            continue
        kss = find_kss(cds_asdomains)
        ks_adapters = find_ks_adapter_domains(cds_asdomains)
        ks_sequences.extend(name_and_concatenate_sequences(genbank, ats, kss))
        ks_adapter_sequences.extend(name_and_concatenate_sequences(genbank, ats, ks_adapters))
    save_sequences(dest_folder, 'AT_KS.fasta', ks_sequences)
    save_sequences(dest_folder, 'AT_KS_adapter.fasta', ks_adapter_sequences)

if __name__ == '__main__':
    main()
