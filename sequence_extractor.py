# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 15:35:28 2022

@author: mmabesoone
"""
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import sys
from pathlib import Path
from collections import defaultdict

def convert_AA_range_to_DNA_range(leading, intermediate, trailing):
    return [leading*3, intermediate*3, trailing*3]

def find_asdomains(seqRecord):
    # Finds aSDomains and groups them by CDS
    grouped_asdomains = defaultdict(list)
    asdomains = [feature for feature in seqRecord.features if feature.type=="aSDomain"]
    cdss = [feature for feature in seqRecord.features if feature.type=="CDS"]
    pairs = [(cds, asdomain) for asdomain in asdomains for cds in cdss if asdomain.location.start in cds]
    for k,v in pairs:
        grouped_asdomains[k].append(v)

    return grouped_asdomains
 
def combine_overlapping_sequence_boundaries(boundaries):
    # Sort the boundaries by their starting points
    sorted_boundaries = sorted(boundaries, key=lambda x: x[0])

    merged_boundaries = []
    current_start, current_end = sorted_boundaries[0]

    for start, end in sorted_boundaries[1:]:
        # If the current boundary overlaps with the next boundary, update the current_end
        if current_end >= start:
            current_end = max(current_end, end)
        else:
            # If there is no overlap, add the current merged boundary and update the current_start and current_end
            merged_boundaries.append([current_start, current_end])
            current_start, current_end = start, end

    # Add the last merged boundary
    merged_boundaries.append([current_start, current_end])

    return merged_boundaries


def extract_motif_from_gb(seqRecord, motif, leading, trailing, intermediate, any_binding):
    '''
    Extract the sequences of a motif from a genbank file

    '''
    gene_dict = find_asdomains(seqRecord)

    spacing = [leading, intermediate, trailing]
    domain_hit_sequences = []
    for gene in gene_dict:
        aSDomains = gene_dict[gene]
        if gene.strand == 1:
            motif_iterator = 0
        else:
            motif_iterator = -1
        #motif_subhits = list()
        position_list = list()
        for aSDomain in aSDomains:
            if any_binding and aSDomain.qualifiers['aSDomain'][0].replace('PKS_','') in\
                    ['PP-binding', 'ACP', 'ACP_beta', 'AMP_binding', 'PCP', 'PP']:
                aSDomain.qualifiers['aSDomain'][0] = 'binding'
            if aSDomain.qualifiers['aSDomain'][0].replace('PKS_', '').replace('cMT','MT').replace('oMT','MT').replace('nMT','MT').lower() == motif[motif_iterator].lower():
                if motif_iterator == 0 or motif_iterator == -1:
                    # When the first domain in the motif is found in the PKS...
                    # Find the starting and ending positions of the motif, taking into account the flanking bases and
                    # strand direction
                    start = max(gene.location.start, aSDomain.location.start - spacing[int(0.5-gene.strand/2)])
                    end = min(gene.location.end, aSDomain.location.end + spacing[int(0.5-gene.strand/2)])
                    position_list.append([start, end])
                    #motif_subhits.append(get_sequence_with_surrounding(gene, aSDomain, seqRecord, leading, intermediate))
                    motif_iterator += 1 * gene.strand
                elif motif_iterator == len(motif)*gene.strand - int(0.5+gene.strand/2):
                    # Add the gene.strand/2 to account for the fact that the negative indexing goes from -1 to -len(motif)
                    start = max(gene.location.start, aSDomain.location.start - spacing[int(0.5-gene.strand/2)])
                    end = min(gene.location.end, aSDomain.location.end + spacing[int(0.5-gene.strand/2)])
                    position_list.append([start, end])
                    # Combine any overlapping regions in the position_list
                    position_list = combine_overlapping_sequence_boundaries(position_list)
                    # Get the sequence of the motif
                    AA_sequence = ''
                    for start, end in position_list[::gene.strand]:
                        # Make a SeqFeature object from the start and end positions
                        DNA_sequence = FeatureLocation(start, end).extract(seqRecord)
                        # Reverse complement the sequence if the gene is on the negative strand
                        if gene.strand == -1:
                            DNA_sequence = DNA_sequence.reverse_complement()
                        # Translate the sequence to protein and add it to the list of domain hits
                        AA_sequence += DNA_sequence.translate()
                    domain_hit_sequences.append((int(position_list[0][0]), int(position_list[-1][1]), AA_sequence))
                    position_list = list()
                    motif_iterator = int(gene.strand/2 - gene.strand**2/2)
                else:
                    start = max(gene.location.start, aSDomain.location.start - spacing[int(0.5-gene.strand/2)])
                    end = min(gene.location.end, aSDomain.location.end + spacing[int(0.5-gene.strand/2)])
                    position_list.append([start, end])
                    motif_iterator += 1 * gene.strand
            else:
                motif_iterator = int(gene.strand/2 - gene.strand**2/2)
                position_list = list()

    return domain_hit_sequences


def convert_domain_hits_to_sequences(domain_hits, file):
    '''
    From the list of domain hits, concatenate the sequences into fasta-style text

    Parameters
    ----------
    domain_hits : list of lists of SeqRecords
        A list of lists containing the individual SeqRecords matching the domain hits.

    Returns
    -------
    DNA_fasta : str
        The DNA sequences of the domain hits in fasta format
    AA_fasta : str
        The AA sequences of the domain hits in fasta format
    '''
    AA_fasta = ''
    for domain_hit in domain_hits:
        header = f'>{Path(file).stem} | start/end in genome: {domain_hit[0]} - {domain_hit[1]}'
        AA_fasta += header + '\n' + str(domain_hit[2].seq)+ '\n'
    return AA_fasta

def extract_domain_sequences(motif, genbank_folder, destination_folder, leading=100, trailing=100, intermediate=15, any_binding=False):
    '''
    Extract the sequences of a motif from a genbank file
    '''
    leading, intermediate, trailing = convert_AA_range_to_DNA_range(leading, intermediate, trailing)
    
    # Prepare the motif
    if any_binding:
        motif = [element if element not in ['PP-binding', 'ACP', 'ACP_beta', 'AMP_binding', 'PCP', 'PP'] else 'binding' for element in motif ]
    motif = [element.replace(' ','_') for element in motif]
    
    # Run the extraction
    files = [file for file in Path(genbank_folder).iterdir() if '.gb' in file.name]
    AA_fasta = open(Path(destination_folder).joinpath('_'.join(motif)+'_AA.fasta'), 'wt')
    hits = 0
    for file in files:
        seqRecord = SeqIO.read(file, 'genbank')
        domain_hits = extract_motif_from_gb(seqRecord, motif, leading, trailing, intermediate, any_binding)
        hits += len(domain_hits)
        AA_sequences = convert_domain_hits_to_sequences(domain_hits, file)
        AA_fasta.write(AA_sequences)
    AA_fasta.close()
    
    return [hits]
        