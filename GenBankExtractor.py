# -*- coding: utf-8 -*-
"""
For our SCA analysis on various didomains of trans-AT PKSs, we need to make
allignments of those. For this, we need to be able to easily extract
combinations of KS-XX domains from the MiBiG database.
The MiBiG database is available for downloading

Mathijs Mabesoone, 12-02-2021
"""

import os
from Bio import SeqIO
# import pandas as pd
from shutil import copy
from shutil import rmtree
# from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
# from tkinter import Tk
from tkinter.filedialog import askopenfile
import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
# from time import sleep
from tqdm import tqdm
import pickle
import traceback


class Gene:
    def __init__(self, gene, featNum, metadata):
        self.start = gene.location.start.real
        self.end = gene.location.end.real
        if 'gene' in gene.qualifiers.keys():
            self.name = gene.qualifiers['gene'][0]
        elif 'locus_tag' in gene.qualifiers.keys():
            self.name = gene.qualifiers['locus_tag'][0]
        else:
            self.name = ''
        self.featNum = featNum
        self.direction = gene.location.strand
        self.BGC = metadata[0]
        self.organism = metadata[1]
        # This can be extended with the specificity of the KS later

    def find_domains(self, features):
        # Extracts the domains in the gene from the GenBank data
        # Returns a

        antiSMASH = features[self.featNum].type == 'CDS'

        featNum = self.featNum + 1
        domains = []
        totFeat = len(features)

        while featNum <= totFeat-1 and end(features[featNum]) <= self.end:
            if antiSMASH and features[featNum].type == 'aSDomain':
                domains.append(Domain(features[featNum], featNum, antiSMASH))
            if features[featNum].type == 'CDS' and \
                    'label' in features[featNum].qualifiers.keys():
                domains.append(Domain(features[featNum], featNum, antiSMASH))
            featNum += 1

        if self.direction == -1:
            domains_rev = [domains[-i-1] for i in range(len(domains))]
            domains = domains_rev

        self.domains = domains

    def annotate_domain_sequence(self):
        # Extracts the strand direction corrected domain sequence from the gene
        # Returns a list with the domain sequences

        domain_sequence = []

        for index in range(len(self.domains)):
            domain_sequence.append(self.domains[index].type)
        self.domain_sequence = domain_sequence

    def contains_domain_motif(self, motif):
        # Checks if a certain domain motiv exists in the gene
        # NOTE: the function removes the KS follower number to check
        # Returns True or False and the positions in the domain sequence where
        # the domain motif starts

        try:
            seq = self.domain_sequence
        except AttributeError:
            raise AttributeError('The domain sequence has not yet been \
                                 obtained for ' + self.name)

        for ind, domain in enumerate(seq):
            if domain[0:4] == 'PKS_':
                domain = domain[4:]

            if 'KS' in domain:
                if 'KS0' in domain:
                    pass
                else:
                    seq[ind] = ''.join(list(filter(
                                            lambda x: x.isalpha(), domain)))

        [positions, acps] = find_motif(seq, motif)

        if len(positions) == 0:
            return [False, [], []]
        else:
            return [True, positions, acps]


def find_motif(sequence, motif):
            # Imagine ABCDEAZC
            # Query motif AxC
            motif_positions = []
            acps=[]
            for i in range(0, len(sequence)-len(motif)+1):
                j = 0
                acp = 0
                while motif[j] == sequence[i+j+acp] or \
                        motif[j] == '*WILDCARD*' or \
                            all([motif[j] == '*ACP_n>0*', sequence[i+j+acp] == 'ACP']):
                    if motif[j] == '*ACP_n>0*':
                        while i+j+acp+1 < len(sequence) and \
                                sequence[i+1+j+acp] == 'ACP':
                            acp += 1
                    j += 1
                    if j >= len(motif) or i+j+acp > len(sequence)-1:
                        break

                if j == len(motif):
                    motif_positions.append(i)
                    acps.append(acp)

                i += 1
            return [motif_positions, acps]

class Domain:
    def __init__(self, domain, featNum, antiSMASH):
        self.start = domain.location.start.real
        self.end = domain.location.end.real
        if antiSMASH:
            self.type = domain.qualifiers['aSDomain'][0]
        else:
            self.type = domain.qualifiers['label'][0]
        if self.type[0:4] == 'PKS_':
            self.type = self.type[4:]
        if self.type == 'AT docking' or self.type == 'transAT docking' or \
                self.type == 'Trans-AT_docking':
            self.type = 'trans-AT docking'
        if self.type == 'TE':
            self.type = 'Thioesterase'
        if self.type == 'ACP_beta' or self.type == 'ACP beta':
            self.type = 'ACP'
        self.type = self.type.replace('-binding', ' binding')
        self.type = self.type.replace('_domain', ' domain')

        if 'specificity' in domain.qualifiers.keys():
            specificity = domain.qualifiers['specificity'][0]
            if specificity[0] == '0':
                specificity = specificity[1:]
            while specificity and specificity[-1] in ['.', ' ']:
                specificity = specificity[0:-1]
            if len(specificity) > 0:
                self.specificity = specificity.strip()

        self.featNum = featNum
        self.direction = domain.location.strand
        # This can be extended with the specificity of the KS later


def new_folder(folderName):
    if not os.path.isdir(folderName):
        os.mkdir(folderName)
    return folderName


def clear_folder(folder):
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
    try:
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            rmtree(file_path)
    except Exception as e:
        ('Failed to delete %s.  Reason: %s' % (file_path, e))


def separate_BGCs(GenBankFile, GenBankDir):
    # Separates the individual GenBank files in the long list of GenBank files

    GenBankFile = os.path.join(GenBankFile)

    folder_indivBGCs = 'BGC Genbanks'
    if os.path.isdir(folder_indivBGCs):
        rmtree(folder_indivBGCs)

    folderName = new_folder('/'.join([os.getcwd(), GenBankDir]))

    lines = open(GenBankFile).readlines()
    lastLine = 0
    numLines = len(lines)

    while lastLine+1 < numLines:
        firstLine = lastLine
        # Find first entry
        while 'LOCUS' not in lines[firstLine]:
            firstLine += 1

        locus = lines[firstLine].split(' ')
        while locus[1] == '':
            del locus[1]
        NatProd = locus[1]
        NatProd = NatProd.replace('/', '-')

        # Find next entry and save
        lastLine = firstLine
        while 'LOCUS' not in lines[lastLine+1]:
            lastLine += 1
            if lastLine + 1 == numLines:
                break

        try:
            fileID = open('/'.join([folderName, NatProd + '.gb']), 'w+')
            fileID.write(''.join(lines[firstLine:lastLine+1]))
            fileID.close()
        except FileNotFoundError:
            raise ValueError(NatProd + ' gave an error in line ' +
                             str(firstLine))
        except PermissionError:
            print('Denied writing permission', 'Writing permission denied on ' +
                       str(os.path.join(os.getcwd(), folder_indivBGCs)) +
                       '. \n Please copy the executable to another folder.')


def obtainSourceData(GenBank):
    # Try to extract the source organism from a parsed GenBank file
    # Returns the organism name, country and db_xref or 'NotFound'
    # !The organism that is annotated as feature precedes over the 'organism'
    # in the header!

    organism = 'organismNotFound'
    country = 'countryNotFound'
    db_xref = 'xrefNotFound'
    accession = 'accessionNotFound'

    if 'accessions' in GenBank.annotations.keys():
        accession = GenBank.annotations['accessions'][0]
    if 'organism' in GenBank.annotations.keys():
        if len(GenBank.annotations['organism']) == 0:
            pass
        else:
            organism = GenBank.annotations['organism']

    for feature in GenBank.features:
        if feature.type == 'source':
            keys = feature.qualifiers.keys()
            if 'organism' in keys:
                organism = feature.qualifiers['organism'][0]
                #if 'strain' in keys:
                #    organism = organism + ', ' + \
                #        str(feature.qualifiers['strain'][0])
            if 'country' in keys:
                country = feature.qualifiers['country'][0]
            if 'db_xref' in keys:
                db_xref = feature.qualifiers['db_xref'][0]

    return [organism, country, db_xref, accession]


def find_next_gene(gb, featNum, totFeat, metadata):
    # Finds the name and feature number of the next gene in the gb file

    # adapt for antiSMASH
    try:
        if 'antiSMASH-Data' in gb.annotations['structured_comment'].keys():
            antiSMASH = True
        else:
            antiSMASH = False
    except KeyError:
        antiSMASH = False

    if antiSMASH:
        while gb.features[featNum].type != 'CDS':
            featNum += 1
            if featNum == totFeat-1:
                break

        if featNum == totFeat-1 and gb.features[featNum].type != 'gene':
            nextGene = None
        else:
            nextGene = Gene(gb.features[featNum], featNum, metadata)

    else:
        while gb.features[featNum].type != 'gene':
            featNum += 1
            if featNum == totFeat-1:
                break

        if featNum == totFeat-1 and gb.features[featNum].type != 'gene':
            nextGene = None
        else:
            nextGene = Gene(gb.features[featNum], featNum, metadata)

    return nextGene


def new_domains(domain_types, gene):
    # This function checks if there are new domain types in the gene
    # Returns the new domains

    new_domains = []
    for domain in gene.domain_sequence:
        if domain[0:4] == 'PKS_':
            domain = domain[4:]
        if 'KS' in domain:
            if 'KS0' in domain:
                pass
            else:
                domain = ''.join(list(filter(
                    lambda x: x.isalpha(), domain)))

        if domain not in domain_types and domain not in new_domains:
            new_domains.append(domain)

    return new_domains


def new_specificity(specificities, gene):
    # This function checks if there are new domain types in the gene
    # Returns the new domains

    new_specificity = []
    for domain in gene.domains:
        try:
            specificity = domain.specificity
            if specificity not in specificities and \
                    specificity not in new_specificity:
                new_specificity.append(specificity)
        except Exception:
            pass

    return new_specificity


def end(feature):
    # Returns the end of a GenBank feature
    return feature.location.end.real


def extract_domains(gb, featNumGene):
    # This function extract the domain sequence from an annotated gene in the
    # GenBank file

    featNum = featNumGene + 1
    domains = dict()
    domainNum = 0

    if 'antiSMASH-Data' in gb.annotations['structured_comment'].keys():
        antiSMASH = True
    else:
        antiSMASH = False

    # adapt for antiSMASH
    while end(gb, featNumGene) >= end(gb, featNum):
        if antiSMASH and gb.features[featNum].type == 'aSDomain':
            domains[domainNum] = Domain(gb.features[featNum], featNum,
                                        antiSMASH)
            domainNum += 1
        elif not antiSMASH and gb.features[featNum].type == 'CDS' and \
                'label' in gb.features[featNum].qualifiers.keys():
            domains[domainNum] = Domain(gb.features[featNum], featNum,
                                        antiSMASH)
            domainNum += 1
        featNum += 1

    return domains


def motif_sequence(gene, pos, gb, basepairs, only_annotations):
    # This function extracht the DNA and protein sequence of the motif in the
    # gene. If the gene is long enough, the function takes 300 bases or 100
    # amino acids additionally on both the beginning and end of the sequence
    # The basepairs list containst [leader, follower, intermediate]
    # If only_annotations = true, only the annoations with spacings with the 
    # length of 'intermediate' bp at either side of non-peripheral domains will
    # be extracted.
    # E.g. basepairs = [100, 100, 50] results in 100-dom1-50-50-dom2-50-...
    # 50-dom_last-100.
    # If the domains are close together than 2x intermediate, the whole linking
    # regions is taken, without overlaps.
    # Returns both the DNA and protein sequence

    if gene.direction == 1:
        start = pos[0]
        end = pos[1]
        genestart = gene.start
        genestop = gene.end
    else:
        start = pos[0]
        end = pos[1]
        genestart = gene.end
        genestop = gene.start

    leader = basepairs[0]
    follower = basepairs[1]
    intermediate = 0
    if len(basepairs) == 3:
        intermediate = basepairs[2]

    if only_annotations:
        domain_positions = range(start, end+gene.direction, gene.direction)
        DNAseq = ''
        AAseq = ''
        for position in domain_positions:
            # Only add the leader and followers at the periphery of the domain
            # motif. In between, only extract the annotations.
            # If first domain, take the leading spacing in account
            if position == domain_positions[0]:
                if gene.direction == 1:
                    startpos = max(gene.domains[position].start -
                                   gene.direction*leader, genestart)
                    previous_end = 0
                else:
                    startpos = min(gene.domains[position].end -
                                   gene.direction*leader, genestart)
                    previous_end = gene.end
            else:
                startpos = gene.domains[position].start - \
                    gene.direction*intermediate
                previous_end = endpos
            while (gene.domains[position].start - startpos) % 3 != 0:
                startpos += gene.direction*1
            if (startpos - previous_end)*gene.direction < 0:
                startpos = previous_end

            # If last domain, take following spacing in account
            if position == domain_positions[-1]:
                if gene.direction == 1:
                    endpos = min(gene.domains[position].end +
                                 gene.direction*follower, len(gb), genestop)
                else:
                    endpos = max(gene.domains[position].start +
                                 gene.direction*follower, 0, genestop)

            else:
                endpos = gene.domains[position].end + \
                         gene.direction*intermediate

            while (endpos - gene.domains[end].end) % 3 != 0:
                endpos -= gene.direction*1

            if gene.direction == 1:
                dnaseq = gb[startpos:endpos].seq
            if gene.direction == -1:
                dnaseq = gb[endpos:startpos].seq
                dnaseq = dnaseq.reverse_complement()
            if len(dnaseq) % 3 != 0:
                print('The sequence of ' + gb.name + ', ' + gene.name +
                                 ', domain ' + str(pos[0]) + 
                                 ' is not a multiple of 3. Length is: ' +
                                 str(len(DNAseq)) + '. start|end: ' + str(startpos) + ', ' + str(endpos)+
                                 str(DNAseq) + '\nI will proceed anyway, but check the database!')
            aaseq = str(dnaseq.translate())
            dnaseq = str(dnaseq)
            DNAseq += dnaseq
            AAseq += aaseq

        if len(AAseq) > 0 and AAseq[-1] == '*':
            AAseq = AAseq[0:-1]
            print('I detected a stop codon for ' + gb.name + ', ' + gene.name +
                  ', domain ' + str(pos[0]) + '. Removing it from the AA code, but keeping it in the DNA.')

    else:
        if gene.direction == 1:
            startpos = max(gene.domains[start].start -
                           gene.direction*leader, genestart)
        else:
            startpos = min(gene.domains[start].end -
                           gene.direction*leader, genestart)
        while (gene.domains[start].start - startpos) % 3 != 0:
            startpos += 1*gene.direction

        if gene.direction == 1:
            endpos = min(gene.domains[end].end +
                         gene.direction*follower, len(gb), genestop)
        else:
            endpos = max(gene.domains[end].start +
                         gene.direction*follower, 0, genestop)

        while (endpos - gene.domains[end].end) % 3 != 0:
            endpos -= 1*gene.direction

        if gene.direction == 1:
            motif_slice = gb[startpos:endpos]
        if gene.direction == -1:
            motif_slice = gb[endpos:startpos]

        DNAseq = motif_slice.seq
        if gene.direction is -1:
            DNAseq = DNAseq.reverse_complement()

        if len(DNAseq) % 3 != 0:
            print('The sequence of ' + gb.name + ', ' + gene.name +
                             ', domain ' + str(pos[0]) + 
                             ' is not a multiple of 3. Length is: ' +
                             str(len(DNAseq)) + '. start|end: ' + str(startpos) + ', ' + str(endpos)+
                             str(DNAseq) + '\nI will proceed anyway, but check the database!')

        AAseq = str(DNAseq.translate())
        if AAseq[-1] == '*':
            AAseq = AAseq[0:-1]
            print('I detected a stop codon for ' + gb.name + ', ' + gene.name +
                  ', domain ' + str(pos[0]) + '. Removing it from the AA code, but keeping it in the DNA.')
        DNAseq = str(DNAseq)

    return [DNAseq, AAseq]


def scan_genbanks(genbank_folder):
    # This function scans all the gen bank files in genbank_folder and returns
    # a dictionary containing all module types, BGC names and organisms

    props = dict()
    keys = ['organisms', 'BGCs', 'domain_types', 'specificities']
    for key in keys:
        props[key] = []

    props['domain_types'].extend(['*WILDCARD*','*ACP_n>0*'])
    props['specificities'].extend(['NotAnnotated'])

    genbanks = os.listdir(genbank_folder)
    for genbank in genbanks:
        try:
            if '.gb' not in genbank:
                continue
            gb = SeqIO.read('/'.join([genbank_folder, genbank]), 'genbank')
        except ValueError:
            continue

        BGCname = gb.name

        [organism, country, db_xref, accession] = obtainSourceData(gb)

        metadata = [BGCname, organism]
        totFeat = len(gb.features)
        featNum_nextGene = 0
        genes = []

        # Find all genes in the BGC and find all domains in the genes
        while featNum_nextGene < totFeat-1:
            gene = find_next_gene(gb, featNum_nextGene, totFeat, metadata)
            if gene is None:
                break

            gene.find_domains(gb.features)
            gene.annotate_domain_sequence()

            props['domain_types'].extend(new_domains(
                props['domain_types'], gene))
            props['specificities'].extend(new_specificity(
                props['specificities'], gene))
            num_domains = len(gene.domains)
            if num_domains == 0:
                featNum_nextGene = gene.featNum + 1
            else:
                featNum_nextGene = gene.domains[
                    len(gene.domains)-1].featNum + 1

        if organism not in props['organisms']:
            props['organisms'].append(organism)
        if BGCname not in props['BGCs']:
            props['BGCs'].append(BGCname)

        genes.append(gene)

    for key in props.keys():
        props[key] = sorted(props[key])

    return props


def load_database_properties(indiv_genbank_dir):
    # This function loads the database properties
    # Returns a dictionary with the properties

    try:
        if not os.path.exists(indiv_genbank_dir):
            GenBankFile = askopenfile(title='Select the location of the '
                                      + 'GenBank database')
            if GenBankFile == 'None':
                return
            GenBankFile = GenBankFile.name
            while GenBankFile[-3:] != '.gb':
                GenBankFile = askopenfile(title='Please select a .gb file')
            if GenBankFile is None:
                return
            separate_BGCs(GenBankFile, indiv_genbank_dir)
    except PermissionError:
        print('Denied writing permission', 'Writing permission denied on ' +
                   str(os.path.join(os.getcwd(), indiv_genbank_dir)) +
                   '. \n Please copy the executable to another folder.')
    else:
        GenBankFile = None

    try:
        with open('/'.join([os.getcwd(), '_Settings', 'settings.pkl']), 'rb')\
                as f:
            properties = pickle.load(f)

        if len(properties) < 5:
            raise IOError('Too few properties stored in the settings file.')

    except IOError:
        properties = scan_genbanks(indiv_genbank_dir)
        properties['genbank_dir'] = indiv_genbank_dir
        properties['database_path'] = GenBankFile

        settingsFolder = new_folder('_Settings')
        f = open('/'.join(['_Settings', 'settings.pkl']), 'wb')
        pickle.dump(properties, f, pickle.HIGHEST_PROTOCOL)

    return properties


def load_database_properties_api(indiv_genbank_dir):
    # This function loads the database properties
    # Returns a dictionary with the properties

    try:
        with open('/'.join([os.getcwd(), 'Settings', 'settings.pkl']), 'rb')\
                as f:
            properties = pickle.load(f)

        if len(properties) < 5:
            raise IOError('Too few properties stored in the settings file.')

    except IOError:
        properties = scan_genbanks(indiv_genbank_dir)
        properties['genbank_dir'] = indiv_genbank_dir
        try:
            properties['database_path'] = None
        except Exception:
            properties['database_path'] = None

        settingsFolder = new_folder('_Settings')
        f = open('/'.join(['_Settings', 'settings.pkl']), 'wb')
        pickle.dump(properties, f, pickle.HIGHEST_PROTOCOL)

    return properties


def reload_database_properties(indiv_genbank_dir):
    # This function loads the database properties
    # Returns a dictionary with the properties

    if os.listdir(indiv_genbank_dir):
        use_old_gbs = messagebox.askquestion('Found existing genbanks',
                           'I found already existing individual GenBank files in the GenBank directory. \n Do you want to use these existing individual GenBank files instead of a new database?')
        if use_old_gbs == 'yes':
            GenBankFile = 'None'
        elif use_old_gbs == 'no':
            GenBankFile = askopenfile(title=
                                      'Select the location of the GenBank' +
                                      ' database')
            if GenBankFile is None:
                return
            GenBankFile = GenBankFile.name
            while GenBankFile[-3:] != '.gb':
                GenBankFile = askopenfile('Please select a .gb file')
            clear_folder(indiv_genbank_dir)
            separate_BGCs(GenBankFile, indiv_genbank_dir)
    else:
        GenBankFile = askopenfile(title='Select the location of the GenBank' +
                                  ' database')
        if GenBankFile is None:
            return
        GenBankFile = GenBankFile.name
        while GenBankFile[-3:] != '.gb':
            GenBankFile = askopenfile('Please select a .gb file')
        separate_BGCs(GenBankFile, indiv_genbank_dir)
    properties = scan_genbanks(indiv_genbank_dir)
    properties['genbank_dir'] = indiv_genbank_dir
    properties['database_path'] = GenBankFile

    settingsFolder = new_folder('Settings')
    with open('/'.join([settingsFolder, 'settings.pkl']), 'wb') as f:
        pickle.dump(properties, f, pickle.HIGHEST_PROTOCOL)

    return properties


def obtain_sequence_and_header(gb, gene, query, start_positions, acps,
                               basepairs, metadata, specificity_filter,
                               only_annotations):
    # This function extracts the AA and DNA sequence and constructs the header
    # When the specificities are not in the domains, the funtion returns None

    [BGCname, organism, country, db_xref, accession] = metadata

    len_motif = len(query)

    DNAseqs = []
    AAseqs = []
    headers = []

    for i, start in enumerate(start_positions):
        specificities = []
        for domain in gene.domains[start:start+len_motif+acps[i]]:
            if domain.type[0:2] == 'KS':
                if 'specificity' in dir(domain):
                    spec = domain.specificity
                    specificities.append(spec)

                else:
                    spec = 'NotAnnotated'
                    specificities.append(spec)

        if not all([i in specificity_filter for i in specificities]):
            continue

        [DNAseq, AAseq] = \
            motif_sequence(gene,
                           [start, start+len_motif-1], gb,
                           basepairs, only_annotations)
        header = '>' + ' | '.join([BGCname, 'gene: ' + gene.name,
                                   'start, end in domain sequence: ['
                                   + str(start+1) + ', ' + str(start+len_motif)
                                   + ']', organism, country, db_xref,
                                   accession, ', '.join(specificities)]) + '\n'
        DNAseqs.append(DNAseq)
        AAseqs.append(AAseq)
        headers.append(header)

    return [DNAseqs, AAseqs, headers]


def main_extractor(GenBankDir, query, basepairs, BGC_filter, organism_filter,
                   specificity_filter, properties, fname, dest_folder,
                   only_annotations):
    # This function does the main extraction based on the limitations supplied
    # Writes the files and exits

    seqCounter = 0
    genes = dict()
    headers = []
    DNAseqs = []
    AAseqs = []

    # Analyze all GenBank files of all relevant BGCs for the presence of the
    # motif
    # Store all data in the headers and seqs lists
    for BGC in os.listdir(GenBankDir):
        if '.gb' not in BGC:
            continue

        gb = SeqIO.read('/'.join([GenBankDir, BGC]), 'genbank')

        try:
            if 'antiSMASH-Data' in gb.annotations['structured_comment'].keys():
                antiSMASH = True
            else:
                antiSMASH = False
        except KeyError:
            antiSMASH = False

        BGCname = gb.name
        if BGCname not in BGC_filter:
            continue
        if 'region' in BGC:
            BGCname += BGC.split('.')[-2]

        [organism, country, db_xref, accession] = obtainSourceData(gb)
        if organism not in organism_filter:
            continue

        metadata = [BGCname, organism, country, db_xref, accession]

        totFeat = len(gb.features)
        featNum_nextGene = 0
        genes[BGCname] = []

        # Find all genes in the BGC and find all domains in the genes
        while featNum_nextGene < totFeat-1:
            gene = find_next_gene(gb, featNum_nextGene, totFeat, metadata)

            if gene is None:
                break
            try:
                gene.find_domains(gb.features)
                gene.annotate_domain_sequence()

                [motif_in_gene, start_positions, acps] = \
                    gene.contains_domain_motif(query)

                if motif_in_gene:
                    results = obtain_sequence_and_header(gb, gene, query,
                                                         start_positions, acps,
                                                         basepairs, metadata,
                                                         specificity_filter,
                                                         only_annotations)
                    DNAseq_list = results[0]
                    AAseq_list = results[1]
                    header_list = results[2]
                    if DNAseq_list:
                        headers.extend([h for h in header_list])
                        DNAseqs.extend([seq for seq in DNAseq_list])
                        AAseqs.extend([seq for seq in AAseq_list])

                # Store the gene and move on to the next one
                genes[BGCname].append(gene)
                num_domains = len(gene.domains)
            except Exception as e:
                traceback.print_exc()
                print(e)
                # print(featNumGene)
                print(gene)
                print('Encountered an error in ' + gene.name + ' of the ' + BGCname + ' cluster\nThe error is:')
                print('Proceeding nonetheless...\n')

            if num_domains == 0:
                featNum_nextGene = gene.featNum + 1
            else:
                featNum_nextGene = gene.domains[len(gene.domains)-1].featNum+1

    DNAfold = new_folder(dest_folder)
    AAfold = new_folder(dest_folder)

    if not fname:
        fname = '_'.join(query)
    fname = ''.join(x for x in fname if x.isalnum() or x in ['_', '-', ','])
    print('filename is ' + fname)

    try:
        if 'DNAfile' not in locals():
            DNAfile = open(DNAfold + '/' + fname +
                           '_DNA.txt', 'w+')

        if 'AAfile' not in locals():
            AAfile = open(AAfold + '/' + fname +
                          '_AA.txt', 'w+')

        numSeqs = len(DNAseqs)

        for index in range(len(DNAseqs)):
            if len(DNAseqs[index]) < 3:
                continue
            DNAfile.write(headers[index])
            DNAfile.write(DNAseqs[index]+'\n')
            AAfile.write(headers[index])
            AAfile.write(str(AAseqs[index])+'\n')

        DNAfile.close()
        AAfile.close()
    except PermissionError:
        print('Denied writing permission', 'Writing permission denied.' +
           '. \n Please copy the executable to another folder.')

    return [numSeqs, headers, DNAseqs, AAseqs, genes]
