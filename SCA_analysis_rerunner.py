# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 13:29:34 2021

@author: mmabesoone
"""
import os
import sys
sys.path.append('../_analysis')
sys.path.append('../_sca')


motifs = ['ACP_KS', 'ACP_KS_trans-ATdocking', 'KR_ACP_KS', 'KS', 'KS_trans-ATdocking',
          'KS_trans-ATdocking_KR']

folder = 'Z:/Scripts/S-005_SCA_Pipeline/_analysis'
pysca_folder = 'Z:/Scripts/S-005_SCA_Pipeline/_sca/'
os.chdir(folder)


for motif in motifs:
    print('Analyzing '+motif)
    folder2 = 'Z:/Scripts/S-005_SCA_Pipeline/Results/Results_2021-04-21_FullSequence/' + motif
    os.system('python Z:/Scripts/S-005_SCA_Pipeline/_analysis/KSanalysis_forReanalysis.py ' + motif + ' ' + folder2 + ' ' + pysca_folder)

motifs = ['ACP_KS', 'ACP_KS_trans-ATdocking', 'KR_ACP_KS', 'KS', 'KS_trans-ATdocking',
          'KS_trans-ATdocking_KR', 'Condensation', 'KR', 'DH', 'ER', 'AT',]

for motif in motifs:
    print('Analyzing '+motif)
    folder2 = 'Z:/Scripts/S-005_SCA_Pipeline/Results/Results_2021-04-21/' + motif
    os.system('python Z:/Scripts/S-005_SCA_Pipeline/_analysis/KSanalysis_forReanalysis.py ' + motif + ' ' + folder2 + ' ' + pysca_folder)

cis_at_motifs = ['ACP_KS_AT', 'KS', 'KS_AT', 'KS_AT_KR']

for motif in cis_at_motifs:
    print('Analyzing '+motif)
    folder2 = 'Z:/Scripts/S-005_SCA_Pipeline/Results_cis-AT/Results_2021-04-12/' + motif
    os.system('python Z:/Scripts/S-005_SCA_Pipeline/_analysis/KSanalysis_forReanalysis.py ' + motif + ' ' + folder2 + ' ' + pysca_folder)
