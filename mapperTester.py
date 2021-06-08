# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 15:08:10 2021

@author: mmabesoone
"""
for i in range(len(headers_full)):
    if Dseq['hd'][model] == headers_full[i]:
        index = i
        header = headers_full[index]
        break
print(model)
print(index)
print(Dseq['alg'][model])
print('......')
ats_sca = Dseq['ats']
pos_in_sca = 0
pos_in_ats = 0
full_seq_model = sequences_full[index]
print(full_seq_model)
sca_seq = Dseq['alg'][model]
full_seq = ''
new_ats = []
insert_count = 0
for pos_in_seq, let in enumerate(full_seq_model):
    if let == '-':
        insert_count += 1
    else:
        full_seq += let
        if let == sca_seq[pos_in_sca]:
            if pos_in_sca in ats_sca:
                new_ats.append(pos_in_sca + insert_count)
                pos_in_ats += 1
                print(pos_in_sca)
            pos_in_sca += 1

mapping = dict()
mapping['-'] = None
nxt = 0
for i, let in enumerate(full_seq):
    if let == sca_seq[nxt]:
        x = ' - ' + let + ' ' + str(nxt)
        mapping[nxt] = i
        nxt += 1
