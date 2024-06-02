import pandas as pd
import sys
sys.path.append('.')
import functions.sequencefunctions as sf
from datetime import time

vdj = pd.read_csv('data/tanno_A1naive_with_stitchr_seq.csv.gz', index_col=0).dropna(subset = ['TCRa', 'TCRb'], how='any')
print(vdj.shape)
vdj0 = vdj.drop_duplicates(subset=['TCRa', 'TCRb'])
print('keep only complete sequences: ', vdj.shape)
vdj1 = vdj0.copy()
print('finding cdr3s...')
vdj1['cdr3a_IMGTgaps'] = [x if pd.isna(x) else sf.get_cdr3_seq_with_gaps(x.strip()) for x in vdj1.TCRa]
vdj1['cdr3b_IMGTgaps'] = [x if pd.isna(x) else sf.get_cdr3_seq_with_gaps(x.strip()) for x in vdj1.TCRb]

print('run checks...')
assert (vdj1.dropna(subset=['cdr3a_IMGTgaps']).apply(lambda x: x.cdr3a_IMGTgaps.replace('-', '') in x.TCRa, axis=1)).all()
assert (vdj1.dropna(subset=['cdr3b_IMGTgaps']).apply(lambda x: x.cdr3b_IMGTgaps.replace('-', '') in x.TCRb, axis=1)).all()

vdj2 = vdj1.dropna(subset=['cdr3a_IMGTgaps', 'cdr3b_IMGTgaps'])
# remove cdr3s that do not start with C
vdj3 = vdj2.loc[(vdj2['cdr3a_IMGTgaps'].str[0] == 'C') & (vdj2['cdr3b_IMGTgaps'].str[0] == 'C')]
vdj3['len_cdr3a'] = [len(x) for x in vdj3['cdr3a_IMGTgaps']]
vdj3['len_cdr3b'] = [len(x) for x in vdj3['cdr3b_IMGTgaps']]
# remove sequences that are stupidly long
print(vdj3.len_cdr3a.value_counts())# only 22 with length > 19
print(vdj3.len_cdr3b.value_counts()) # only 124 with length > 19

vdj4 = vdj3.loc[(vdj3.len_cdr3a < 20) & (vdj3.len_cdr3b < 20)]
print('saving without TCR...')
vdj4.to_csv('data/tanno_A1naive_subset_for_MI.csv')

# print('calculate renumbered TCR...')
# vdj4['TCRa_IMGTgaps'] = [sf.get_vregion_seq_with_gaps(x) for x in vdj4.TCRa]
# vdj4['TCRb_IMGTgaps'] = [sf.get_vregion_seq_with_gaps(x) for x in vdj4.TCRb]

# print('saving with TCR...')
# vdj4.to_csv('data/tanno_A1naive_IMGT_TCRs.csv.gz')