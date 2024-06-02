import pandas as pd
import sys
sys.path.append('.')
import functions.sequencefunctions as sf

vdj = pd.read_csv('data/vdj_cleaned_set_with_stitchr_seq.csv', index_col=0).dropna(subset = ['TCRa', 'TCRb'], how='any')
print(vdj.shape)
vdj0 = vdj.drop_duplicates(subset=['TCRa', 'TCRb', 'Epitope'])
print('keep only complete sequences: ', vdj.shape)
counts = vdj0.Epitope.value_counts()
# not_too_big = [x for x in counts.keys() if counts[x]<10000 and counts[x] > 100]
not_too_small = [x for x in counts.index if counts[x] > 100]
vdj1 = vdj0.loc[vdj.Epitope.isin(not_too_small)]
print('finding cdr1s...')
vdj1['cdr1a_IMGTgaps'] = [x if pd.isna(x) else sf.get_cdr1_seq_with_gaps(x.strip()) for x in vdj1.TCRa]
vdj1['cdr1b_IMGTgaps'] = [x if pd.isna(x) else sf.get_cdr1_seq_with_gaps(x.strip()) for x in vdj1.TCRb]
print('finding cdr2s...')
vdj1['cdr2a_IMGTgaps'] = [x if pd.isna(x) else sf.get_cdr2_seq_with_gaps(x.strip()) for x in vdj1.TCRa]
vdj1['cdr2b_IMGTgaps'] = [x if pd.isna(x) else sf.get_cdr2_seq_with_gaps(x.strip()) for x in vdj1.TCRb]
print('finding cdr3s...')
vdj1['cdr3a_IMGTgaps'] = [x if pd.isna(x) else sf.get_cdr3_seq_with_gaps(x.strip()) for x in vdj1.TCRa]
vdj1['cdr3b_IMGTgaps'] = [x if pd.isna(x) else sf.get_cdr3_seq_with_gaps(x.strip()) for x in vdj1.TCRb]

print('run checks...')
assert (vdj1.dropna(subset=['cdr1a_IMGTgaps']).apply(lambda x: x.cdr3a_IMGTgaps.replace('-', '') in x.TCRa, axis=1)).all()
assert (vdj1.dropna(subset=['cdr1b_IMGTgaps']).apply(lambda x: x.cdr3b_IMGTgaps.replace('-', '') in x.TCRb, axis=1)).all()
assert (vdj1.dropna(subset=['cdr2a_IMGTgaps']).apply(lambda x: x.cdr3a_IMGTgaps.replace('-', '') in x.TCRa, axis=1)).all()
assert (vdj1.dropna(subset=['cdr2b_IMGTgaps']).apply(lambda x: x.cdr3b_IMGTgaps.replace('-', '') in x.TCRb, axis=1)).all()
assert (vdj1.dropna(subset=['cdr3a_IMGTgaps']).apply(lambda x: x.cdr3a_IMGTgaps.replace('-', '') in x.TCRa, axis=1)).all()
assert (vdj1.dropna(subset=['cdr3b_IMGTgaps']).apply(lambda x: x.cdr3b_IMGTgaps.replace('-', '') in x.TCRb, axis=1)).all()

vdj2 = vdj1.dropna(subset=['cdr1a_IMGTgaps', 'cdr1b_IMGTgaps', 'cdr2a_IMGTgaps', 'cdr2b_IMGTgaps', 'cdr3a_IMGTgaps', 'cdr3b_IMGTgaps'])
# remove cdr3s that do not start with C
vdj3 = vdj2.loc[(vdj2['CDR3-a'].str[0] == 'C') & (vdj2['CDR3-b'].str[0] == 'C')]
vdj3['len_cdr3a'] = [len(x) for x in vdj3['cdr3a_IMGTgaps']]
vdj3['len_cdr3b'] = [len(x) for x in vdj3['cdr3b_IMGTgaps']]
# remove sequences that are stupidly long
print(vdj3.len_cdr3a.value_counts())# only 18 with length > 19
print(vdj3.len_cdr3b.value_counts()) # only 39 with length > 19

vdj4 = vdj3.loc[(vdj3.len_cdr3a < 20) & (vdj3.len_cdr3b < 20)]
print('saving without TCR...')
vdj4.to_csv('data/vdj_cleaned_subset_for_MI.csv')

print('calculate renumbered TCR...')
vdj4['TCRa_IMGTgaps'] = [sf.get_vregion_seq_with_gaps(x) for x in vdj4.TCRa]
vdj4['TCRb_IMGTgaps'] = [sf.get_vregion_seq_with_gaps(x) for x in vdj4.TCRb]

print('saving with TCR...')
vdj4.to_csv('data/vdj_cleaned_subset_for_MI_IMGT_TCRs.csv')