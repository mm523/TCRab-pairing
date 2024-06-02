#%%
import pandas as pd
# import matplotlib.pyplot as plt
# from collections import Counter
# %%
vdj = pd.read_csv('data/vdj-export-03022023.tsv', sep='\t')
#%%
print('separate alpha/beta and create complexes')
vdja = vdj.loc[vdj.Gene=='TRA']
vdjb = vdj.loc[vdj.Gene=='TRB']

cols = ['Gene', 'CDR3', 'V', 'J', 'CDR3fix']

vdja.columns = [x + '-a' if x in cols else x for x in vdja.columns]
vdjb.columns = [x + '-b' if x in cols else x for x in vdjb.columns]
#%%
ab_asone = pd.merge(vdja.loc[vdja['complex.id'] != 0], vdjb.loc[vdjb['complex.id'] != 0]).drop_duplicates(subset = ['CDR3-a', 'V-a', 'J-a', 'CDR3-b', 'V-b', 'J-b', 'Epitope'])
#%%
# plt.bar(Counter(ab_asone.Epitope).keys(), Counter(ab_asone.Epitope).values())
#%%
print('Keep only epitopes for which I have >100 sequences')
big_eps = [x for x in set(ab_asone.Epitope) if ab_asone.Epitope.value_counts()[x] > 100]
ab_asone_big = ab_asone.loc[ab_asone.Epitope.isin(big_eps)] # down to 19 epitopes
# %%
charac_a = {y for x in ab_asone_big['CDR3-a'] for y in x}
charac_b = {y for x in ab_asone_big['CDR3-b'] for y in x}
AAs = 'ARNDCEQGHILKMFPSTWYV'

assert len(charac_a) <= 20
assert sorted(list(charac_a)) == sorted(list(AAs))
assert len(charac_b) <= 20
assert sorted(list(charac_b)) == sorted(list(AAs))

# %%
ab_asone_big.to_csv('data/vdj_cleaned.csv')