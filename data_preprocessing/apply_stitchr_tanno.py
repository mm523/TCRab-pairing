#%%
import pandas as pd
import os
import tidytcells as tt
# %%
tanno = pd.read_csv('data/Tanno_combined.csv.gz', index_col=0)
ind = 'A1 naive'
tanno = tanno.loc[tanno['sample'] == ind]
for c in ['VL', 'VH', 'JL', 'JH']:
    print(c)
    tanno[c] = tanno.apply(
        lambda row: pd.NA if type(row[c]) != str else tt.tcr.standardise(
            gene_name=row[c],
            species='HomoSapiens',
            precision='allele'
        ),
        axis=1
    )
tanno = tanno.dropna(subset=['VL', 'VH', 'JL', 'JH'])
print('shape before tidytcells: ', tanno.shape, '; shape after tidytcells: ', tanno.shape)
print(tanno.shape)
# %%
import subprocess
# %%
os.chdir('../stitchr/') # assumes stitchr is cloned in directory at the same level as this repo and code is run from outside the data_preprocessing directory
str1 = 'conda activate MutualInformation'
str2 = ' python Scripts/stitchr.py -m AA'
#%%
species_dict = {'HomoSapiens':'HUMAN', 'MusMusculus':'MOUSE'}
for i in range(tanno.shape[0]):
    if i%10 == 0:
        print(i, '/', tanno.shape[0])
    line = tanno.iloc[i]
    va, ja, cdr3a = line['VL'], line['JL'], line['CDRL3_AA']
    vb, jb, cdr3b = line['VH'], line['JH'], line['CDRH3_AA']
    if (not pd.isna(va)) and (not pd.isna(ja)) and (not pd.isna(cdr3a)):
        str3 = ' -v ' + va + ' -j ' + ja + ' -cdr3 ' + cdr3a + ' -s HUMAN' 
        mycommand = str2 + str3
        x = os.popen(mycommand)
        tcra = x.read()
        tanno.at[tanno.index[i], 'TCRa'] = tcra
    else:
        print('skipping alpha')
    if (not pd.isna(vb)) and (not pd.isna(jb)) and (not pd.isna(cdr3b)):
        str3 = ' -v ' + vb + ' -j ' + jb + ' -cdr3 ' + cdr3b + ' -s HUMAN'
        mycommand = str2 + str3
        x = os.popen(mycommand)
        tcrb = x.read()
        tanno.at[tanno.index[i], 'TCRb'] = tcrb
    else:
        print('skipping beta')
# %%
print(tanno.dropna(subset=['TCRa','TCRb'], how='any').shape)
print('saving...')
os.chdir('../TCRab-pairing/')
tanno.to_csv('data/tanno_A1naive_with_stitchr_seq.csv.gz')