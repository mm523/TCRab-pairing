#%%
import pandas as pd
import os
# %%
vdj = pd.read_csv('data/vdj_cleaned.csv')
# %%
import subprocess
# %%
os.chdir('../stitchr/') # assumes stitchr is cloned in directory at the same level as this repo and code is run from outside the data_preprocessing directory
str1 = 'conda activate MutualInformation'
str2 = ' python Scripts/stitchr.py -m AA'
#%%
species_dict = {'HomoSapiens':'HUMAN', 'MusMusculus':'MOUSE'}
for i in range(vdj.shape[0]):
    if i%10 == 0:
        print(i, '/', vdj.shape[0])
    line = vdj.iloc[i]
    va, ja, cdr3a, species = line['V-a'], line['J-a'], line['CDR3-a'], line['Species']
    vb, jb, cdr3b = line['V-b'], line['J-b'], line['CDR3-b']
    if (not pd.isna(va)) and (not pd.isna(ja)) and (not pd.isna(cdr3a)):
        str3 = ' -v ' + va + ' -j ' + ja + ' -cdr3 ' + cdr3a + ' -s ' + species_dict[species]
        mycommand = str2 + str3
        x = os.popen(mycommand)
        tcra = x.read()
        vdj.at[vdj.index[i], 'TCRa'] = tcra
    else:
        print('skipping alpha')
    if (not pd.isna(vb)) and (not pd.isna(jb)) and (not pd.isna(cdr3b)):
        str3 = ' -v ' + vb + ' -j ' + jb + ' -cdr3 ' + cdr3b + ' -s ' + species_dict[species]
        mycommand = str2 + str3
        x = os.popen(mycommand)
        tcrb = x.read()
        vdj.at[vdj.index[i], 'TCRb'] = tcrb
    else:
        print('skipping beta')
# %%
print('saving...')
os.chdir('../TCRab-pairing/')
vdj.to_csv('data/vdj_cleaned_set_with_stitchr_seq.csv')