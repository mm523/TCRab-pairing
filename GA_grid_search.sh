: '
I use two epitopes because I want to check that it is general behaviour. GLC should have lots of signal, YLQ should have a 
bit less. They have 345 and 333 pairs each (so they should be equivalent in sample-size effect)

Default settings:

- vgene to no 
- small_ind to no 
- kNN to 10
- n_repeats to 100
- distance_type to lev

Parameters I should sweep:
- vgene
- kNN
- distance_type
'

# scan over distance and kNN

ep='GLCTLVAML'
ep1='YLQPRTFLL'
folder='data/output/pairwise_distances_no-small-study/cdr3'
input_data='data/vdj_cleaned_subset_for_MI_no-small-study.csv'
output_folder='data/output/GA/grid_search/'
declare -a arrD=('lev' 'tcrdist' 'triplet' 'weightedlev')
declare -a arrK=(1 5 10 20 50 100) 

for k in "${arrK[@]}"
do
   for dist in "${arrD[@]}"
   do
      echo ${conf}
      python PairingVDJdb_GA.py --epitope $ep1 --input $input_data --kNN $k --distance_type $dist --distance_dfa "${folder}/dijA_${dist}_${ep}.csv" --distance_dfb "${folder}/dijB_${dist}_${ep}.csv" --output $output_folder
      python PairingVDJdb_GA.py --epitope $ep --input $input_data --kNN $k --distance_type $dist --distance_dfa "${folder}/dijA_${dist}_${ep1}.csv" --distance_dfb "${folder}/dijB_${dist}_${ep1}.csv" --output $output_folder
   done
 done