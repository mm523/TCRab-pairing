: '
Script to run GA on two benchmarking epitopes to find best k and best dist
'

# scan over distance and kNN

ep='GLCTLVAML'
ep1='YLQPRTFLL'
folder='data/output/pairwise_distances/cdr3'
declare -a arrD=('triplet' 'weightedlev' 'lev' 'tcrdist')
declare -a arrK=(1 5 10 20 50 100) 

for k in "${arrK[@]}"
do
   for dist in "${arrD[@]}"
   do
      echo ${conf}
      python PairingVDJdb_GA.py --epitope $ep1 --input data/vdj_cleaned_subset_for_MI.csv --kNN $k --distance_type $dist --distance_dfa "${folder}/dijA_${dist}_${ep}.csv" --distance_dfb "${folder}/dijB_${dist}_${ep}.csv" --output dist_screen
      python PairingVDJdb_GA.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --kNN $k --distance_type $dist --distance_dfa "${folder}/dijA_${dist}_${ep1}.csv" --distance_dfb "${folder}/dijB_${dist}_${ep1}.csv" --output dist_screen
   done
 done