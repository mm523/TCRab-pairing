: '
Script to run GA on all epitopes
'

# All epitope list:
# scan over distance and kNN

folder='data/output/pairwise_distances/cdr3'
dist='lev'
k=20
declare -a arrE=('GLCTLVAML' 'LLWNGPMAV' 'HGIRNASFI'
                 'CINGVCWTV' 'YLQPRTFLL' 'NLVPMVATV' 'RLRAEAQVK'  
                 'SSPPMFRV' 'LSLRNPILV' 'SSYRRPVGI' 'LTDEMIAQY' 
                 'TTDPSFLGRY' 'SPRWYFYYL' 'KSKRTPMGF' 'ATDALMTGF' 
                 'SSLENFRAYV' 'ELAGIGILTV' 'IVTDFSVIK' 'ASNENMETM')

for ep in "${arrE[@]}"
do
   echo ${ep}
   python PairingVDJdb_GA.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --kNN $k --distance_type $dist --distance_dfa "${folder}/dijA_${dist}_${ep}.csv" --distance_dfb "${folder}/dijB_${dist}_${ep}.csv" --output all_eps
done

declare -a arrE=('AVFDRKSDAK' 'RAKFKQLL' 'GILGFVFTL' 'KLGGALQAK') #
declare -a arri=(0 1 2 3 4 5 6 7 8 9)

for i in "${arri[@]}"
  do
    for ep in "${arrE[@]}"
      do
         echo ${ep}
         python PairingVDJdb_GA.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --kNN $k --distance_type $dist --distance_dfa "${folder}/dijA_${dist}_${ep}_ss${i}_700.csv" --distance_dfb "${folder}/dijB_${dist}_${ep}_ss${i}_700.csv" --output big_eps_700/${i}
      done
  done