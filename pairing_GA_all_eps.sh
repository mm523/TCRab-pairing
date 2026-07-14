# run best with all epitopes

folder='data/output/pairwise_distances_no-small-study/cdr3'
dist='lev'
k=20
output_folder='data/output/GA/'

# all the small epitopes here
declare -a arrE=('GLCTLVAML' 'CINGVCWTV' 'YLQPRTFLL' 'NLVPMVATV' 'RLRAEAQVK' 'HGIRNASFI'
                    'LLWNGPMAV' 'SSPPMFRV' 'LSLRNPILV' 'SSYRRPVGI' 'LTDEMIAQY' 'TTDPSFLGRY' 
                    'SPRWYFYYL''KSKRTPMGF' 'ATDALMTGF' 'SSLENFRAYV' 'ELAGIGILTV' 'IVTDFSVIK'
                    'ASNENMETM') 

for ep in "${arrE[@]}"
  do
    python PairingVDJdb_GA.py --epitope $ep --input data/vdj_cleaned_subset_for_MI_no-small-study.csv --kNN $k --distance_type $dist --distance_dfa "${folder}/dijA_${dist}_${ep}.csv" --distance_dfb "${folder}/dijB_${dist}_${ep}.csv" --output "${output_folder}all_eps_oct2025"
  done

# run best with large epitopes

declare -a arrE=('RAKFKQLL' 'GILGFVFTL' 'AVFDRKSDAK')
declare -a arri=(0 1 2 3 4)


for i in "${arri[@]}"
  do
    for ep in "${arrE[@]}"
      do
         echo ${ep}
         python PairingVDJdb_GA.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --kNN $k --distance_type $dist --distance_dfa "${folder}/dijA_${dist}_${ep}_ss${i}_700.csv" --distance_dfb "${folder}/dijB_${dist}_${ep}_ss${i}_700.csv" --output "${output_folder}${i}"
      done
  done