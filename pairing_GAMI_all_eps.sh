# run best with all epitopes

dist='lev'
k=20
folder='data/output/pairwise_distances_no-small-study/cdr3'
output_folder='data/output/GAMI/'
ga_folder='data/output/GA/'

# all the small epitopes here 
declare -a arrE=('GLCTLVAML' 'CINGVCWTV' 'YLQPRTFLL' 'NLVPMVATV' 'RLRAEAQVK' 'HGIRNASFI'
                    'LLWNGPMAV' 'SSPPMFRV' 'LSLRNPILV' 'SSYRRPVGI' 'LTDEMIAQY' 'TTDPSFLGRY' 
                    'SPRWYFYYL''KSKRTPMGF' 'ATDALMTGF' 'SSLENFRAYV' 'ELAGIGILTV' 'IVTDFSVIK'
                    'ASNENMETM') 

for ep in "${arrE[@]}"
do
   GA="${ga_folder}GA-epitope-${ep}_ind-no_vgene-no_distance_type-lev_kNN-20_n_repeats-100.csv.gz"
   python PairingVDJdb_GAandMI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI_no-small-study.csv --MI_repair yes --GA_thresh 0.95 --Vgene no --output "${output_folder}" --GA_precomputed $GA --MI_repeats 10 --MI_L 0.6
   python PairingVDJdb_GAandMI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI_no-small-study.csv --MI_repair yes --GA_thresh 0.95 --Vgene no --output "${output_folder}" --GA_precomputed $GA --MI_repeats 10 --MI_L 1
done

# run with large epitopes

declare -a arrE=('RAKFKQLL' 'GILGFVFTL' 'AVFDRKSDAK')
declare -a arri=(0 1 2 3 4)

for i in "${arri[@]}"
  do
   for ep in "${arrE[@]}"
   do
      GA="${GA_folder}${i}/GA-epitope-${ep}_ind-no_vgene-no_distance_type-lev_kNN-20_n_repeats-100.csv.gz"
      python PairingVDJdb_GAandMI.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --MI_repair yes --GA_thresh 0.95 --Vgene no --output "${output_folder}${i}" --GA_precomputed "$GA" --MI_repeats 10 --MI_L 0.6
      python PairingVDJdb_GAandMI.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --MI_repair yes --GA_thresh 0.95 --Vgene no --output "${output_folder}${i}" --GA_precomputed "$GA" --MI_repeats 10 --MI_L 1
   done
 done