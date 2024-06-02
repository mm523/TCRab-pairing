: '
Script to run GA+MI-IPA on all epitopes
'

# All epitope list:

folder='data/output/pairwise_distances/cdr3'
dist='lev'
k=20
declare -a arrE=('GLCTLVAML' 'LLWNGPMAV' 'HGIRNASFI'
                 'CINGVCWTV' 'YLQPRTFLL' 'NLVPMVATV' 'RLRAEAQVK'  
                 'SSPPMFRV' 'LSLRNPILV' 'SSYRRPVGI' 'LTDEMIAQY' 
                 'TTDPSFLGRY' 'SPRWYFYYL' 'KSKRTPMGF' 'ATDALMTGF' 
                 'SSLENFRAYV' 'ELAGIGILTV' 'IVTDFSVIK' 'ASNENMETM'
                 )

for ep in "${arrE[@]}"
do
   GA="data/output/pairing_GA/all_eps/GA-epitope-${ep}_ind-no_vgene-no_distance_type-lev_kNN-20_n_repeats-100.csv.gz"
   # with learning in MI-IPA
   python PairingVDJdb_GAandMI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --MI_repair yes --GA_thresh 0.95 --Vgene no --output all_eps_L06 --GA_precomputed $GA --MI_repeats 10 --MI_L 0.6
   # with no learning in MI-IPA
   python PairingVDJdb_GAandMI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --MI_repair yes --GA_thresh 0.95 --Vgene no --output all_eps_L1 --GA_precomputed $GA --MI_repeats 10 --MI_L 1
done

declare -a arrE=('GILGFVFTL'  'AVFDRKSDAK' 'RAKFKQLL')
declare -a arri=(0 1 2 3 4 5 6 7 8 9)

for i in "${arri[@]}"
  do
   for ep in "${arrE[@]}"
   do
      GA="/home/marti/OneDrive-personal/Desktop/OneDrive - University College London/Desktop/PhD/AB-interaction-project/pairing_output_folder/pairing_GA/big_eps_700/${i}/GA-epitope-${ep}_ind-no_vgene-no_distance_type-lev_kNN-20_n_repeats-100.csv.gz"
      # with learning in MI-IPA
      python PairingVDJdb_GAandMI.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --MI_repair yes --GA_thresh 0.95 --Vgene no --output big_eps_700_L06/${i} --GA_precomputed "$GA" --MI_repeats 10 --MI_L 0.6
      # with no learning in MI-IPA
      python PairingVDJdb_GAandMI.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --MI_repair yes --GA_thresh 0.95 --Vgene no --output big_eps_700_L1/${i} --GA_precomputed "$GA" --MI_repeats 10 --MI_L 1
   done
 done