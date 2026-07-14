# run best with all epitopes

l=0.6
l1=1
w=0.6
step=6
conf='none'
V='no'
output_folder='/home/marti/PairingOutputNoSmallStudy/MIIPA/'

# all the small epitopes here 
declare -a arrE=('GLCTLVAML' 'CINGVCWTV' 'YLQPRTFLL' 'NLVPMVATV' 'RLRAEAQVK' 'HGIRNASFI'
                    'LLWNGPMAV' 'SSPPMFRV' 'LSLRNPILV' 'SSYRRPVGI' 'LTDEMIAQY' 'TTDPSFLGRY' 
                    'SPRWYFYYL''KSKRTPMGF' 'ATDALMTGF' 'SSLENFRAYV' 'ELAGIGILTV' 'IVTDFSVIK'
                    'ASNENMETM') 

for ep in "${arrE[@]}"
  do
    python PairingVDJdb_MI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI_no-small-study.csv --L $l --weight $w --step $step --confidence $conf --Vgene $V --output "${output_folder}"
    python PairingVDJdb_MI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI_no-small-study.csv --L 1 --weight $w --step $step --confidence $conf --Vgene $V --output "${output_folder}"
    # best performance
    python PairingVDJdb_MI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI_no-small-study.csv --L $l --weight $w --step $step --confidence $conf --Vgene $V --output "${output_folder}" --prop_test 0 --n_repeat 1
  done

# run best with large epitopes

declare -a arrE=('RAKFKQLL' 'GILGFVFTL' 'AVFDRKSDAK')
declare -a arri=(0 1 2 3 4)

for i in "${arri[@]}"
  do
    for ep in "${arrE[@]}"
      do
        python PairingVDJdb_MI.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --L $l --weight $w --step $step --confidence $conf --Vgene $V --output "${output_folder}${i}"
        python PairingVDJdb_MI.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --L 1 --weight $w --step $step --confidence $conf --Vgene $V --output "${output_folder}${i}"
        # best performance
        python PairingVDJdb_MI.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --L $l --weight $w --step $step --confidence $conf --Vgene $V --output "${output_folder}${i}" --prop_test 0 --n_repeat 1
      done
  done