: '
Script to run MI-IPA on all epitopes
'

l=0.6
l1=1
w=0.6
step=6
conf='none'
V='no'

declare -a arrE=('GLCTLVAML' 'CINGVCWTV' 'YLQPRTFLL' 'NLVPMVATV' 'RLRAEAQVK' 'HGIRNASFI' 
                  'LLWNGPMAV' 'SSPPMFRV' 'LSLRNPILV' 'SSYRRPVGI' 'LTDEMIAQY' 'TTDPSFLGRY' 
                  'SPRWYFYYL' 'KSKRTPMGF' 'ATDALMTGF' 'SSLENFRAYV' 'ELAGIGILTV' 'IVTDFSVIK' 
                  'ASNENMETM') 

for ep in "${arrE[@]}"
  do
    python PairingVDJdb_MI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --L $l --weight $w --step $step --confidence $conf --Vgene $V --output all_eps_noconf
    python PairingVDJdb_MI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --L 1 --weight $w --step $step --confidence $conf --Vgene $V --output all_eps_noconf
  done

declare -a arrE=('AVFDRKSDAK' 'GILGFVFTL' 'RAKFKQLL')
declare -a arri=(0 1 2 3 4 5 6 7 8 9)

for i in "${arri[@]}"
  do
    for ep in "${arrE[@]}"
      do
        python PairingVDJdb_MI.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --L $l --weight $w --step $step --confidence $conf --Vgene $V --output big_eps_noconf_700/${i}
        python PairingVDJdb_MI.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --L 1 --weight $w --step $step --confidence $conf --Vgene $V --output big_eps_noconf_700/${i}
      done
  done

# run best theoretical performance

### TO ADD: best perf for all eps

l=0.6
l1=1
w=0.6
step=6
conf='none'
V='no'

declare -a arrE=('GLCTLVAML' 'CINGVCWTV' 'YLQPRTFLL' 'NLVPMVATV' 'RLRAEAQVK' 'HGIRNASFI' 
                  'LLWNGPMAV' 'SSPPMFRV' 'LSLRNPILV' 'SSYRRPVGI' 'LTDEMIAQY' 'TTDPSFLGRY' 
                  'SPRWYFYYL' 'KSKRTPMGF' 'ATDALMTGF' 'SSLENFRAYV' 'ELAGIGILTV' 'IVTDFSVIK' 
                  'ASNENMETM') 

for ep in "${arrE[@]}"
  do
    python PairingVDJdb_MI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --L $l --weight $w --step $step --confidence $conf --Vgene $V --output all_eps_best/${i} --prop_test 0 --n_repeat 1
  done


declare -a arrE=('AVFDRKSDAK' 'RAKFKQLL' 'GILGFVFTL') 
declare -a arri=(0 1 2 3 4 5 6 7 8 9)

for ep in "${arrE[@]}"
  do
    for i in "${arri[@]}"
      do
        python PairingVDJdb_MI.py --epitope $ep --input data/big_epitopes_subsamples_700/${ep}_ss${i}.csv --L $l --weight $w --step $step --confidence $conf --Vgene $V --output big_eps_best_noconf_700/${i} --prop_test 0 --n_repeat 1
      done
  done