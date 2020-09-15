#/bin/bash

fdata=~/projects/vRNAsite/scripts/data_fasta
vdata=~/projects/vRNAsite/scripts/data_vRNAsite
sdata=~/projects/vRNAsite/scripts/data_SPLASH

cmpA=~/projects/vRNAsite/scripts/compare_A
cmpUR=~/projects/vRNAsite/scripts/compare_UR

vRNAsite=~/projects/vRNAsite/scripts/vRNAsite

# do vRNAsite
declare -a experiments=("PR8-Udorn_PB1+NA" "PR8-Wyo03_NA+PB1" "PR8-Wyo03_NAudsub+Ud_PB1" "PR8-Wyo03_NAudsub+PB1" "Wyo03")
for exp in "${experiments[@]}"
do
    echo "DO folding"
    out=$vRNAsite/$exp
    inf=$fdata/${exp}_FASTA.fa
    pcl=$vRNAsite/$exp/${exp}_heats.pcl
    tsv=$vRNAsite/$exp/${exp}.tsv
    nice -n 10 vRNAsite.py -pfx $out -fsa $inf -thr 48 -ovr -wgm -vrh -vrs -vre -18.0 -plh -pls 2.0 -bkh -bks 0.5 -mat -cdl
    cp $tsv ${vdata}/${exp}_vRNAsite.tsv
done

# do compare
declare -a experiments=("PR8-Udorn_PB1+NA" "PR8-Wyo03_NA+PB1" "PR8-Wyo03_NAudsub+Ud_PB1" "PR8-Wyo03_NAudsub+PB1" "Wyo03")
for exp in "${experiments[@]}"
do
    echo "DO comparison_A"
    out=$cmpA/$exp
    vrs=$vdata/${exp}_vRNAsite.tsv
    sps=$sdata/${exp}_SPLASH.tsv
    nice -n 10 compare_x2.py -pfx $out -vrs $vrs -sps $sps -ort A -ovr -mle 5
    echo "DO comparison_UR"
    out=$cmpUR/$exp
    nice -n 10 compare_x2.py -pfx $out -vrs $vrs -sps $sps -ort UR -ovr -mle 5
done

# do mfe compare
mdir=~/projects/vRNAsite/scripts/SPLASH_mfe_rank_comparison
declare -a experiments=("WSN-Rep1" "WSN-Rep2" "PR8-Rep1" "PR8-Rep2" "Udorn-Rep1" "Udorn-Rep2" "PR8-Udorn_PB1+NA" "PR8-Wyo03_NA+PB1" "PR8-Wyo03_NAudsub+Ud_PB1" "PR8-Wyo03_NAudsub+PB1" "Wyo03")
for exp in "${experiments[@]}"
do
    echo "DO comparison_A"
    out=$mdir/$exp
    vrs=$vdata/${exp}_vRNAsite.tsv
    sps=$sdata/${exp}_SPLASH.tsv
    nice -n 10 compare_x2.py -pfx $out -vrs $vrs -sps $sps -ort A -ovr -mle 5 -cdr 1 -msp
    echo "DO comparison_UR"
    out=$mdir/$exp
    nice -n 10 compare_x2.py -pfx $out -vrs $vrs -sps $sps -ort UR -ovr -mle 5 -cdr 1 -msp
done
