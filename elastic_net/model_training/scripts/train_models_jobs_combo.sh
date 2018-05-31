#!/bin/bash

for tiss in `cat mesa_pops2.txt`
do
    echo $tiss
for peer in `cat mesa_peers.txt`
do
    echo $peer
for covar in `cat mesa_covars.txt`
do
    echo $covar
    for chrom in {18..22}
    do
        echo -e "    $chrom"
        qsub pop_train_models_combo.pbs  -N mesa_training_${tiss}_${chrom} -v tiss=${tiss},chrom=${chrom},peer=${peer},covar=${covar}
        sleep .5
    done
done
done
done
