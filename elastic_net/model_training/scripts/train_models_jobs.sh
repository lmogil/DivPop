#!/bin/bash

for tiss in `cat mesa_pops2.txt`
do
    echo $tiss
    for chrom in {22..22}
    do
        echo -e "    $chrom"
        qsub pop_train_models.pbs  -N mesa_training_${tiss}_${chrom} -v tiss=${tiss},chrom=${chrom}
        sleep .5
    done
done

