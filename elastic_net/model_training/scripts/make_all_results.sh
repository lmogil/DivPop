#!/bin/bash
# Script to concatenate result files split by chromosome together.

tissue=$1
allResults=$2

i=0
for resultsfile in $(ls /home/lauren/PredictDB_Pipeline_GTEx_v7/new_output/${tissue}_nested_cv_chr*_model_summaries*); do
        if [ $i -eq 0 ] ; then
                head -n 1 $resultsfile > $allResults
                i=1
        fi
        tail -n +2 $resultsfile >> $allResults
done

