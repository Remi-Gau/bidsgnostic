#! /bin/bash

set -eux -o pipefail

################################################################################
# bidsngostic group

datasets=$(find bids-examples -mindepth 1 -maxdepth 1 -type d)

for dataset in $datasets; do

   if [ -e "bids-examples/.git" ]; then
       continue
   fi
   if [ -e "bids-examples/.github" ]; then
       continue
   fi   

    ds_name=$(basename "$dataset")

    : "Running on $dataset"

    mkdir -p ./derivatives/${ds_name}

    bidsgnostic_layout  "$PWD/$dataset" "$PWD/derivatives/${ds_name}" group --plot_by suffix

done

################################################################################
# bidsngostic participant

datasets=$(find bids-examples -mindepth 1 -maxdepth 1 -type d -name "ds*")

DATASET_TO_SKIP=(".github .git ds000001-fmriprep ds000246 ds000247 ds000117 ds006 ds007 ds008")

for dataset in $datasets; do

    ds_name=$(basename "$dataset")

    if [[ " ${DATASET_TO_SKIP[*]} " =~ " ${ds_name} " ]]; then
        continue    
    fi   

    : "Running on $dataset"

    mkdir -p ./derivatives/${ds_name}

    bidsgnostic  "$PWD/$dataset" "$PWD/derivatives/${ds_name}" participant --cores all --force-output

done

datasets=$(find bids-examples -mindepth 1 -maxdepth 1 -type d -name "ieeg*")

DATASET_TO_SKIP=(".github .git")

for dataset in $datasets; do

    ds_name=$(basename "$dataset")

    if [[ " ${DATASET_TO_SKIP[*]} " =~ " ${ds_name} " ]]; then
        continue    
    fi   

    : "Running on $dataset"

    mkdir -p ./derivatives/${ds_name}

    bidsgnostic  "$PWD/$dataset" "$PWD/derivatives/${ds_name}" participant --cores all --force-output

done


datasets=$(find bids-examples -mindepth 1 -maxdepth 1 -type d -name "eeg*")

DATASET_TO_SKIP=(".github .git")

for dataset in $datasets; do

    ds_name=$(basename "$dataset")

    if [[ " ${DATASET_TO_SKIP[*]} " =~ " ${ds_name} " ]]; then
        continue    
    fi   

    : "Running on $dataset"

    mkdir -p ./derivatives/${ds_name}

    bidsgnostic  "$PWD/$dataset" "$PWD/derivatives/${ds_name}" participant --cores all --force-output

done



