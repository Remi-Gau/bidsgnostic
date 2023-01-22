#! /bin/bash

set -eux -o pipefail

################################################################################
# bidsgnostic group

datasets=$(find bids-examples -mindepth 1 -maxdepth 1 -type d)

DATASET_TO_SKIP=(".github .git ds116 ds000248")

for dataset in $datasets; do

    ds_name=$(basename "$dataset")

    if [[ " ${DATASET_TO_SKIP[*]} " =~ " ${ds_name} " ]]; then
        continue    
    fi    

    : "Running on $dataset"

    mkdir -p ./derivatives/${ds_name}

    bidsgnostic_layout  "$PWD/$dataset" "$PWD/derivatives/${ds_name}" group --plot_by suffix

done
