#!/bin/bash

# chmod +x Shell_Script_Template.sh
# cd ~/Programs/Novelty_analysis_KA/MoSeqAnalysis/
# ./Shell_Script_Template.sh

# STEP 1: EXTRACTION
# to quickly get path to all .dat files nested under folder, enter this into terminal:
# find -type d -printf '%d\t%P\n' | sort -r -nk1 | cut -f2-
# to get all .mp4 files:
# find -type f -regex '.*.mp4' -printf '%d\t%P\n' | sort -r -nk1 | cut -f2-
# to delete all .mp4 files under current folder:
# find -type f -name '*.mp4' -delete

# no fiber:
moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Dataset_name/Dataset_name_MoSeq/ path_to_file /depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl --bg-roi-dilate 75 75 --frame-trim 900 0

# fiber:
moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Dataset_name/Dataset_name_MoSeq/ path_to_file /depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl --bg-roi-dilate 75 75 --frame-trim 900 0 --use-tracking-model True --cable-filter-iters 1


# STEP 2: PCA
# for small datasets
moseq2-pca train-pca -c 6 -n 1 &> train-pca_output
moseq2-pca apply-pca -c 6 -n 1
moseq2-pca compute-changepoints -c 6 -n 1

# for larger datasets:
moseq2-pca train-pca -c 6 -n 1 -m 25GB &> train-pca_output

# fiber:
moseq2-pca train-pca -c 6 -n 1 --missing-data &> train-pca_output


# STEP 3: MODEL
# set kappa to number of frames listed during train-pca step (saved in train-pca_output)
moseq2-model learn-model --kappa 9999999 --save-model _pca/pca_scores.h5 my_model.p

# for larger datasets:
moseq2-model learn-model --kappa 9999999 --save-model _pca/pca_scores.h5 my_model.p --whiten each &> learn-model_output


# STEP 4: VISUALIZATION
moseq2-viz generate-index 
moseq2-viz make-crowd-movies --max-syllable 1000 --sort False moseq2-index.yaml my_model.p

# to extract uuids from moseq_index file:
# open moseq2-index.yaml in Visual Studio
# CTRL F, enter this into box (path:\r?$)|(.*00.h5)|(.*yaml)|(\s uuid.*)
# CTRL + SHIFT + L to select all occurrances found
# CTRL + C --> CTRL + V into new file
# run extract_uuid matlab script


# EXTRAS
#moseq2-viz add-group -k SubjectName -v "Alcohol" -v "Amine" -v "Ketone" -v "Fred" -v "Harry" -v "Neville" -g "group1" moseq2-index.yaml
#moseq2-viz add-group -k SubjectName -v "Aldehyde" -v "Ester" -v "Thiol" -v "George" -v "Hermione" -v "Ron" -g "group2" moseq2-index.yaml
#moseq2-viz plot-usages moseq2-index.yaml my_model.p --group group1 --group group2
#moseq2-viz plot-usages moseq2-index.yaml my_model.p
#moseq2-viz plot-scalar-summary moseq2-index.yaml
#moseq2-viz plot-transition-graph moseq2-index.yaml my_model.p --group group1 --group group2
