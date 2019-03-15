#!/bin/bash

# chmod +x Shell_Script_Template.sh
# find -type d -printf '%d\t%P\n' | sort -r -nk1 | cut -f2-
# cd ~/Programs/Novelty_analysis_KA/MoSeqAnalysis/
# ./Shell_Script_Template.sh

#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/7day_preexposure_MoSeq/S6_Neville/181125/session_20181125122500/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl


moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Test/session_20180926174440/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Test/session_20180924185040/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Test/session_20180922154525/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl

#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Miami/190126/session_20190126115413/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Miami/190125/session_20190125135727/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Miami/190124/session_20190124153413/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Kennebunk/190216/session_20190216133155/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Kennebunk/190214/session_20190214171921/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Juneau/190220/session_20190220112655/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Juneau/190219/session_20190219150148/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Juneau/190218/session_20190219145922/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Juneau/190218/session_20190218121232/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Ithaca/190217/session_20190217130937/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Ithaca/190216/session_20190216144245/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Ithaca/190215/session_20190215135915/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Houston/190218/session_20190218113654/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Houston/190217/session_20190217134349/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Houston/190216/session_20190216151902/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Gardner/190217/session_20190217123133/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Gardner/190216/session_20190216140726/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Gardner/190215/session_20190215132339/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Francisco/190220/session_20190220120254/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#oseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Francisco/190219/session_20190219153503/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl
#moseq2-extract extract /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Francisco/190218/session_20190218143552/depth.dat --flip-classifier /home/alex/moseq2/flip_classifier_k2_c57_10to13weeks.pkl

#moseq2-pca train-pca -c 6 -n 1
#moseq2-pca apply-pca -c 6 -n 1
#moseq2-pca compute-changepoints -c 6 -n 1
#moseq2-model learn-model --kappa 539637 --save-model _pca/pca_scores.h5 my_model.p
#moseq2-model learn-model --kappa 4610752 --save-model _pca/pca_scores.h5 my_model.p

#moseq2-viz generate-index 
#moseq2-viz add-group -k SubjectName -v "DRD_37mo_seq" -g "group2" moseq2-index.yaml
moseq2-viz add-group -k SubjectName -v "Alcohol" -v "Amine" -v "Ketone" -v "Fred" -v "Harry" -v "Neville" -g "group1" moseq2-index.yaml
moseq2-viz add-group -k SubjectName -v "Aldehyde" -v "Ester" -v "Thiol" -v "George" -v "Hermione" -v "Ron" -g "group2" moseq2-index.yaml
#moseq2-viz plot-usages moseq2-index.yaml my_model.p --group group1 --group group2
#moseq2-viz plot-usages moseq2-index.yaml my_model.p
#moseq2-viz plot-scalar-summary moseq2-index.yaml
#moseq2-viz plot-transition-graph moseq2-index.yaml my_model.p --group group1 --group group2
#moseq2-viz make-crowd-movies --max-syllable 1000 --sort False moseq2-index.yaml my_model.p
