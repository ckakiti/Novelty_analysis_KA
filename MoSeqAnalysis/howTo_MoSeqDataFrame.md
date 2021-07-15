# [IN PROGRESS] How to work with data from MoSeqDataFrame file

## You'll need these files:
- MiceIndex.mat (generated with [extract_uuid.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/extract_uuid.m))
- MoSeqDataFrame.mat (generated with [ModelDataTransfer.py](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/ModelDataTransfer.py))
- depth_ts.txt (for each session, primary file generated when collecting data)
- rgb_ts.txt (for each session, primary file generated when collecting data)

## About MiceIndex.mat:
- MATLAB structure containing information about each mouse and each session that went into MoSeq analysis batch
  - animal name
  - novelty condition (e.g. S, C for stimulus, contextual)
  - ExpDay (sub-structure with rows for session date and MSid)

## About MoSeqDataFrame.mat
- MATLAB structure containing output from MoSeq analysis
- relevant fields:
