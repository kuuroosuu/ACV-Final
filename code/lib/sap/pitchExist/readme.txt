goCollectData.m: Collect training/test data
goSelectInput.m: Select 2 inputs
goPlotDataAll.m: Data plots of all kinds
goPlotSelectedInput.m: 2D scatter data plot of the selected inputs
goQcDesign.m: QC design
goGmmcDesign.m: GMMC design


For HMM:
1. Run getTransProb.m to Compute transition probability and store it as transProb.mat
2. Run goGmmcDesign.m to get GMM for each state (Input index should be set to [1 2 3] to select all inputs)

Function list:
suvParamSet.m: Set up parameters for SUV classification
wave2feature.m: Wave to feature conversion
hmm4suv.m: A function for using HMM for SUV classification
wave2suv.m: Decode a wave file into SU and V