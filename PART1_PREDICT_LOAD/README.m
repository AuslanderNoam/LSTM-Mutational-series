% % % PART1_PREDICT_LOAD - predict mutational load from a time series of
% mutations, and plotting relevant figures. 




% % % % % A. Plot Figures

% 1. FIGURE1AD - plots figure1A-D
% 2. FIGURE1EG - plots figure1E-G
% 3. FIGURE1HL - plots figure1H-L (clinical analysis of test set)





% % % % % B. Scripts to run:

% 1. Run_PRED_LOAD_LSTM
%Training LSTMs to predict the mutational load, colon&lung 

% 2. Run_PRED_LOAD_LINEAR_CLASSIFIERS
%Training linear classifiers to predict the mutational load, colon&lung 

% 3. Run_PRED_LOAD_SHUFFLE_DATA
%Training classifiers to predict the mutational load based on randonly selected K mutations, colon&lung 





% % % % % C. Functions:

% 1. TrainLSTM
%%%%% Training LSTM to predict mutational load

% 2. testLSTM
%%%%% Get test results fron trained LSTM for mutational load prediction





% % % % % D. Data structures
% COAD2.mat - TCGA COAD training data
% MGI2.mat - COAD test1
% DFCI2.mat - COAD test2
% LUAD2.mat - TCGA LUAD training data
% BROAD2.mat - LUAD test

% MGI_CLIN - COAD test1 clinical data
% DFCI_CLIN - COAD test2 clinical data
% BROAD_CLIN - LUAD test clinical data


%%%**Other data structures are saved result structures


