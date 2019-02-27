% % % 1. PREDICT_PASSENGERS_FROM_DRIVERS - using the sequence of drivers
% to predict each passanger, for lung and colon.
% !!!!!!!!!!!!!!WARNING!!!! THESE CODES REQUIRES A VERY LONG RUNNING TIME (may take A
% few days on HPC). The results for both lung and colon are saved as
% 'RES.mat' in the COLON, LUNG directories



% % % % % B. Scripts to run:

% 1. RUN_FOR_COLON - Predict passengers from drivers for colon cancer
% For each driver genes (i=1-42) train and test LSTMs predicting all
% possible passenger genes from the sequence of drivers up to i


% 2. RUN_FOR_LUNG - Predict passengers from drivers for lung cancer
% For each driver genes (i=1-26) train and test LSTMs predicting all
% possible passenger genes from the sequence of drivers up to i




% % % % % C. Functions:

% 1. RUNPREDPASS_COLON
% Training LSTM to predict each passanger from driver mutations i for
% colon

% 2. RUNPREDPASS_LUNG
% Training LSTM to predict each passanger from driver mutations i for
% lung

% 3. TrainLSTM - Training LSTM to predict passesnger mutation fron a sequence of drivers up to dim1

% 3. getAUCLSTM - Get test AUC fron trained LSTM to predict passesnger mutation fron a sequence of drivers up to dim1






% % % % % D. Data structures
% COAD2.mat - TCGA COAD passanger mutations (training labels) data
% COADd.mat - TCGA COAD ordered driver mutations (training features) data

% MGI2.mat - COAD passanger mutations (test1 labels) data
% MGId.mat - COAD ordered driver mutations (test1 features) data

% DFCI2.mat - COAD passanger mutations (test2 labels) data
% DFCId.mat - COAD ordered driver mutations (test2 features) data

% LUAD2.mat - TCGA LUAD passanger mutations (training labels) data
% LUADd.mat - TCGA LUAD ordered driver mutations (training features) data

% BROAD2.mat - LUAD passanger mutations (test labels) data
% BORADd.mat - LUAD ordered driver mutations (test features) data

