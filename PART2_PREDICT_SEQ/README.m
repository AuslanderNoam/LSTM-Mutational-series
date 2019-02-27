% % % 2. PART2_PREDICT_SEQ - predict the next mutation in the time-sequence
% step-by-step, and plotting relevant figures. 




% % % % % A. Plot Figures

% FIGURE2AC - plots figure2A-C





% % % % % B. Scripts to run:

% 1. PRED_ALL_SEQ_PART1
%Training LSTMs to predict the next mutation in the sequence, colon&lung
%for the first 1000 mutations (last 1000 score-ordered)

% 2. PRED_ALL_SEQ_PART2
%Training LSTMs to predict the next mutation in the sequence, colon&lung
%for the mutations 1001-end 





% % % % % C. Functions:

% 1. TrainLSTM
%%%%% raining LSTM to predict mutation dim1 fron a sequence up to dim1-1

% 2. TrainLSTM2
%%%%% raining LSTM to predict mutation dim1 fron a sequence up to dim1-1 (using 'Trail' last mutation each time)

% 3. getAUCLSTM
%%%%% Get test AUC fron trained LSTM next mutation prediction

% 4. getAUCLSTM2
%%%%% Get test AUC fron trained LSTM next mutation prediction (using 'Trail' last mutation each time)





% % % % % D. Data structures
% COAD2.mat - TCGA COAD training data
% MGI2.mat - COAD test1
% DFCI2.mat - COAD test2
% LUAD2.mat - TCGA LUAD training data
% BROAD2.mat - LUAD test

% COLONDRIVERS - colon cancer drivers
% LUADDRIVERS - lung cancer drivers


%%%**Other data structures are saved result structures


