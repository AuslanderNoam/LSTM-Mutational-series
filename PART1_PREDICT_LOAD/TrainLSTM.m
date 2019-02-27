function [tnet,XTS,YTS,sm] = TrainLSTM(dim1,MUTD2)
%%%%% Training LSTM to predict mutational load

%%%%% INPUT:
%%% dim1 - dimension, number of mutations to use for training
%%% MUTD2 - data structure for training, the order discrete mutational data

%%%%% OUTPUT:
%%% tnet - the trained network
%%% XTS - the training data, disctere mutational data) 
%%% YTS - the training labels, high(1) vs. low(0) mutational load
%%% sm - the training mutational load, non discrete


MUTDR = MUTD2; 
MUTD2.tab = MUTD2.tab(1:dim1,:); %%Cut the data for training, use dim1 mutations
[~,mm] = size(MUTD2.tab);


sm = sum(MUTDR.tab>0)'; %%Get load


T2={};L2={};
%%%prepare training data format for LSTM
for j = 1:mm
    T2{j} = MUTD2.tab(1:dim1,j)';
end

%%%prepare label data format for LSTM
L1 = zeros(size(sm));
L1(sm>quantile(sm,0.5))=1;
L2 = categorical(L1);

maxEpochs = 100; %%%Can change 

miniBatchSize = 27;
inputSize = 1;
numHiddenUnits = 5;
numClasses = length(unique(L1));
layers = [ ...
    sequenceInputLayer(inputSize)
    lstmLayer(numHiddenUnits,'OutputMode','last')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];


options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Verbose',0);



XTS = T2;
YTS = L2;
tnet = trainNetwork(XTS,YTS,layers,options);    %%%training the network
