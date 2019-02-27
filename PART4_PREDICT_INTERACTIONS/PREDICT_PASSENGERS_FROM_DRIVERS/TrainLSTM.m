function [tnet,XTS,YTS] = TrainLSTM(dim1,MUTD2,L1)
%%%%% TrainLSTM - Training LSTM to predict passesnger mutation fron a sequence of drivers up to dim1

%%%%% INPUT:
%%% dim1 - dimension, current mutation to predict in the sequence
%%% MUTD2 - data structure for training, the order discrete mutational data
%%% L1 - labels, occurrence of the dim1 mutation in training set

%%%%% OUTPUT:
%%% tnet - the trained network
%%% XTS - the training data, disctere mutational data) 
%%% YTS - the training labels, occurrence of the dim1 mutation in training set

MUTDR = MUTD2;
MUTD2.tab = MUTD2.tab(1:dim1,:);%%%Use driver mutations up to dim1
[nn,mm] = size(MUTD2.tab);


T2={};L2={};
%%%prepare training data format for LSTM
for j = 1:mm
    T2{j} = MUTD2.tab(1:dim1,j)';

end

%%%prepare label data format for LSTM
L2 = categorical(L1);


maxEpochs = 10;%%%LOW. Increase to improve performance
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
class(YTS)
tnet = trainNetwork(XTS,YTS,layers,options);%%%training the network
