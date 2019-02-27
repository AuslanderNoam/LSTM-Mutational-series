function [AUC,scores] = getAUCLSTM2(tnet,MUTD2,dim1,L1)
%%%%% Get test AUC fron trained LSTM next mutation prediction
%%%%% (using 'Trail' last mutation each time)

%%%%% INPUT:
%%% tnet - the trained network
%%% MUTD2 - data structure for training, the order discrete mutational data
%%% dim1 - dimension, current mutation to predict in the sequence
%%% L1 - labels, occurrence of the dim1 mutation in training set

%%%%% OUTPUT:
%%% AUC - the resulting test AUC
%%% scores - the classification scores


MUTDR = MUTD2;

Trail = 100;
MUTD2.tab = MUTD2.tab(1:dim1-1,:); %%%Use mutations up to dim1-1

[nn,mm] = size(MUTD2.tab>0);

T2={};L2={};
%%%prepare training data format for LSTM
for j = 1:mm
    T2{j} = MUTD2.tab(dim1-Trail:dim1-1,j)';  %%%get Trail last mutations up to dim1-1
end


%%%prepare label data format for LSTM
L2 = categorical(L1);



c1=categorical(1);
c0=categorical(0);
ll=[];scc=[];
[YPred,scores] = classify(tnet,T2);
[n,m] = size(scores);
for i = 1:n
    scc(i) = scores(i,2);
    ll(i)=L2(i);
end

[X,Y,T,AUC] = perfcurve(ll',scc',1);