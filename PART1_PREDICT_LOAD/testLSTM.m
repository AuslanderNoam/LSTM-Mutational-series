function [AUC,scores,sm] = testLSTM(tnet,MUTD2,dim1)
%%%%% Get test results fron trained LSTM for mutational load prediction

%%%%% INPUT:
%%% tnet - the trained network
%%% MUTD2 - data structure for testing, the order discrete mutational data
%%% dim1 - dimension, number of mutations to use for testing

%%%%% OUTPUT:
%%% AUC - the resulting test AUC
%%% scores - the classification scores
%%% sm - the test mutational load, non discrete
MUTDR = MUTD2;

MUTD2.tab = MUTD2.tab(1:dim1,:); %%Cut the test data, use dim1 mutations

[nn,mm] = size(MUTD2.tab>0);

sm = sum(MUTDR.tab)'; %%Get load

T2={};L2={};
%%%prepare training data format for LSTM
for j = 1:mm
    T2{j} = MUTD2.tab(1:dim1,j)';
end

%%%prepare label data format for LSTM
L1 = zeros(size(sm));
L1(sm>quantile(sm,0.5))=1;
L2 = categorical(sm);



ll=[];scc=[];
[YPred,scores] = classify(tnet,T2); %%%perform testing
[n,m] = size(scores);

for i = 1:n
    scc(i) = scores(i,1);
    ll(i)=L2(i);
end


[X,Y,T,AUC] = perfcurve(ll',scc',1);