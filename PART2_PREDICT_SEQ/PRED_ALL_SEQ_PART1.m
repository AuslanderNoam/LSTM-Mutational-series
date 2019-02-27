%%%% Training LSTMs to predict the next mutation in the sequence, colon&lung
%%%% for the first 1000 mutations (last 1000 score-ordered)

%%%%%%% Predict Sequence for COAD data

%%%Load data structures
load('COAD2.mat') %%training
load('MGI2.mat')  %%test1
load('DFCI2.mat') %%test2

RES_COAD = struct
for i = 2:1000
    i
    LT = COAD2.tab(i,:)'; %%Label - the ith mutation
    [tnet] = TrainLSTM(i,COAD2,LT); %%Train LSTM to predict ith mutation 
    
    L1 = MGI2.tab(i,:)';
    [RES_COAD.AUC1(i),sc] = getAUCLSTM(tnet,MGI2,i,L1); %%%Test on test set 1
    RES_COAD.sc1(:,i)=sc(:,2);
    RES_COAD.L1(:,i) = L1;
    
    L2 = DFCI2.tab(i,:)';
    [RES_COAD.AUC2(i),sc] = getAUCLSTM(tnet,DFCI2,i,L2); %%%Test on test set 2
    RES_COAD.sc2(:,i)=sc(:,2);
    RES_COAD.L2(:,i) = L2;


end


%%%%%%% Predict Sequence for LUAD data

%%Load data structures
load('BROAD2.mat')  %%training
load('LUAD2.mat')   %%test


RES_LUAD = struct;
for i =  2:1000
    i
    LT = LUAD2.tab(i,:)'; %%Label - the ith mutation
    [tnet] = TrainLSTM(i,LUAD2,LT); %%Train LSTM to predict ith mutation 
    
    
    L2 = BROAD2.tab(i,:)';
    [RES_LUAD.AUC(i),sc] = getAUCLSTM(tnet,BROAD2,i,L2); %%%Test on test set
    RES_LUAD.sc(:,i)=sc(:,2);
    RES_LUAD.L(:,i) = L2;


end