% 1. Run_PRED_LOAD_LSTM
%Training LSTMs to predict the mutational load, colon&lung 


%%%%%%% Run For COAD

%%%Load data structures
load('COAD2.mat') %%training
load('MGI2.mat')  %%test1
load('DFCI2.mat') %%test2


COAD_RES = struct

for i = 2:101 
    i
    [tnet] = TrainLSTM(i,COAD2); %%%Train LSTM with i last-ordered mutations
    
    [COAD_RES.AUC1(i),sc,sumut1] = testLSTM(tnet,MGI2,i); %%%Test on test set 1
    [COAD_RES.rh1(i),COAD_RES.p1(i)] = corr(sc(:,2),sumut1);
    COAD_RES.sc1(:,i)=sc(:,2);
    COAD_RES.sumut1 = sumut1;
    
    [COAD_RES.AUC2(i),sc,sumut2] = testLSTM(tnet,DFCI2,i); %%%Test on test set 2
    [COAD_RES.rh2(i),COAD_RES.p2(i)] = corr(sc(:,2),sumut2);
    COAD_RES.sc2(:,i)=sc(:,2);
    COAD_RES.sumut2 = sumut2;
end





%%%%%%% Run For LUAD

%%Load data structures
load('BROAD2.mat')  %%training
load('LUAD2.mat')   %%test


LUAD_RES = struct
for i = 2:101
    i
    [tnet] = TrainLSTM(i,LUAD2);    %%%Train LSTM with i last-ordered mutations
  
    [LUAD_RES.AUC1(i),sc,sumut1] = testLSTM(tnet,BROAD2,i); %%%Test on test set
    [LUAD_RES.rh1(i),LUAD_RES.p1(i)] = corr(sc(:,2),sumut1,'type','Spearman');
    LUAD_RES.sc1(:,i)=sc(:,2);
    LUAD_RES.sumut1 = sumut1;
  
end
