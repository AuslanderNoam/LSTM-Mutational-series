% 3. Run_PRED_LOAD_SHUFFLE_DATA
%Training classifiers to predict the mutational load based on randonly selected K mutations, colon&lung 


%%%%%%%% Run For COAD

%%%Load data structures
load('COAD2.mat') %%training
load('MGI2.mat')  %%test1
load('DFCI2.mat') %%test2


%%%get labels
sm = sum(COAD2.tab);%%training
Label = sm>median(sm);


sm1 = sum(DFCI2.tab);%%test1
LabelT1 = sm1>median(sm1);

sm2 = sum(MGI2.tab);%%test2
LabelT2 = sm2>median(sm2);



COAD22=COAD2;
DFCI22=DFCI2;
MGI22=MGI2;

for i = 1:100
    i
    
    
    %%%Do permutation
    rp = randperm(length(COAD2.gene));
    COAD2.tab=COAD2.tab(rp,:);
    DFCI2.tab=DFCI22.tab(rp,:);
    MGI2.tab=MGI22.tab(rp,:);
    
    
    %%%Cut data (use i random mutations)
    Train = COAD2.tab(1:i,:);
    Test1 = DFCI2.tab(1:i,:);
    Test2 = MGI2.tab(1:i,:);
    
    %%%KNN 
    Mdl = fitcknn(Train',Label,'NumNeighbors',5,'Standardize',1);
    
    [label,score1,cost1] = predict(Mdl,Test1');
    [X1,Y1,T,SHUFF_COAD.AUC1KNN(i)] = perfcurve(LabelT1,score1(:,2),1);
    [SHUFF_COAD.RH1KNN(i)] = corr(sm1',score1(:,2),'type','Spearman');
     
    [label,score2,cost2] = predict(Mdl,Test2');
    [X1,Y1,T,SHUFF_COAD.AUC2KNN(i)] = perfcurve(LabelT2,score2(:,2),1);
    
    
    
    [SHUFF_COAD.RH2KNN(i)] = corr(sm2',score2(:,2),'type','Spearman');
    
    
    
    %%%SVM
    Mdl2 = fitcsvm(Train',Label);
    
    [label,score11] = predict(Mdl2,Test1');
    [X1,Y1,T,SHUFF_COAD.AUC1SVM(i)] = perfcurve(LabelT1,score11(:,2),1);
    [SHUFF_COAD.RH1SVM(i)] = corr(sm1',score11(:,2),'type','Spearman');
    
    [label,score22] = predict(Mdl2,Test2');
    [X1,Y1,T,SHUFF_COAD.AUC2SVM(i)] = perfcurve(LabelT2,score22(:,2),1);
    [SHUFF_COAD.RH2SVM(i)] = corr(sm2',score22(:,2),'type','Spearman');
    
    
    
     %%%logistic
     Mdl3 = fitclinear(Train',Label,'Learner','logistic');
    
    [label,score111] = predict(Mdl3,Test1');
    [X1,Y1,T,SHUFF_COAD.AUC1LIN(i)] = perfcurve(LabelT1,score111(:,2),1);
    [SHUFF_COAD.RH1LIN(i)] = corr(sm1',score111(:,2),'type','Spearman');
    
    [label,score222] = predict(Mdl3,Test2');
    [X1,Y1,T,SHUFF_COAD.AUC2LIN(i)] = perfcurve(LabelT2,score222(:,2),1);
    [SHUFF_COAD.RH2LIN(i)] = corr(sm2',score222(:,2),'type','Spearman');
    
    
    
    %%%LSTM
    [tnet] = TrainLSTM(i,COAD2);
    [SHUFF_COAD.AUC1LSTM(i),sc,sumut1] = testLSTM(tnet,MGI2,i);
    [SHUFF_COAD.RH1LSTM(i)] = corr(sc(:,2),sumut1);
    [SHUFF_COAD.AUC2LSTM(i),sc,sumut2] = testLSTM(tnet,DFCI2,i);
    [SHUFF_COAD.RH2LSTM(i),RES.p2(i)] = corr(sc(:,2),sumut2);
    
end


save('SHUFF_COAD.mat','SHUFF_COAD')




% 
load('LUAD2.mat')
load('BROAD2.mat')

sm = sum(LUAD2.tab);
Label = sm>median(sm);


sm = sum(BROAD2.tab);
LabelT1 = sm>median(sm);
LUAD22=LUAD2;
BROAD22=BROAD2;


for i = 1:100
    
    %%%Do permutation
    rp = randperm(length(LUAD2.gene));
    LUAD2.tab=LUAD22.tab(rp,:);
    BROAD2.tab=BROAD22.tab(rp,:);
    
    i
    
    %%Cut data (use i random mutations)
    Train = LUAD2.tab(1:i,:);
    Test1 = BROAD2.tab(1:i,:);
    
    
    %%%KNN 
    Mdl = fitcknn(Train',Label,'NumNeighbors',5,'Standardize',1);
    [label,score1,cost1] = predict(Mdl,Test1');
    [X1,Y1,T,SHUFF_LUAD.AUCKNN(i)] = perfcurve(LabelT1,score1(:,2),1);
    [SHUFF_LUAD.RHKNN(i)] = corr(sm',score1(:,2),'type','Spearman');
  
    %%%SVM
    Mdl2 = fitcsvm(Train',Label);
    [label,score11,cost1] = predict(Mdl2,Test1');
    [X1,Y1,T,SHUFF_LUAD.AUCSVN(i)] = perfcurve(LabelT1,score11(:,2),1);
    [SHUFF_LUAD.RHSVM(i)] = corr(sm',score11(:,2),'type','Spearman');
   
    %%%logistic
    Mdl3 = fitclinear(Train',Label,'Learner','logistic');
    [label,score111] = predict(Mdl3,Test1');
    [X1,Y1,T,SHUFF_LUAD.AUCLIN(i)] = perfcurve(LabelT1,score111(:,2),1);
    [SHUFF_LUAD.RHLIN(i)] = corr(sm',score111(:,2),'type','Spearman');
    
    %%%LSTM
    [tnet] = TrainLSTM(i,LUAD2);
    [SHUFF_LUAD.AUCLSTM(i),sc,sumut1] = testLSTM(tnet,BROAD2,i);
    [SHUFF_LUAD.RHLSTM(i)] = corr(sc(:,2),sumut1,'type','Spearman');
    
end

save('SHUFF_LUAD.mat','SHUFF_LUAD')


