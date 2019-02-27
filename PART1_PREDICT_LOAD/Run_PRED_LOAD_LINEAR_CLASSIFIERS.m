% 2. Run_PRED_LOAD_LINEAR_CLASSIFIERS
%Training linear classifiers to predict the mutational load, colon&lung 

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


for i = 2:101
    i
    
    %%%Cut data (i first ordered mutations)
    Train = COAD2.tab(1:i,:);
    Test1 = DFCI2.tab(1:i,:);
    Test2 = MGI2.tab(1:i,:);
    
    
    %%%KNN 
    Mdl = fitcknn(Train',Label,'NumNeighbors',5,'Standardize',1);
    [label,score1,cost1] = predict(Mdl,Test1');
    [X1,Y1,T,LINCLASS_COAD.AUC1KNN(i)] = perfcurve(LabelT1,score1(:,2),1);
    [LINCLASS_COAD.RH1KNN(i)] = corr(sm1',score1(:,2),'type','Spearman');
     
    [label,score2,cost2] = predict(Mdl,Test2');
    [X1,Y1,T,LINCLASS_COAD.AUC2KNN(i)] = perfcurve(LabelT2,score2(:,2),1);
    [LINCLASS_COAD.RH2KNN(i)] = corr(sm2',score2(:,2),'type','Spearman');
    
    
    
    %%%SVM
    Mdl2 = fitcsvm(Train',Label);
    [label,score11] = predict(Mdl2,Test1');
    [X1,Y1,T,LINCLASS_COAD.AUC1SVM(i)] = perfcurve(LabelT1,score11(:,2),1);
    [LINCLASS_COAD.RH1SVM(i)] = corr(sm1',score11(:,2),'type','Spearman');
    
    [label,score22] = predict(Mdl2,Test2');
    [X1,Y1,T,LINCLASS_COAD.AUC2SVM(i)] = perfcurve(LabelT2,score22(:,2),1);
    [LINCLASS_COAD.RH2SVM(i)] = corr(sm2',score22(:,2),'type','Spearman');
    
    
    
    %%%logistic
     Mdl3 = fitclinear(Train',Label,'Learner','logistic');
    
    [label,score111] = predict(Mdl3,Test1');
    [X1,Y1,T,LINCLASS_COAD.AUC1LIN(i)] = perfcurve(LabelT1,score111(:,2),1);
    [LINCLASS_COAD.RH1LIN(i)] = corr(sm1',score111(:,2),'type','Spearman');
    
    [label,score222] = predict(Mdl3,Test2');
    [X1,Y1,T,LINCLASS_COAD.AUC2LIN(i)] = perfcurve(LabelT2,score222(:,2),1);
    [LINCLASS_COAD.RH2LIN(i)] = corr(sm2',score222(:,2),'type','Spearman');
end


save('LINCLASS_COAD.mat','LINCLASS_COAD')





%%%%%%%% Run For LUAD

%%%Load data structures
load('BROAD2.mat')  %%training
load('LUAD2.mat')   %%test

sm = sum(LUAD2.tab);
Label = sm>median(sm);


sm = sum(BROAD2.tab);
LabelT1 = sm>median(sm);



for i = 2:101
    i
    
    %%%Cut data (i first ordered mutations)
    Train = LUAD2.tab(1:i,:);
    Test1 = BROAD2.tab(1:i,:);
    
    
    %%%KNN 
    Mdl = fitcknn(Train',Label,'NumNeighbors',5,'Standardize',1);
    [label,score1,cost1] = predict(Mdl,Test1');
    [X1,Y1,T,LINCLASS_LUAD.AUCKNN(i)] = perfcurve(LabelT1,score1(:,2),1);
    [LINCLASS_LUAD.RHKNN(i)] = corr(sm',score1(:,2),'type','Spearman');
  
    %%%SVM 
    Mdl2 = fitcsvm(Train',Label);
    [label,score11,cost1] = predict(Mdl2,Test1');
    [X1,Y1,T,LINCLASS_LUAD.AUCSVN(i)] = perfcurve(LabelT1,score11(:,2),1);
    [LINCLASS_LUAD.RHSVM(i)] = corr(sm',score11(:,2),'type','Spearman');
   
    %%%logistic
    Mdl3 = fitclinear(Train',Label,'Learner','logistic');
    [label,score111] = predict(Mdl3,Test1');
    [X1,Y1,T,LINCLASS_LUAD.AUCLIN(i)] = perfcurve(LabelT1,score111(:,2),1);
    [LINCLASS_LUAD.RHLIN(i)] = corr(sm',score111(:,2),'type','Spearman');
end


save('LINCLASS_LUAD.mat','LINCLASS_LUAD')