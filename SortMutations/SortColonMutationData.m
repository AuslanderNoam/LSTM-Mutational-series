%%1. Get Order score of each mutation

P1 = [];
P2 = [];
for i = 1:length(MUTD.gene)-1
    i
    for j = i:length(MUTD.gene)
        
        T1 = MUTD.tab(i,:)&MUTD.tab(j,:); %%i occuring with j
        T2 = MUTD.tab(i,:)&~MUTD.tab(j,:);%%i occuring without j
        P1(i,j) = nnz(T1);
        P2(i,j) = nnz(T2);
        
    end
end

P11 = P1;
P22 = P2;

%%Make symetrical matrix
for i = 1:length(MUTD.gene)-1
    i
    
    P11(i+1:end,i) = P1(i,i+1:end-1)';
    P22(i+1:end,i) = P2(i,i+1:end-1)';
end


%%Calc ratios of each mutation 
I2 = P11./P22;

%%Remove Nans and Inf
I2(I2==Inf) = 0;
I2(isnan(I2)) = 0;

%%Order score is sum of ratios
ss2=sum(I2);
MUTD.score=ss2;


%%2. Sort mutational data by order score
MUTD.tab2 = double(MUTD.tab>0); %%%%DECIDE


sm = MUTD.score';

[sv,si] = sort(sm,'descend');

%%TCGA data
COAD2.gene=MUTD.gene(si);
COAD2.tab=MUTD.tab2(si,:);
COAD2.sample=MUTD.sample;

sm2 = sum(COAD2.tab>0)';
[sv2,si2] = sort(sm2);
COAD2.sample=COAD2.sample(si2);
COAD2.tab=COAD2.tab(:,si2);


%%MGI data
load('MGI.mat')
NN = 18955; %%%Number of mutations to account for
g = COAD2.gene(1:NN);




MGI2 = struct;
MGI2.sample = MGI.sample;
MGI2.gene = COAD2.gene(1:NN);
tab = [];
for i = 1:NN
    ind = find(strcmp(MGI.genes,COAD2.gene(i)));
    if ~isempty(ind)
        tab(i,:) = MGI.tab(ind,:);
    end
end
MGI2.tab = tab;
MGI2.tab(18955,:)=0;

sm3 = sum(MGI2.tab>0)';
[sv3,si3] = sort(sm3);
MGI2.sample=MGI2.sample(si3);
MGI2.tab=MGI2.tab(:,si3);


%%DFCI data
load('DFCI.mat')
g = COAD2.gene(1:NN);

DFCI2 = struct;
DFCI2.sample = DFCI.sample;
DFCI2.gene = COAD2.gene(1:NN);

for i = 1:NN
    ind = find(strcmp(DFCI.genes,COAD2.gene(i)));
    try
    DFCI2.tab(i,:) = DFCI.tab(ind,:);
    end
end

DFCI2.tab(18955,:)=0;


sm4 = sum(DFCI2.tab>0)';
[sv4,si4] = sort(sm4);
DFCI2.sample=DFCI2.sample(si4);
DFCI2.tab=DFCI2.tab(:,si4);


%%Get mutations that appear at least once in all datasets
MI = sum(DFCI2.tab')>0&sum(MGI2.tab')>0&sum(COAD2.tab')>0;

DFCI2.gene = DFCI2.gene(MI==1);
DFCI2.tab = DFCI2.tab(MI==1,:);
MGI2.gene = MGI2.gene(MI==1);
MGI2.tab = MGI2.tab(MI==1,:);
COAD2.gene = COAD2.gene(MI==1);
COAD2.tab = COAD2.tab(MI==1,:);

