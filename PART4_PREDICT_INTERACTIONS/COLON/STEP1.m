%%%% STEP1 - generate tables of AUCs connecting between drivers and predicted passangers, 
%%%% and identify passangets that are repeatedly predicted from drivers genes



%%%% Load data structures
load('COAD2.mat')
load('COADd.mat')
load('RESF.mat')
load('STRING.mat')

drg = COADd.gene;

x = find(sum(RESF.AUC')>0);


%%% get table of AUCs ordered by considered genes
[TABSC3] = GetTableAUC2(RESF.AUC,x,0.925);



s1=sum(COAD2.tab(x,:)');
%%% remove genes that are very infrequently mutated
TABSC3(s1<3,:) = 0;

ids = 1;
[n,m] = size(RESF.AUC);


RESF.ss1=sum(TABSC3');
RESF.ss2=sum(COAD2.tab');

%%% passenger genes to consider
RESF.ggg=find(RESF.ss1>3&RESF.ss2>8);
RESF.TABSC3 = TABSC3;

