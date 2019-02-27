%%%% STEP1 - generate tables of AUCs connecting between drivers and predicted passangers, 
%%%% and identify passangets that are repeatedly predicted from drivers genes



%%%% Load data structures
load('BROAD2.mat')
load('BROADd.mat')
load('LUADd.mat')
load('LUAD2.mat')
load('RESF.mat')
load('STRING.mat')

drg = LUADd.gene;


x = find(sum(RESF.AUC')>0);


%%% get table of AUCs ordered by considered genes
[TABSC3] = GetTableAUC2(RESF.AUC,x,0.925);


s1=sum(LUAD2.tab(x,:)');
%%% remove genes that are very infrequently mutated
TABSC3(s1<3,:) = 0;

ids = 1;


[n,m] = size(RESF.AUC);

RESF.ss1=sum(TABSC3');
RESF.ss2=sum(LUAD2.tab');

%%% passenger genes to consider
RESF.ggg=find(RESF.ss1>3&RESF.ss2>8);
RESF.TABSC3 = TABSC3;

