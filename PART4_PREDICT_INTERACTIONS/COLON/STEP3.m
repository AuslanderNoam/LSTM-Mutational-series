%%%% STEP3 - perform enrichment analysis for predicted interactions with
%%%% STRING interactions


%%Load dataset
load('COAD2.mat')
load('COADd.mat')
load('RESF.mat')
load('STRING.mat')
drg = COADd.gene;

%%%% Get string interactions table for considered genes
x = find(sum(RESF.AUC')>0);
[TABIN] = GetTableString(drg,COAD2.gene(x),STRING);
TABIN=TABIN';

%%% Get discrete predicted interactions table between drivers and passangers
[TAB4] = GetTableINTER(RESF.SS1,RESF.ggg,TABIN);
s4 = sum(TAB4>0);
TAB4(:,s4<=15)=0;


[n,m] = size(RESF.AUC);


%%% Perform hyper-geomentric enrichment between predicted and STRING
%%% interactions
d1=[];d2=[];d3=[];HPV3=[];
for i = 1:m
    
    d2(i)=nnz(TABIN(:,i));
    d1(i)=nnz(TABIN(:,i)&(TAB4(:,i)));
    d3(i)=nnz(TAB4(:,i));
    HPV3(i) = 1-hygecdf(d1(i),length(COAD2.gene),d2(i),d3(i));
    
end
 


%%% Plot enrichment (Figure 3)
HPV3(d1==0|d2==0) = 1;
t1 = d2;
t2 = d3;
t3 = d1;

load('mapR.mat')
load('mapB.mat')

subplot(4,1,1)
h1 = heatmap(COADd.gene,'STRING',t1)
h1.Colormap = mapB;
caxis([0,100])

subplot(4,1,2)
h1 = heatmap(COADd.gene,'PREDICTED',t2)
h1.Colormap = mapB;
caxis([0,100])

subplot(4,1,3)
h1 = heatmap(COADd.gene,'INTERSECT',t3)
h1.Colormap = mapB;
caxis([0,10])


subplot(4,1,4)
h1 = heatmap(COADd.gene,'-Log(P)', -log(HPV3))
h1.Colormap = mapR;
caxis([0,2.9])


