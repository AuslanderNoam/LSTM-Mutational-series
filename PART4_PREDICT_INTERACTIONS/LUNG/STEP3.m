%%%% STEP3 - perform enrichment analysis for predicted interactions with
%%%% STRING interactions


%%Load dataset
load('LUADd.mat')
load('LUAD2.mat')
load('RESF.mat')
load('STRING.mat')

drg = LUADd.gene;

%%%% Get string interactions table for considered genes
x = find(sum(RESF.AUC')>0);
[TABIN,numI] = GetTableString(drg,LUAD2.gene(x),STRING);


TABIN=TABIN';

%%% Get discrete predicted interactions table between drivers and passangers
[TAB4] = GetTableINTER(RESF.SS1,RESF.ggg,TABIN);
s4=sum(TAB4>0);
TAB4(:,s4<15)=0;

[n,m] = size(RESF.AUC);


%%% Perform hyper-geomentric enrichment between predicted and STRING
%%% interactions
for i = 1:m
    
    d2(i)=nnz(TABIN(:,i));
    d1(i)=nnz(TABIN(:,i)&(TAB4(:,i)));
    d3(i)=nnz(TAB4(:,i));
    HPV3(i) = 1-hygecdf(d1(i),length(LUAD2.gene),d2(i),d3(i));
    
end
 

HPV3(d1==0|d2==0) = 1;
t1 = d2;
t2 = d3;
t3 = d1;


%%% Plot enrichment (Figure 3)

load('mapO.mat')
load('mapG.mat')

subplot(4,1,1)
h1 = heatmap(LUADd.gene,'STRING',t1)
h1.Colormap = mapG;
caxis([0,100])

subplot(4,1,2)
h1 = heatmap(LUADd.gene,'PREDICTED',t2)
h1.Colormap = mapG;
caxis([0,100])

subplot(4,1,3)
h1 = heatmap(LUADd.gene,'INTERSECT',t3)
h1.Colormap = mapG;
caxis([0,10])


subplot(4,1,4)
h1 = heatmap(LUADd.gene,'-Log(P)', -log(HPV3))
h1.Colormap = mapO;
caxis([0,2.9])


