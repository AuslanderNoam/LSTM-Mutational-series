%%%% STEP4 - survival analysis for modules of drivers and their predicted
%%%% interactions


load('RESF.mat')
load('LUAD1.mat')
load('TAB4.mat')
load('LUAD_clin.mat');
load('BROAD_CLIN.mat')
load('LUAD2.mat')
load('LUADd.mat')
load('BROAD1.mat')
load('BROAD2.mat')
load('BROADd.mat')

LUAD1.sample2=strrep(LUAD1.sample,'-01','');
[a,b,c] = intersect(LUAD_clin.sample,LUAD1.sample2);


%%%Orginize survival structures
SURV.survival = LUAD_clin.survival(b);
SURV.death = LUAD_clin.death(b);
SURV.sample = a;
SURV.tab = LUAD1.tab(:,c);
SURV.gene = LUAD1.gene;


[a,b,c] = intersect(BROAD_CLIN.sample,BROAD1.sample);
BSURV.survival = BROAD_CLIN.PFS_months(b);
BSURV.survival(isnan(BSURV.survival)) = 0;

BSURV.death = BROAD_CLIN.cencored(b);
BSURV.death(isnan(BSURV.death)) = 0;

BSURV.sample = a;
BSURV.tab = BROAD1.tab(:,c);
BSURV.gene = BROAD1.gene;






txt1 = 'high signature score';
txt2 = 'low signature score';


%%%Driver module 1
driveSig1 = [ {'DNMT3A'} {'BRAF'}    {'LPHN2'}  {'TP53'}   {'NF1'} {'STK11' }]

%%%Driver module 2
driveSig2 = [ {'ATM'}       {'MMP2'}      {'SMARCA4'} {'SMAD4'}     ]


%%%Plot survival for Driver module 1
[P1,IDS1,T2] = PlotSurv1(SURV,RESF,LUAD2,LUADd,TAB4,driveSig1);
[P11,IDS11,T22] = PlotSurv1(BSURV,RESF,BROAD2,BROADd,TAB4,driveSig1);


[P2,PASS1] = PlotSurv2(SURV,LUAD2,LUADd,driveSig1,IDS1,T2);
% 
% % 
figure()
subplot(1,2,1)
[P3] = PlotSurv3(SURV,LUAD2,IDS1,T2);
text(1,1,['Log-rank P = ',num2str(P3)])
subplot(1,2,2)
[P33] = PlotSurv3(BSURV,BROAD2,IDS11,T22)

text(1,1,['Log-rank P = ',num2str(P33)])



% 
%%%Plot survival for Driver module 2
[P1,IDS1,T2] = PlotSurv1(SURV,RESF,LUAD2,LUADd,TAB4,driveSig2);
[P11,IDS11,T22] = PlotSurv1(BSURV,RESF,BROAD2,BROADd,TAB4,driveSig2);

[P22,PASS2] = PlotSurv2(SURV,LUAD2,LUADd,driveSig2,IDS1,T2);

figure()
subplot(1,2,1)
[P31] = PlotSurv3(SURV,LUAD2,IDS1,T2);
text(1,1,['Log-rank P = ',num2str(P31)])
subplot(1,2,2)
[P331] = PlotSurv3(BSURV,BROAD2,IDS11,T22)
text(1,1,['Log-rank P = ',num2str(P331)])