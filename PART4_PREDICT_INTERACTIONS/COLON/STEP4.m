%%%% STEP4 - survival analysis for modules of drivers and their predicted
%%%% interactions


load('RESF.mat')
load('COAD1.mat')
load('TAB4.mat')
load('COAD_clin.mat');
load('COAD2.mat')
load('COADd.mat')

COAD1.sample2=strrep(COAD1.sample,'-01','');

[a,b,c] = intersect(COAD_clin.sample,COAD1.sample);

%%%Orginize survival structure
SURV.survival = COAD_clin.survival(b);

SURV.death = COAD_clin.death(b);

SURV.sample = a;

SURV.tab = COAD1.tab(:,c);

SURV.gene = COAD1.gene;


txt1 = 'high signature score'
txt2 = 'low signature score'




%%%Driver module 1
driveSig1 = [ {'ATM'}    {'ARID1A'}  {'BCOR'}       {'CHD4'}     {'CHD9'}   {'SMAD4'}    {'KRAS'}        {'LPHN2'}    {'NF1'}     {'TRIO'}    {'TP53'} ]

%%%Driver module 2
driveSig2 = [{'APC'}    {'ATRX'}       {'MAP3K4'}    {'MECOM'}    {'MED12' }  {'AXIN2' } {'CEP290'}]

%%%Driver module 3
driveSig3 = [ {'DNMT3A'} {'CREBBP'} {'CNOT1' } {'WT1'   } ]





%%%Plot survival for Driver module 1

[P1,IDS1,T2,L] = PlotSurv1(SURV,RESF,COAD2,COADd,TAB4,driveSig1);


[P21] = PlotSurv2(SURV,COAD2,COADd,driveSig1,IDS1,T2);

figure()
[P31,PASS1] = PlotSurv3(SURV,COAD2,IDS1,T2);
text(1,1,['Log-rank P = ',num2str(P31)])


%%%Plot survival for Driver module 2

[P1,IDS1,T2,L] = PlotSurv1(SURV,RESF,COAD2,COADd,TAB4,driveSig2);
% 
[P22] = PlotSurv2(SURV,COAD2,COADd,driveSig2,IDS1,T2);
% 
figure()
[P32,PASS2] = PlotSurv3(SURV,COAD2,IDS1,T2);
text(1,1,['Log-rank P = ',num2str(P32)])
% 
% 




%%%Plot survival for Driver module 3

[P1,IDS1,T2,L] = PlotSurv1(SURV,RESF,COAD2,COADd,TAB4,driveSig3);
% 

[P23] = PlotSurv2(SURV,COAD2,COADd,driveSig3,IDS1,T2); 
figure()
[P33,PASS3] = PlotSurv3(SURV,COAD2,IDS1,T2);
text(1,1,['Log-rank P = ',num2str(P33)])

