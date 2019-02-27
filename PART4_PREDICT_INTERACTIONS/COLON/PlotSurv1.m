function [P1,IDS1,T2,L] = PlotSurv1(SURV,RESF,PASSd,DRIVEd,TAB4,driveSig)
%%%%Survival analysis based on interacting/non-interacting drivers

%%%%Input: 
%SURV - survival data
%RESF - Result structure
%PASSd - passanger genes mutational data
%DRIVEd - driver genes mutational data
%TAB4 - computed table of predicted interactions
%driveSig - driver module 

%%%%Output: 
%P1 - logrank P-value
%IDS - indices of sampels with mutations in driver mudule (larger than
%median)
%T2 - predicted interaction table for drivers in module
%L - number of interactors with the module 

P1=[];
txt1 = 'Interacting mutations';
txt2 = 'Not Interacting mutations';

gn2 = DRIVEd.gene;
[a,ff] = intersect(DRIVEd.gene,driveSig);
gn2(ff) = [];

mut11 = zeros(1,length(SURV.sample));
for i = 1:length(driveSig)
    m1 = find(strcmp(SURV.gene,driveSig(i)));
	mut11 = mut11+SURV.tab(m1,:);
end



mut22 = zeros(1,length(SURV.sample));
for i = 1:length(gn2)
    m2 = find(strcmp(SURV.gene,gn2(i)));
	mut22 = mut22+SURV.tab(m2,:);
end



IDS1=mut11>quantile(mut11,0.5);
IDS2=mut22>quantile(mut22,0.5);

x1 = [SURV.survival(IDS1),SURV.death(IDS1)];
x2 = [SURV.survival(IDS2),SURV.death(IDS2)];

l1 = length(x1);
l2 = length(x2);

[P1] = logrank2(x1,x2,txt1,txt2);


T2=TAB4(:,ff);
[nn,mm] = size(T2);
ss=sum(T2'>0);

PASS = PASSd.gene(ss>=mm);
L=length(PASS);


