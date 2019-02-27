function [P3] = PlotSurv3(SURV,PASSd,IDS1,T2)
%%%%Survival analyses for samples with high mutation rate in module
%%%%for score counting numbber of interactors with the module

%%%%Input: 
%SURV - survival data
%PASSd - passanger genes mutational data
%IDS - indices of sampels with mutations in driver mudule%driveSig - driver module 
%T2 - predicted interaction table for drivers in module

%%%%Output: 
%P2 - logrank P-values
%PASS - passenger interactors

txt1 = 'Interacting mutations & mutated (P)'
txt2 = 'Interacting mutations & ~mutated (P)'


[nn,mm] = size(T2)
ss=sum(T2'>0);

PASS = PASSd.gene(ss>=mm);

mut2 = zeros(1,length(SURV.sample));
for j = 1:length(PASS)
    m1 = find(strcmp(SURV.gene,PASS(j)));
    mut2 = mut2+SURV.tab(m1,:);
end
  
IDS4=mut2>median(mut2)&IDS1;
IDS5=mut2<=median(mut2)&IDS1;

x1 = [SURV.survival(IDS4),SURV.death(IDS4)];
x2 = [SURV.survival(IDS5),SURV.death(IDS5)];

l1 = length(x1)
l2 = length(x2)

[P3] = logrank(x1,x2,txt1,txt2);
