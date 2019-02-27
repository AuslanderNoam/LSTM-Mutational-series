function [P2,PASS] = PlotSurv2(SURV,PASSd,DRIVEd,driveSig,IDS1,T2)
%%%%Survival analyses for samples with high mutation rate in module
%%%%for individual interactors with the module

%%%%Input: 
%SURV - survival data
%PASSd - passanger genes mutational data
%DRIVEd - driver genes mutational data
%IDS - indices of sampels with mutations in driver mudule%driveSig - driver module 
%T2 - predicted interaction table for drivers in module

%%%%Output: 
%P2 - logrank P-values

txt1 = 'Interacting mutations & mutated (P)';
txt2 = 'Interacting mutations & ~mutated (P)';

P2=[];l1=[];l2=[];
cnt = 1;qq=[];

[a,ff] = intersect(DRIVEd.gene,driveSig);
[nn,mm] = size(T2);
ss=sum(T2'>0);

PASS = PASSd.gene(ss>=mm)
hfig = figure();

P11=[];P22=[];
cnt = 1;
for j = 1:length(PASS)
%     j
    
    gi = PASS(j);

    mut2 = zeros(1,length(SURV.sample));
    m1 = find(strcmp(SURV.gene,gi));
    mut2 = SURV.tab(m1,:);

    
    IDS2=mut2>0&IDS1;
    IDS3=mut2==0&IDS1;

    x1 = [SURV.survival(IDS2),SURV.death(IDS2)];
    
    x2 = [SURV.survival(IDS3),SURV.death(IDS3)];



    l1(j) = length(x1);
    l2(j) = length(x2);
    ll(j) = l1(j)+l2(j);
    try
    [P11(j)] = logrank2(x1,x2,txt1,txt2);
    end
    if l1(j)>=1&l2(j)>=1&P11(j)<0.06
        subplot(4,3,cnt);

        [P2(cnt)] = logrank(x1,x2,txt1,txt2);
        cnt = cnt+1;
        title(gi)
    end
    
end
    