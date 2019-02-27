%%%%PLOT figure1A-D

%%%Load results structures
load('LUAD_RES.mat')
load('COAD_RES.mat')



%%%Figure 1A
subplot(2,6,1:6)
plot(COAD_RES.AUC1)
hold on 
plot(COAD_RES.AUC2)
hold on 
plot(LUAD_RES.AUC1)
ylabel('AUC')
legend('COAD: Geneintech 2012, n = 72','COAD: DFCI cell reports 2016, n = 619','LUAD: BROAD Cell 2018, n = 183')
xlim([0,101])



%%%Figure 1B
subplot(2,6,7:8)
scatter(COAD_RES.sc1(:,20),log(COAD_RES.sumut1),'filled','MarkerFaceColor',[0 .9 .6])
legend(['Rho = ',num2str(COAD_RES.rh1(20)),'P = ',num2str(COAD_RES.p1(20))])
xlabel('Assigned score')
ylabel('Log Mutation count')

p = polyfit(COAD_RES.sc1(:,20),log(COAD_RES.sumut1),1); 
v = polyval(p,COAD_RES.sc1(:,20));
hold on
plot(COAD_RES.sc1(:,20),v,'r');
 
    
%%%Figure 1C
subplot(2,6,9:10)
scatter(COAD_RES.sc2(:,20),log(COAD_RES.sumut2),'filled','MarkerFaceColor',[0 .9 .6])
legend(['Rho = ',num2str(COAD_RES.rh2(20)),'P = ',num2str(COAD_RES.p2(20))])
xlabel('Assigned score')
ylabel('Log Mutation count')

p = polyfit(COAD_RES.sc2(:,20),log(COAD_RES.sumut2),1); 
v = polyval(p,COAD_RES.sc2(:,20));
hold on
plot(COAD_RES.sc2(:,20),v,'r');
 

%%%Figure 1D
subplot(2,6,11:12)
scatter(LUAD_RES.sc1(:,20),log(LUAD_RES.sumut1),'filled','MarkerFaceColor',[0 .9 .6])
legend(['Rho = ',num2str(LUAD_RES.rh1(20)),'P = ',num2str(LUAD_RES.p1(20))])
xlabel('Assigned score')
ylabel('Log Mutation count')

p = polyfit(LUAD_RES.sc1(:,20),log(LUAD_RES.sumut1+0.01),1); 
v = polyval(p,LUAD_RES.sc1(:,20));
hold on
plot(LUAD_RES.sc1(:,20),v,'r');

xlim([min(LUAD_RES.sc1(:,20)), max(LUAD_RES.sc1(:,20))])

