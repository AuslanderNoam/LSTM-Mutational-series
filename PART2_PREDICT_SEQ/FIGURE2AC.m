
%%%%%%%PLOT Figure 2A-C


%%%Load results
load('RES_COAD.mat')
load('RES_LUAD.mat')
load('COAD2.mat')
load('LUADDRIVERS.mat')


%%Fig 2A
subplot(1,5,1:3)
histogram(RES_COAD.AUC1,100)
hold on
histogram(RES_COAD.AUC2,100)
hold on
histogram(RES_LUAD.AUC,100)
legend('COAD test 1','COAD test 2','LUAD test')


load('COLONDRIVERS.mat')
load('LUAD2.mat')

[a,b,c] = intersect(COAD2.gene,COLONDRIVERS);

%%Fig 2B
subplot(1,5,4)

mm=((RES_COAD.AUC1+RES_COAD.AUC2)/2);
mm2=mm;
mm2(b)=[];
X = [mm2,mm(b)];
L = [ones(size(mm2)),zeros(size(mm(b)))];
L1={};
L1(L==0) = {'Drivers'};
L1(L==1) = {'Others'};
hold on
scatter(ones(size(mm2)),mm2,30,'y')
scatter(2*ones(size(mm(b))),mm(b),30,'g')
boxplot(X,L1)
RES_LUAD.AUC(b),RES_LUAD.AUC;

[p1,h1] = ranksum(mm(b),mm2,'tail','left');
text(1,1,['P-value = ',num2str(p1)])
ylabel('AUC')

%%Fig 2C
subplot(1,5,5)

[a,b,c] = intersect(LUAD2.gene,LUADDRIVERS);

mm=RES_LUAD.AUC;
mm2=mm;
mm2(b)=[];
X = [mm2,mm(b)];
L = [ones(size(mm2)),zeros(size(mm(b)))];
L1={};
L1(L==0) = {'Drivers'};
L1(L==1) = {'Others'};
hold on
scatter(ones(size(mm2)),mm2,30,'y')
scatter(2*ones(size(mm(b))),mm(b),30,'g')
boxplot(X,L1)
RES_LUAD.AUC(b),RES_LUAD.AUC;
[p1,h1] = ranksum(mm(b),mm2,'tail','left');
text(1,1,['P-value = ',num2str(p1)])
ylabel('AUC')