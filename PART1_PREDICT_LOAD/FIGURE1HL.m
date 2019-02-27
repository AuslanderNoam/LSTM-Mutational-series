%%%%PLOT figure1A-D

%%%Load results structures
load('BROAD_CLIN.mat')
load('BROAD2.mat')
load('LUAD_RES.mat')




%%%Figure 1H
subplot(1,5,1)


[a,b,c] = intersect(MGI2.sample,MGI_CLIN.sample);
MSI = MGI_CLIN.MSI(c);
score = COAD_RES.sc1(b,20);


X = [score(strcmp(MSI,'MSS'));score(strcmp(MSI,'MSI'))]
L=[zeros(nnz(strcmp(MSI,'MSS')),1);ones(nnz(strcmp(MSI,'MSI')),1)]
TL = {}
TL(L==0) = {'MSS'}
TL(L==1) = {'MSI'}

[p5,h5] = ranksum(score(strcmp(MSI,'MSS')),score(strcmp(MSI,'MSI')),'tail','left');

B5 = [X,L];
boxplot(B5(:,1),TL)
ylabel('score')
text(0,0.55,['P=',num2str(p5)])
text(1,0.5,['n=',num2str(nnz(L==0))])
text(2,0.5,['n=',num2str(nnz(L==1))])


%%%%%



load('DFCI_CLIN.mat')
load('COAD_RES.mat')
load('DFCI2.mat')
load('MGI2.mat')
load('MGI_CLIN.mat')

[a,b,c] = intersect(DFCI2.sample,DFCI_CLIN.sample);

score = COAD_RES.sc2(b,20);
grade = DFCI_CLIN.grade(c);
MSI = DFCI_CLIN.MSI(c);
CIMP = DFCI_CLIN.CIMP(c);

%%%Figure 1I

subplot(1,5,2)

X = [score(strcmp(grade,'Well-Moderate'));score(strcmp(grade,'Poor'))]
L = [zeros(nnz(strcmp(grade,'Well-Moderate')),1);ones(nnz(strcmp(grade,'Poor')),1)]
TL = {}
TL(L==0) = {'Well-Moderate'}
TL(L==1) = {'Poor-Moderate'}

[p1,h1] = ranksum(score(strcmp(grade,'Well-Moderate')),score(strcmp(grade,'Poor')),'tail','left');
% 

B2 = [X,L];
boxplot(B2(:,1),TL)
ylabel('score')
text(0,0.55,['P=',num2str(p1)])
text(1,0.5,['n=',num2str(nnz(L==0))])
text(2,0.5,['n=',num2str(nnz(L==1))])

%%%Figure 1J

subplot(1,5,3)

X = [score(strcmp(MSI,'MSS'));score(strcmp(MSI,'MSI-high'))]
L=[zeros(nnz(strcmp(MSI,'MSS')),1);ones(nnz(strcmp(MSI,'MSI-high')),1)]
TL = {}
TL(L==0) = {'MSS'}
TL(L==1) = {'MSI'}

[p2,h2] = ranksum(score(strcmp(MSI,'MSS')),score(strcmp(MSI,'MSI-high')),'tail','left');

B3 = [X,L];
boxplot(B3(:,1),TL)
ylabel('score')
text(0,0.55,['P=',num2str(p2)])
text(1,0.5,['n=',num2str(nnz(L==0))])
text(2,0.5,['n=',num2str(nnz(L==1))])

 
%%%Figure 1K

subplot(1,5,4)

X = [score(strcmp(CIMP,'CIMP-0/low'));score(strcmp(CIMP,'CIMP-high'))]
L = [zeros(nnz(strcmp(CIMP,'CIMP-0/low')),1);ones(nnz(strcmp(CIMP,'CIMP-high')),1)]
TL = {}
TL(L==0) = {'CIMP low'}
TL(L==1) = {'CIMP high'}

[p3,h3] = ranksum(score(strcmp(CIMP,'CIMP-0/low')),score(strcmp(CIMP,'CIMP-high')),'tail','left');

B4 = [X,L];
boxplot(B4(:,1),TL)
ylabel('score')
text(0,0.55,['P=',num2str(p3)])
text(1,0.5,['n=',num2str(nnz(L==0))])
text(2,0.5,['n=',num2str(nnz(L==1))])



% % % % Figure 1L


%%%Order survival data
[a,b,c] = intersect(BROAD_CLIN.sample,BROAD2.sample);

SURV.sample = a;
SURV.sc = LUAD_RES.sc1(c,:);
SURV.sc2 = sum(BROAD2.tab(:,c));

SURV.PFS = BROAD_CLIN.PFS_months(b);
SURV.cenc = BROAD_CLIN.cencored(b);
SURV.cenc(isnan(SURV.cenc)) = 1;

sc = SURV.sc(:,20);

x1=SURV.PFS(sc>=median(sc));
mean(x1(~isnan(x1)));
y1 = SURV.cenc(sc>=median(sc));


x2=SURV.PFS(sc<median(sc));
mean(x2(~isnan(x2)));
y2 = SURV.cenc(sc<median(sc));



subplot(1,5,5)
[p,h] = ranksum(x2,x1,'tail','right');


X = [x1;x2]
L = [zeros(size(x1));ones(size(x2))]
TL = {}
TL(L==0) = {'High score'}
TL(L==1) = {'Low score'}

B1 = [X,L]
boxplot(B1(:,1),TL)
ylim([-10,120])
ylabel('PFS months')
text(0.5,40,['P=',num2str(p)])

text(1,100,['n=',num2str(length(x1))])
text(2,100,['n=',num2str(length(x2))])



