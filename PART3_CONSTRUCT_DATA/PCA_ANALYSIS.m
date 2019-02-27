
%%%%%%%Perform PCA analysis and plot results


%%%Load colon results
load('SIMCOAD.mat')
load('DALL_colon.mat')
load('DFCI2.mat')


%%%Plot Colon PCA
subplot(1,2,1)

X1 = SIMCOAD;
X2 = DALL_colon.tab;
XX=[X1,X2];
[coeff,score,latent] = pca(XX');

PC1 = score(:,1);
PC2 = score(:,2);
PC3 = score(:,3);

scatter3(PC1(1:253),PC2(1:253),PC3(1:253),90,'filled')
hold on
scatter3(PC1(253+1:253+72),PC2(253+1:253+72),PC3(253+1:253+72),90,'filled')
hold on
scatter3(PC1(253+72+1:253+72+619),PC2(253+72+1:253+72+619),PC3(253+72+1:253+72+619),90,'filled')
hold on
scatter3(PC1(253+72+619+1:253+72+619+100),PC2(253+72+619+1:253+72+619+100),PC3(253+72+619+1:253+72+619+100),90,'filled')

legend('COAD n = 245','MGI n = 72','DFCI n = 619','Simultated n = 100')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')


%%%Load lung results
load('SIMLUAD.mat')
load('DALL_lung.mat')


%%%Plot lung PCA
subplot(1,2,2)
X1 = SIMLUAD;
X2 = DALL_lung.tab;
XX=[X1,X2];
[coeff,score,latent] = pca(XX');

PC1 = score(:,1);
PC2 = score(:,2);
PC3 = score(:,3);

idx = kmeans(XX',2);

scatter3(PC1(1:506),PC2(1:506),PC3(1:506),90,'filled')
hold on
scatter3(PC1(506+1:506+183),PC2(506+1:506+183),PC3(506+1:506+183),90,'filled')

scatter3(PC1(506+183+1:506+183+100),PC2(506+183+1:506+183+100),PC3(506+183+1:506+183+100),90,'filled')
hold on

legend('LUAD n = 245','BROAD n = 72','Simultated n = 100')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')