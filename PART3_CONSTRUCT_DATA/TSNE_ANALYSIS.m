

%%%%%%%Perform tSNE analysis and plot results


%%%Load colon results
load('SIMCOAD.mat')
load('DALL_colon.mat')
load('DFCI2.mat')

%%%Plot Colon tSNE
subplot(1,2,1)

X1 = SIMCOAD;
X2 = DALL_colon.tab;
XX=[X1,X2];

Y = tsne(XX');

scatter(Y((1:253),1),Y((1:253),2),90,'filled')
hold on 
scatter(Y((253+1:253+72),1),Y((253+1:253+72),2),90,'filled')
hold on
scatter(Y((253+72+1:253+72+619),1),Y((253+72+1:253+72+619),2),90,'filled')
hold on 
scatter(Y((253+72+619+1:253+72+619+100),1),Y((253+72+619+1:253+72+619+100),2),90,'filled')


legend('COAD n = 245','MGI n = 72','DFCI n = 619','Simultated n = 100')
xlabel('Dimension 1')
ylabel('Dimension 2')

%%%Load lung results
load('SIMLUAD.mat')
load('DALL_lung.mat')

subplot(1,2,2)

%%%Plot lung tSNE
X1 = SIMLUAD;
X2 = DALL_lung.tab;
XX=[X1,X2];

Y = tsne(XX');

scatter(Y((1:506),1),Y((1:506),2),90,'filled')
hold on 
scatter(Y((506+1:506+183),1),Y((506+1:506+183),2),90,'filled')
hold on
scatter(Y((506+183+1:506+183+100),1),Y((506+183+1:506+183+100),2),90,'filled')
hold on 

legend('LUAD n = 245','BROAD n = 72','Simultated n = 100')
xlabel('Dimension 1')
ylabel('Dimension 2')

