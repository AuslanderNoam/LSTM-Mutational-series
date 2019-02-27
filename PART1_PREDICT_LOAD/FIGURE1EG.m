%%%%%Plot Figure 1E-G


load('COAD_RES.mat')
load('LINCLASS_COAD.mat')
load('LUAD_RES.mat')
load('LINCLASS_LUAD.mat')
load('SHUFF_COAD.mat')
load('SHUFF_LUAD.mat')

cnt = 1;
rng = 0.5:0.01:1;
RNN1=[];RNN2=[];SVM1=[];SVM2=[];KNN1=[];KNN2=[];SVM3=[];KNN3=[];KNN3=[];
RNN11=[];RNN22=[];SVM11=[];SVM22=[];KNN11=[];KNN22=[];SVM33=[];RNN33=[];KNN33=[];



D1 = [COAD_RES.rh1(1:50);LINCLASS_COAD.RH1SVM(1:50);LINCLASS_COAD.RH1KNN(1:50);LINCLASS_COAD.RH1LIN(1:50);SHUFF_COAD.RH1LSTM(1:50);SHUFF_COAD.RH1SVM(1:50);SHUFF_COAD.RH1KNN(1:50);SHUFF_COAD.RH1LIN(1:50)];
D2 = [COAD_RES.rh2(1:50);LINCLASS_COAD.RH2SVM(1:50);LINCLASS_COAD.RH2KNN(1:50);LINCLASS_COAD.RH2LIN(1:50);SHUFF_COAD.RH2LSTM(1:50);SHUFF_COAD.RH2SVM(1:50);SHUFF_COAD.RH2KNN(1:50);SHUFF_COAD.RH2LIN(1:50)];
D3 = [LUAD_RES.rh1(1:50);LINCLASS_LUAD.RHSVM(1:50);LINCLASS_LUAD.RHKNN(1:50);LINCLASS_LUAD.RHLIN(1:50);SHUFF_LUAD.RHLSTM(1:50);SHUFF_LUAD.RHSVM(1:50);SHUFF_LUAD.RHKNN(1:50);SHUFF_LUAD.RHLIN(1:50)];

%%%%%1E

subplot(1,3,1)
plot(D1')
legend('LSTM','SVM','KNN','Logistic','LSTM shuffle','SVM shuffle','KNN shuffle','Logistic shuffle')
xlabel('Number of mutations')
ylabel('Spearman Rho')
xlim([2,50])

%%%%%1F

subplot(1,3,2)
plot(D2')
xlabel('Number of mutations')
ylabel('Spearman Rho')
xlim([2,50])


%%%%%1G

subplot(1,3,3)
plot(D3')
xlabel('Number of mutations')
ylabel('Spearman Rho')
xlim([2,50])
legend('LSTM','SVM','KNN','Logistic','LSTM shuffle','SVM shuffle','KNN shuffle','Logistic shuffle')
