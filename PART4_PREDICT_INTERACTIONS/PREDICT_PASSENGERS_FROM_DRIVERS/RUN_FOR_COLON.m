%%%%%Predict passengers from drivers for colon cancer
%%%%%For each driver genes (i=1-42) train and test LSTMs predicting all
%%%%%possible passenger genes from the sequence of drivers up to i


for i = 1:42
    RUNPREDPASS_COLON(i);
end

load('COADd.mat')
load('COAD2.mat')

RESF.SC2=[];RESF.SC3=[];
RESF=struct;
for i = 1:42
    load(['RES',num2str(i),'.mat'])
    s1=sum(RES.AUC2);
    id = find(s1);
    RESF.AUC2(:,i) = RES.AUC2(:,id(1));
    RESF.AUC3(:,i) = RES.AUC3(:,id(1));
    RESF.SC2(:,i,1:length(DFCI2.sample)) = RES.SC2(:,id(1),:);
    RESF.SC3(:,i,1:length(MGI2.sample)) = RES.SC3(:,id(1),:);
end