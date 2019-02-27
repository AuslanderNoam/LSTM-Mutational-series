%%%%%Predict passengers from drivers for lung cancer
%%%%%For each driver genes (i=1-26) train and test LSTMs predicting all
%%%%%possible passenger genes from the sequence of drivers up to i


for i = 1:26
    RUNPREDPASS_LUNG(i);
end

load('BROAD2.mat')


RESF.SC=[];RESF.AUC=[];
RESF=struct;
for i = 1:26
    load(['RES',num2str(i),'.mat'])
    s1=sum(RES.AUC);
    id = find(s1);
    RESF.AUC(:,i) = RES.AUC(:,id(1));
    RESF.SC(:,i,1:length(BROAD2.sample)) = RES.SC(:,id(1),:);
end