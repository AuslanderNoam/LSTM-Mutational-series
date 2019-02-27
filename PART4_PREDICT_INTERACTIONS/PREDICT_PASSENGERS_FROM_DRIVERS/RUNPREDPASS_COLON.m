function [] = RUNPREDPASS_COLON(mts)
%%%%% Training LSTM to predict each passanger from driver mutations i
%%%%% (saves numbered data structures)

load('COADd.mat')
load('COAD2.mat')
load('MGI2.mat')
load('MGId.mat')
load('DFCI2.mat')
load('DFCId.mat')
RES = struct;
% mtg = str2num(mts); %%% FOR HPC running
mtg = mts;

for i = 1:length(COAD2.gene)
    i
    L1 = COAD2.tab(i,:)';
    [tnet,XTS,YTS] = TrainLSTM(mtg,COADd,L1);

    
    L3 = DFCI2.tab(i,:)';
    [RES.AUC2(i,mtg),scores] = getAUCLSTM(tnet,DFCId,mtg,L3);
    RES.SC2(i,mtg,1:length(DFCI2.sample)) = scores(:,2);
    
    L3 = MGI2.tab(i,:)';
    [RES.AUC3(i,mtg),scores] = getAUCLSTM(tnet,MGId,mtg,L3);
    RES.SC3(i,mtg,1:length(MGI2.sample)) = scores(:,2);
end

save(['RES',mts,'.mat'],'RES')