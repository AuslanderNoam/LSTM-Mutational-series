function [] = RUNPREDPASS_LUNG(mts)
%%%%% Training LSTM to predict each passanger from driver mutations i
%%%%% (saves numbered data structures)

load('LUAD2.mat')
load('LUADd.mat')
load('BROADd.mat')
load('BROAD2.mat')
RES = struct;
% mtg = str2num(mts); %%% FOR HPC running
mtg = mts;

for i = 1:length(LUAD2.gene)
    i
    L1 = LUAD2.tab(i,:)';
    [tnet,XTS,YTS] = TrainLSTM(mtg,LUADd,L1);

    
    L3 = BROAD2.tab(i,:)';
    [RES.AUC(i,mtg),scores] = getAUCLSTM(tnet,BROADd,mtg,L3);
    RES.SC(i,mtg,1:length(BROAD2.sample)) = scores(:,2);
end

save(['RES',mts,'.mat'],'RES')