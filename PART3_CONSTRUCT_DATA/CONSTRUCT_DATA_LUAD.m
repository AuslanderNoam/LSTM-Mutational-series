% % % CONSTRUCT_DATA_LUAD - constructing simulated data for lung cancer
% one step (mutation) at a time

%%%Load dataset structures:
load('LUAD2.mat')
load('BROAD2.mat')

%%%Combine datasets together
DALL.gene = LUAD2.gene;
DALL.sample = [LUAD2.sample;BROAD2.sample];
DALL.tab = [LUAD2.tab,BROAD2.tab];

S300 = [zeros(55,1);ones(45,1)]'; %%By the frequency of the first mutation
MPRC = sum(DALL.tab')/length(DALL.sample);%% Get all frequencies

%%%%First part: reconstructing mutation 1-300 (using all previous
%%%%mutations).Generates S300-the simulated data of the first 300 mutations
NM = MPRC*100;
for i = 2:300
    i
    
    LT = DALL.tab(i,:)';
    [tnetSeq,sm] = TrainSeq(i,DALL,LT);
    for j = 1:100
        INS = S300(1:i-1,j)';
        [YPred(j),scores] = classify(tnetSeq,INS);
        sc(j) = scores(2);
    end
    
    [sv,si] = sort(sc,'descend');
    MI = zeros(1,100);
    MI(si(1:ceil(NM(i)))) = 1;
    
    S300 = [S300;MI];
end


%%%%Second part: reconstructing mutation 301-end (using a tail of 100
%%%%previos mutations). Generates SIMCOAD-the full simulated mutational data
Tail = 100;%%% Flexible to changed, increases running time


SIMCOAD = S300;
MPRC = sum(DALL.tab')/length(DALL.sample);
NM = MPRC*100;
for i = 301:length(DALL.gene)
    i
    
    LT = DALL.tab(i,:)';
    [tnetSeq,sm] = TrainSeq2(i,DALL,LT);
    for j = 1:100
        INS = SIMCOAD(i-1-Tail:i-1,j)';
        [YPred(j),scores] = classify(tnetSeq,INS);
        sc(j) = scores(2);
    end
    
    [sv,si] = sort(sc,'descend');
    MI = zeros(1,100);
    MI(si(1:ceil(NM(i)))) = 1;
    
    SIMCOAD = [SIMCOAD;MI];
end

save('SIMLUAD.mat','SIMLUAD')