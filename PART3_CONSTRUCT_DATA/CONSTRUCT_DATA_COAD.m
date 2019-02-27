% % % CONSTRUCT_DATA_COAD - constructing simulated data for colon cancer
% one step (mutation) at a time

%%%Load dataset structures:
load('COAD2.mat')
load('DFCI2.mat')
load('MGI2.mat')


%%%Combine datasets together
DALL_colon.gene = COAD2.gene;
DALL_colon.sample = [COAD2.sample;MGI2.sample;DFCI2.sample];
DALL_colon.tab = [COAD2.tab,MGI2.tab,DFCI2.tab];




S300 = [zeros(55,1);ones(45,1)]'; %%By the frequency of the first mutation

MPRC = sum(DALL_colon.tab')/length(DALL_colon.sample); %% Get all frequencies
NM = MPRC*100;

%%%%First part: reconstructing mutation 1-300 (using all previous
%%%%mutations).Generates S300-the simulated data of the first 300 mutations
for i = 2:300
    i
    
    LT = DALL_colon.tab(i,:)';
    [tnetSeq,sm] = TrainSeq(i,DALL_colon,LT);
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

Tail = 100; %%% Flexible to changed, increases running time

DALL_colon.gene = MUTD2.gene;
DALL_colon.sample = [MUTD2.sample;MGI2.sample;DFCI2.sample];
DALL_colon.tab = [MUTD2.tab,MGI2.tab,DFCI2.tab];

RES = struct
SIMCOAD = S300;
for i = 301:length(DALL_colon.gene)
    i
    
    LT = DALL_colon.tab(i,:)';
    [tnetSeq,sm] = TrainSeq2(i,DALL_colon,LT);
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

