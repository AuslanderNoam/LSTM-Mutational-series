%%%% STEP2 - for each of the passenger genes that can be repeatedly
%%%% predicted from the landscape of driver genes (RESF.ggg), correlate the scores
%%%% predicting the passenger with the mutations in each driver and
%%%% find driver mutations that are significantly correlated with the
%%%% scores predicting the passenger mutations, thus are likely to
%%%% contribute to its prediction


load('RESF.mat')
load('DFCI2.mat')

rh2=[];p2=[];RH1=[];RH2=[];szi=[];
for i = 1:length(RESF.ggg)
    i
    gi = RESF.ggg(i);
    SC2(:,:) = RESF.SC2(i,:,:);
    [n1,m1] = size(SC2);
    
    mm = find(RESF.TABSC3(RESF.ggg(i),:));
    
    szi(i) = length(mm);
    for j = 1:length(mm)
        for k = 1:n1
            
            [rh2(j,k),p2(j,k)] = corr(DFCI2.tab(k,:)', SC2(mm(j),:)');
            
        end
        
    end
    
    RH1(i,:) = sum(rh2<0&p2<0.05);
    rh2=[];p2=[];
end


SS1=[]
for i = 1:n1
    
    SS1(:,i) = RH1(:,i)./szi';
end


RESF.SS1 = SS1;

