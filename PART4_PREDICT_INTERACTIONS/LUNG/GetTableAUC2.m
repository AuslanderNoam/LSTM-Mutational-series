function [TAB] = GetTableAUC2(AUC,inds,qauc)
%%%%Get AUC table of ordered by considered genes

%%%%Input: 
%AUC - test AUCs for each driver, passanger
%inds - passenger genes
%qAUC - quantile of top AUCs to consider

%%%%Output: 
%TAB - binary table of AUCs that are above thresholds



%%WORKS
% qdf = 0.7;
TAUC = 0.85;

T1 = AUC(inds,:);
TAB = zeros(size(T1));

[n,m] = size(T1);

mm = max(T1')';

TAB(:,1) = double(T1(:,1)>quantile(T1(:,1),qauc)&T1(:,1)>TAUC);

TAB(:,end) = double(T1(:,end)>quantile(T1(:,end),qauc)&T1(:,end)>TAUC);


for i = 2:m-1
    i
    
    
%     X1 = double(T1(:,i)>quantile(T1(:,i),qauc)&T1(:,i)-mm==0);
%     X2 = double(T1(:,i)>quantile(T1(:,i),qauc)&T1(:,i)-mm==0);
    X1 = double(T1(:,i)>quantile(T1(:,i),qauc)&T1(:,i)>TAUC);
    X2 = double(T1(:,i)>quantile(T1(:,i),qauc)&T1(:,i)>TAUC);

    V = X1'&X2';

    TAB(:,i) = V';
end


