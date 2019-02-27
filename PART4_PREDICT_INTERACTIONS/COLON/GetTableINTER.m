function [INTAB] = GetTableINTER(tab,inds,TI)
%%% Get discrete predicted interactions table between drivers and passangers
%%%%Input: 
%tab - predicted interactions scores
%inds - indices of selected passengers that are interacting with major drivers
%TI - string interactions table

%%%%Output: 
%INTAB - discrete predicted interactions table between major drivers and
%considered passengers

[n,m] = size(tab);


INTAB = zeros(size(TI));

size(INTAB)
for i = 1:m
    xi = tab(:,i);
    
    INTAB(inds,i) = xi;
    
end


MM=median((INTAB(INTAB>0)))
INTAB(INTAB<=MM) = 0;