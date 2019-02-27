function [TAB,numI] = GetTableString(drive,pass,STRING)
%%%%Get string interactions table for considered genes

%%%%Input: 
%drive - driver genes names
%pass - passanger genes names
%STRING - processed data of STRING interactions

%%%%Output: 
%TAB - table of string interactions for major drivers and considered
%passengers

TAB = zeros(length(drive),length(pass));
for i = 1:length(drive)
    i1 = find(strcmp(STRING.gene1,drive(i)));
    g1 = STRING.gene2(strcmp(STRING.gene1,drive(i)));


    i2 = find(strcmp(STRING.gene2,drive(i)));
    g2 = STRING.gene1(strcmp(STRING.gene2,drive(i)));

    
    gg = [g1;g2];
    numI(i) = length(unique(gg));
    ii = [i1;i2];
    [a1,b1,c1] = intersect(gg,pass);
    if ~ isempty(a1)
        TAB(i,c1) = STRING.score(ii(b1));
    end
        
end
        