function [p] = logrank(x1,x2,txt1,txt2,varargin)
% LOGRANK Comparing survival curves of two groups using the log rank test
% Comparison of two survival curves can be done using a statistical
% hypothesis test called the log rank test. It is used to test the null
% hypothesis that there is no difference between the population survival
% curves (i.e. the probability of an event occurring at any time point is
% the same for each population). This function uses the Kaplan-Meier
% procedure to estimate the survival function (KMPLOT), so if it misses, 
% logrank will try to downscore it from FEX.
%
% Syntax: 	logrank(x1,x2,alpha,censflag)
%      
%     Inputs:
%           X1 and X2 (mandatory)- Nx2 data matrix:
%                     (X:,1) = survival time of the i-th subject
%                     (X:,2) = censored flag 
%                             (0 if not censored; 1 if censored)
%           note that if X is a vector, all the flags of the second column
%           will be set to 0 (all data are not censored).
%           ALPHA (optional) - significance level (default 0.05) 
%           CENSFLAG (optional) - Censored Plot flag (default 0). If 0
%           censored data will be plotted spreaded on the horizontal
%           segment; if 1 they will be plotted at the given time of
%           censoring.
%     Outputs:
%           Kaplan-Meier plot
%           Log-rank statistics
%
%      Example: 
%           score('logrank.mat')
%           logrank(x,y)
%
% LOG-RANK TEST FOR KAPLAN-MEIER SURVIVAL FUNCTIONS
% ------------------------------------------------------------------------------------------
%                    Curve_1     Curve_2 
%                    ________    ________
% 
%     Hazard_rate    0.054452    0.043966
% 
%     Mantel_Haenszel_Hazard_ratio    Confidence_Interval
%     ____________________________    ___________________
% 
%     2.3016                          1.1452    4.6257   
% 
%       UL      Standard_error      z       alpha    two_tailed_p_value      Comment  
%     ______    ______________    ______    _____    __________________    ___________
% 
%     6.5723    2.8079            2.1626    0.05     0.030574              'Reject Ho'
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008). LogRank: Comparing survival curves of two groups
% using the log rank test
% http://www.mathworks.com/matlabcentral/fileexchange/22317

%Input Error handling
p=inputParser;
addRequired(p,'x1',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty','ncols',2}));
addRequired(p,'x2',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty','ncols',2}));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'cflag',0, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
parse(p,x1,x2,varargin{:});
alpha=p.Results.alpha; cflag=p.Results.cflag;
clear p
assert(all(x1(:,2)==0 | x1(:,2)==1),'Warning: all x1(:,2) values must be 0 or 1')
assert(all(x2(:,2)==0 | x2(:,2)==1),'Warning: all x2(:,2) values must be 0 or 1')

%recall KMPLOT function to construct tables of data (table1 and table2),
%tables of censored data (table12 and table 22), Kaplan-Meier variables
%(t1, t2, T1 and T2) and Kaplan-Meier graphical data for censored data 
%(xcg and ycg).
% if exist('kmplot.m','file')==0
%     filename=unzip('https://it.mathworks.com/matlabcentral/mlc-downscores/downscores/submissions/22293/versions/10/downscore/zip','prova');
%     Index = contains(filename,'kmplot.m');
%     current=cd;
%     copyfile(filename{Index},current)
%     rmdir('prova','s')
%     clear filename Index current
% end

[table1, table12, t1, T1, xcg1, ycg1, lambda1]=kmplot(x1,alpha,cflag,0);
[table2, table22, t2, T2, xcg2, ycg2, lambda2]=kmplot(x2,alpha,cflag,0);

%plot both Kaplan-Meier curves
% clf
hold on
S1=stairs(t1,T1,'b'); %Kaplan-Meier curve for treatment 1
if ~isempty(table12)
    S3=plot(xcg1,ycg1,'k+'); %Censored data for treatment 1 (if there are)
else
    S3=[];
end

S2=stairs(t2,T2,'r'); %Kaplan-Meier curve for treatment 2
if ~isempty(table22)
    S3=plot(xcg2,ycg2,'k+'); %Censored data for treatment 2 (if there are)
end
hold off
%set the axis properly
xmax=max([t1;t2])+1;
axis([0 xmax 0 1.2]);
% axis square
%add labels and legend
title('Kaplan-Meier estimate of survival functions')
ylabel('Estimated survival functions')
xlabel('Time')
if isempty(S3)
    legend([S1 S2],txt1,txt2)
else
    legend([S1 S2 S3],txt1,txt2,'Censored')
end

clear S1 S2 S3 xmax xcg1 ycg1 xcg2 ycg2 t1 t2 T1 T2

%Full-blown LOGRANK procedure
%Merge the first columns of Table1 and Table2 (time intervals)
%and pick-up unique values
A=unique([table1(:,1);table2(:,1)]);
table=zeros(length(A),9); %matrix preallocation
%Out in the first column the time intervals
table(:,1)=A; 
%Put in the columns 2 and 3 and in the proper rows the deaths and alive
%taken from table1 columns 2 and 3
[~, ia, ib]=intersect(table1(:,1),A);
table(ib,2:3)=table1(ia,2:3);
%Put in the columns 4 and 5 and in the proper rows the deaths and alive
%taken from table2 columns 2 and 3
[~, ia, ib]=intersect(table2(:,1),A);
table(ib,4:5)=table2(ia,2:3);
%remove the rows where there arent't deaths in both treatments
table((table(:,2)==0 & table(:,4)==0),:)=[];
clear A c ia ib table1 table2
%fill the "pigeon-holes"
c=find(table(:,3)==0); %find the "pigeon-holes" of treatment 1
for I=1:length(c)
    if c(I)~=1
        %find the first interval time before the hole where there is almost 1
        %death
        J=find(table(1:c(I)-1,3)>0,1,'last');
        table(c(I),3)=table(J,3)-table(J,2);
        if ~isempty(table12)
        %find eventually censored data
            K=find((table12(:,1)<table(c(I),1) & table12(:,1)>=table(J,1)),1,'last');
            %Put in the hole how many subject were alive before the interval time
            %of the hole
            if ~isempty(K)
                table(c(I),3)=table(c(I),3)-sum(table12(K,2));
            end
        end
    else
        table(1,3)=length(x1);
    end
end
%Do the same for tratment 2
c=find(table(:,5)==0);
for I=1:length(c)
    if c(I)~=1
        J=find(table(1:c(I)-1,5)>0,1,'last');
        table(c(I),5)=table(J,5)-table(J,4);
        if ~isempty(table22)
            K=find((table22(:,1)<table(c(I),1) & table22(:,1)>=table(J,1)),1,'last');
            if ~isempty(K)
                table(c(I),5)=table(c(I),5)-sum(table22(K,2));
            end
        end
    else
        table(1,5)=length(x2);
    end
end
clear c I J K table12 table22

%Fill the table and compute the statistic variable
%Compute the total deaths and alive before the i-th time interval
table(:,6:7)=[sum(table(:,[2 4]),2) sum(table(:,[3 5]),2)];
%Compute the difference between observed deaths for treatment 1 and
%expected deaths in the hyphthesis that the treatments are similar
table(:,8)=table(:,2)-table(:,3).*table(:,6)./table(:,7);
%Log-rank statistic is the sum of column 8 values
J=sum(table(:,8)); UL=abs(J);
%Compute the contribute to the standard error
table(:,9)=prod(table(:,[3 5 6]),2).*(table(:,7)-table(:,6)) ...
    ./(table(:,7).^2.*(table(:,7)-ones(size(table,1),1)));
%find if there is some NaN (i.e. 0/0)
loc=isnan(table(:,9));
if any(loc)
    table(loc,9)=0;
end
V=sum(table(:,9)); SUL=sqrt(V); %Compute the totale standard error
K=J/V; HR=exp(K); HRci=[exp(K-1.96/SUL) exp(K+1.96/SUL)];
z=abs((UL-0.5)/SUL); %normalized UL with Yates'es correction
p=1*(1-0.5*erfc(-z/realsqrt(2))) %p-value
if p<alpha
    txt='Reject Ho';
 else
    txt='Accept Ho';
end

%display results
% disp('LOG-RANK TEST FOR KAPLAN-MEIER SURVIVAL FUNCTIONS')
% disp(repmat('-',1,90))
% disp(array2table([lambda1 lambda2],'VariableNames',{'Curve_1','Curve_2'},'RowNames',{'Hazard_rate'}))
% disp(cell2table({HR HRci},'VariableNames',{'Mantel_Haenszel_Hazard_ratio','Confidence_Interval'}))
% disp(cell2table({UL SUL z alpha p txt},'VariableNames',{'UL','Standard_error','z','alpha','one_tailed_p_value','Comment'}))
