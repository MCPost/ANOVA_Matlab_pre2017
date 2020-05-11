function EMM = EstMargMeans(Data, T, Info_Eff, X, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input Arguments:
% 
%  Data  -  Nxp matrix. N number of all observation. p contains one column
%           for Subjects (Sub), one column for the dependent variable (DV), 
%           and at least one column for a between-subjects (BW) or
%           within-subjects (WI) factor. Matrix has to be in long format,
%           so each row represents one distinct observation of the DV
%  T     -  table outputted from ANOVA function
%  X     -  design matrix outputted from ANOVA function
%  Info_eff - cell struct of information outputted from ANOVA function
%  Sub   -  integer indicating the column of Data containing the Subjects    
%  DV    -  integer indicating the column of Data containing the DV
%  BW    -  integer or vector of integers indicating the column of Data
%           containing all the between-subjects effects
%  WI    -  integer or vector of integers indicating the column of Data
%           containing all the within-subjects effects
%  CI_fun - (optional) character string indicating whether CI shall be
%           computed using a t-distribution ('t'; default) or a normal
%           distribution ('norm')
%  CI_alpha - number between 0 and 1 indicating the alpha values for two- 
%           sided CI (default: 0.05)
%
% Output Arguments:
%
%  EMM   - cell struct containing estimated marginal means, standard errors
%          and CI for all main effects and interactions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


list_input = varargin;
list_input = reshape(list_input,2,[])';

if(sum(strcmp(list_input(:,1),'Sub')) == 1 && sum(strcmp(list_input(:,1),'DV')) == 1)
    if(isnumeric(list_input{strcmp(list_input,'Sub'),2}) && isnumeric(list_input{strcmp(list_input,'DV'),2})... 
        && length(list_input{strcmp(list_input,'Sub'),2}) == 1 && length(list_input{strcmp(list_input,'DV'),2}) == 1)
        Sub = list_input{strcmp(list_input,'Sub'),2};
        DV = list_input{strcmp(list_input,'DV'),2};
    else
        error('Error: Sub Index and DV Index need to be single integer!')
    end
else
    error('Error: Sub Index and/or DV Index missing!');
end

if(sum(strcmp(list_input(:,1),'BW')) == 1)
    BW = list_input{strcmp(list_input(:,1), 'BW'),2};
    if(sum(strcmp(list_input(:,1),'BW_names')) == 1)
        BW_names = list_input{strcmp(list_input(:,1), 'BW_names'),2};
    end
else
    BW = [];
end
if(sum(strcmp(list_input(:,1),'WI')) == 1)
    WI = list_input{strcmp(list_input(:,1), 'WI'),2};
    if(sum(strcmp(list_input(:,1),'WI_names')) == 1)
        WI_names = list_input{strcmp(list_input(:,1), 'WI_names'),2};
    end
else
    WI = [];
end
if(isempty(WI) && isempty(BW))
    error('Error:  BW and WI both empty!')
end

if(Sub == DV || ~isempty(intersect([Sub DV],[BW WI])) || ~isempty(intersect(BW, WI)))
    error('Error: Sub, DV, WI and BW have to be independent integers or vectors of integers!')
end

if(sum(strcmp(list_input(:,1),'CI_fun')) == 1)
    if(strcmp(list_input(strcmp(list_input(:,1),'CI_fun'),2),'norm'))
        CI_fun = @(x,y) norminv(x);
    else
        CI_fun = @(x,y) tinv(x,y);
    end
else
    CI_fun = @(x,y) tinv(x,y);
end

if(sum(strcmp(list_input(:,1),'CI_alpha')) == 1)
    CI_alpha = list_input(strcmp(list_input(:,1),'CI_fun'),2);
else
    CI_alpha = 0.05;
end

BW_ind = BW - min(BW) + 1;
if(isempty(BW_ind))
    WI_ind = WI - min(WI) + 1;
    BW_lev_pr = 1;
    WI_lev_pr = cellfun(@(c) length(unique(c)), num2cell(Data(:,WI), 1));
else
    if(isempty(WI))
        WI_ind = [];
        BW_lev_pr = cellfun(@(c) length(unique(c)), num2cell(Data(:,BW), 1));
        WI_lev_pr = 1;
    else
        WI_ind = WI - min(WI) + 1 + BW_ind;
        BW_lev_pr = cellfun(@(c) length(unique(c)), num2cell(Data(:,BW), 1));
        WI_lev_pr = cellfun(@(c) length(unique(c)), num2cell(Data(:,WI), 1));
    end
end

Design = Data(:,[BW WI DV]);

D_mean = unique(Data(:,[BW WI]),'rows','stable');
Mean_Big = zeros(size(D_mean,1),1);
for i = 1:size(D_mean,1)
    Mean_Big(i) = nanmean(Design(sum(abs(bsxfun(@minus, Design(:,1:end-1), D_mean(i,:))),2) == 0,end));
end

[~,ind_se] = sortrows(D_mean,WI_ind);

comb = {};
for i = 1:size(D_mean,2)
    comb = [comb; num2cell(nchoosek(1:size(D_mean,2),i),2)];
end

ind = find(cellfun(@isempty, regexp(Info_Eff(:,1),'_Err')));
X_full = [];
ind_full = {};
for i = 2:length(ind)
    ind_full{i-1,1} = size(X_full,2) + (1:size(Info_Eff{ind(i),3},2));
    X_full = [X_full, Info_Eff{ind(i),3}];
end
X_mean = unique(X_full,'rows','stable');
X_se = X_mean(ind_se,:);

df_res = T.df(strcmp(T.Effect,'S_Err'));

e = (eye(size(X,1)) - X*((X'*X)\X'))*Data(:,DV);
e_w = reshape(e,[prod(WI_lev_pr) length(e)/prod(WI_lev_pr)])';
VCOV = e_w'*e_w*(1/df_res);

BW_design = X(:,cell2mat(Info_Eff(2:find(~cellfun(@isempty, regexp(Info_Eff(:,1),'S_Err')))-1,2)'));
BW_design(BW_design == -1) = 0;

X2 = unique([Data(:,Sub), ones(size(BW_design,1),1),BW_design],'rows','stable');
X2(:,1) = []; 
[~,R] = qr(X2);
cov_unscaled = inv(R'*R);

M = kron(VCOV,cov_unscaled);

if(isempty(BW_ind))
    E = eye(prod(WI_lev_pr));
else
    E = kron(eye(prod(WI_lev_pr)), unique([ones(size(BW_design,1),1),BW_design],'rows','stable'));
end

EMM = cell(length(ind_full),3);
for i = 1:length(ind_full)
    EMM{i,1} = T{ind(i+1),1};
    EMM{i,2} = unique(Design(:,comb{i}),'rows','stable');
    EMM{i,3}(:,1) = (X_mean(:,ind_full{i})'*Mean_Big)./sum(X_mean(:,ind_full{i}),1)';
    EMM{i,3}(:,2) = sqrt(diag(((X_se(:,ind_full{i})'*E)/(size(E,1)/size(X_se(:,ind_full{i})',1)))*M*((X_se(:,ind_full{i})'*E)/(size(E,1)/size(X_se(:,ind_full{i})',1)))'));
    EMM{i,3}(:,3) = EMM{i,3}(:,1) + CI_fun(CI_alpha/2,df_res)*EMM{i,3}(:,2);
    EMM{i,3}(:,4) = EMM{i,3}(:,1) + CI_fun(1-CI_alpha/2,df_res)*EMM{i,3}(:,2);
end

end