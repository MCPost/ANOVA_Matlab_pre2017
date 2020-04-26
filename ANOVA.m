function [T, Info_Eff, X] = ANOVA(Data, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input Arguments:
% 
%' Data  -  Nxp matrix. N number of all observation. p contains one column
%           for Subjects (Sub), one column for the dependent variable (DV), 
%           and at least one column for a between-subjects (BW) or
%           within-subjects (WI) factor. Matrix has to be in long format,
%           so each row represents one distinct observation of the DV
%' Sub   -  integer indicating the column of Data containing the Subjects    
%' DV    -  integer indicating the column of Data containing the DV
%' BW    -  integer or vector of integers indicating the column of Data
%           containing all the between-subjects effects
%' WI    -  integer or vector of integers indicating the column of Data
%           containing all the within-subjects effects
%' type  -  integer indicating Type of Sum of Squares. Options are 1 for
%           Type I, 2 for Type II and 3 for Type III (default is Type III)
%' BW_names - cell vector of strings containing the names of the different
%             between-subject factors. Must be the same length and exact
%             succession of BW (default will be integers)
%' WI_names - cell vector of strings containing the names of the different
%             within-subject factors. Must be the same length and exact
%             succession of WI (default will be integers
%
% Output Arguments:
%
%' T        - table containing the results
%' Info_Eff - cell containing the design and hypotheses matrix for each
%             effect in the ANOVA
%' X        - the design matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


list_input = varargin;
list_input = reshape(list_input,2,[])';

if(sum(strcmp(list_input(:,1),'Sub')) == 1 && sum(strcmp(list_input(:,1),'AV')) == 1)
    if(isnumeric(list_input{strcmp(list_input,'Sub'),2}) && isnumeric(list_input{strcmp(list_input,'AV'),2})... 
        && length(list_input{strcmp(list_input,'Sub'),2}) == 1 && length(list_input{strcmp(list_input,'AV'),2}) == 1)
        Sub = list_input{strcmp(list_input,'Sub'),2};
        AV = list_input{strcmp(list_input,'AV'),2};
    else
        error('Error: Sub Index and AV Index need to be single integer!')
    end
else
    error('Error: Sub Index and/or AV Index missing!');
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

if(Sub == AV || ~isempty(intersect([Sub AV],[BW WI])) || ~isempty(intersect(BW, WI)))
    error('Error: Sub, AV, WI and BW have to be independent integers or vectors of integers!')
end

if(sum(strcmp(list_input(:,1), 'type')) == 1)
    type = list_input{strcmp(list_input(:,1), 'type'),2};
else
    type = 3;
end

if(sum(strcmp(list_input(:,1), 'X_coding')) == 1)
    switch list_input{strcmp(list_input(:,1), 'X_coding'),2}
        case 'SR'
            method = 'SR';
        case 'Treatment'
            method = 'Treatment';
        otherwise
            method = 'SR';
            lastwarn('Wrong X_coding argument. Use "SR" for sigma-restricted or "Treatment" for treatment coding! "SR" now used by default!')
    end
else
    method = 'SR';
end


Min_Design = unique(Data(:,[Sub BW  WI]),'rows');
if(size(Min_Design,1) < size(Data,1))
    Meas = accumarray(Data(:,[Sub fliplr(WI)]),Data(:,AV),[],@mean);
    for i = 1:size(Meas,1)
        Min_Design(Min_Design(:,1) == i,length([Sub BW  WI])+1) = Meas(i,:)';
    end
else
    Min_Design(:,end+1) = Data(:,AV);
end
BW_min = 2:1+length(BW);
WI_min = 1+length(BW)+(1:length(WI));
AV_min = size(Min_Design,2);


[X, y, Info_Eff, S_ind] = get_X_y_Info_S(Min_Design, 1, AV_min, BW_min, WI_min, method);
N = size(S_ind,2);

Err_list = [1; find(~cellfun(@isempty, regexp(Info_Eff(:,1),'Err')))];
if(isempty(WI))
    VariableNames = {'Effect','SS','df','MS','F','p'};
else
    VariableNames = {'Effect','SS','df','MS','F','p','e_GG','p_GG','e_HF','p_HF'};
    
    A = bsxfun(@times, Info_Eff{Err_list(end-1)+1,3}, (y - X*((X'*X)\X')*y));
    a = [];
    for j = 1:size(A,2)
        a(:,j) = A(Info_Eff{Err_list(end-1)+1,3}(:,j) ~= 0,j);
    end
    a = a - ones(size(a,1),size(a,1))'*a.*(1/size(a,1));
    Z = a'*a;
    
    All_WI = unique(X(:,cell2mat(Info_Eff(Err_list(2:end-1)+1,2)')),'rows','stable');

end


switch type
    case 1
        VariableNames{2} = 'SS_Type_I';
    case 2
        VariableNames{2} = 'SS_Type_II';
    case 3
        VariableNames{2} = 'SS_Type_III';
end

if(sum(strcmp(list_input(:,1), 'effectsize')) == 1)
    switch list_input{strcmp(list_input(:,1), 'effectsize'),2}
        case 'eta'
            VariableNames{end+1} = 'part_eta2';
            ES = 'eta';
        case 'omega'
            VariableNames{end+1} = 'part_omega2';
            ES = 'omega';
        otherwise
            lastwarn('Wrong effectsize argument. Use "omega" for partial omega^2 or "eta" for partial eta^2')
    end
end


T = Info_Eff(:,1);
for i = 1:size(Info_Eff)
    
    if(i == 1)
        
        T{i,2} = y'*X(:,1)*((X(:,1)'*X(:,1))\X(:,1)')*y;
        T{i,3} = 1;
        T{i,4} = T{i,2};
        Int = 0;
        
    elseif(sum(i == Err_list) == 0)
        
        switch type
            case 1
                X1 = X(:,[1 cell2mat(Info_Eff(Err_list(find(i > Err_list,1,'first'))+1:i-1,2)')]);
            case 2
                non_cont_list = ~(cellfun(@length, regexp(Info_Eff(:,1), strrep(Info_Eff{i,1},'_','|'))) == length(regexp(Info_Eff{i,1},'_'))+1);
                X1 = X(:,cell2mat(Info_Eff(non_cont_list,2)'));
            case 3
                X1 = X(:,[1 cell2mat(Info_Eff([Err_list(find(i > Err_list,1,'first'))+1:i-1 i+1:Err_list(find(i < Err_list,1,'first'))-1],2)')]);
        end
        X2 = X(:,Info_Eff{i,2});
        
        HX1 = X1*((X1'*X1)\X1');
        X2rr = (eye(size(HX1,1)) - HX1)*X2;
        HX2rr = X2rr*((X2rr'*X2rr)\X2rr');
        
        T{i,2} = y'*HX2rr*y;
        T{i,3} = size(X2rr,2);
        T{i,4} = T{i,2} / T{i,3};
        
        if(~isempty(Info_Eff{i,3}) && Err_list(2) < i)
            
            M = orth(All_WI(:,1:length(Info_Eff{i,2})))';
            All_WI(:,1:length(Info_Eff{i,2})) = [];
            
            S = M*Z*M';
            k = size(S,1);
            
            T(i:Err_list(Err_list > i)-1,7) = num2cell(repmat( sum(diag(S))^2 / (k*sum(sum(S.^2))), length(i:Err_list(Err_list > i)-1),1));
            %T(i:Err_list(Err_list > i)-1,9) = num2cell(repmat( min(((N - T{Err_list(2)-1,3} - 1)*k*T{i,7} - 2) / (k*(N - T{Err_list(2)-1,3} - k*T{i,7})), 1), length(i:Err_list(Err_list > i)-1),1));
            T(i:Err_list(Err_list > i)-1,9) = num2cell(repmat( min(((N)*k*T{i,7} - 2) / (k*(N - T{Err_list(2)-1,3} - k*T{i,7})), 1), length(i:Err_list(Err_list > i)-1),1));
            
        end
        
    else
        
        X3 = X(:,[1 cell2mat(Info_Eff(Err_list(find(Err_list == i)-1)+1:i-1,2)')]);
        
        HX3 = X3*((X3'*X3)\X3');
        Z2rr = orth((eye(size(HX3,1)) - HX3)*Info_Eff{i,3});
        HZ2rr = Z2rr*((Z2rr'*Z2rr)\Z2rr');
        
        T{i,2} = y'*(HZ2rr)*y;
        T{i,3} = size(Z2rr,2);
        T{i,4} = T{i,2} / T{i,3};
        
        T(Err_list(find(Err_list == i)-1)+Int:i-1,5) = num2cell(cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,4)) ./ T{i,4});
        T(Err_list(find(Err_list == i)-1)+Int:i-1,6) = num2cell(fcdf(cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,5)), cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,3)),T{i,3},'upper'));
        if(~strcmp(Info_Eff{i,1},'S_Err'))
            T(Err_list(find(Err_list == i)-1)+Int:i-1,8) = num2cell(fcdf(cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,5)), cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,3)) .* cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,7)), T{i,3} .* cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,7)), 'upper'));
            T(Err_list(find(Err_list == i)-1)+Int:i-1,10) = num2cell(fcdf(cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,5)), cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,3)) .* cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,9)), T{i,3} .* cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,9)), 'upper'));
        end
        
        if(exist('ES','var'))
            if(strcmp(ES, 'eta'))
                T(Err_list(find(Err_list == i)-1)+Int:i-1,length(VariableNames)) = num2cell(cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,2)) ./ (cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,2)) + T{i,2}));
            else
                T(Err_list(find(Err_list == i)-1)+Int:i-1,length(VariableNames)) = num2cell( (cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,3)) .* (cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,4)) - T{i,4})) ./ (cell2mat(T(Err_list(find(Err_list == i)-1)+Int:i-1,2)) + sum(cell2mat(T(strcmp(T(:,1),'Err'),2))) + cell2mat(T(strcmp(T(:,1),'S_Err'),4))) );
            end
        end
        
        Int = 1;
    end
    
end

names_ind_bw = [];
if(exist('BW_names','var'))
    if(length(BW_names) == size(BW,2))
        names_ind_bw = 1:length(BW_names);
        for i = names_ind_bw
            T(:,1) = regexprep(T(:,1), num2str(names_ind_bw(i)), BW_names(i));
        end
    else
        error('Error: BW_names has %i entries, but BW %i factors!',length(BW_names),size(BW,2))
    end
end
if(exist('WI_names','var'))
    if(length(WI_names) == size(WI,2))
        names_ind_wi = (1:length(WI_names)) + length(names_ind_bw);
        for i = names_ind_wi-length(names_ind_bw)
            T(:,1) = regexprep(T(:,1), num2str(names_ind_wi(i)), WI_names(i));
        end
    else
        error('Error: WI_names has %i entries, but WI %i factors!',length(WI_names),size(WI,2))
    end
end


T = cell2table(T,'VariableNames',VariableNames);

end


function [X, y, Info_Eff, S_ind] = get_X_y_Info_S(Data, Sub, AV, BW, WI, method)

X = ones(size(Data,1),1);
y = Data(:,AV);
S = Data(:,Sub);

S_ind = dummyvar(S);
S_ind(:,sum(S_ind,1) == 0) = [];
S_ind2 = S_ind;
if(exist('method','var'))
    switch method
        case 'SR'
            S_ind2 = bsxfun(@minus, S_ind(:,1:end-1), S_ind(:,end));
        case 'Treatment'
            S_ind2(:,end) = [];
    end
end

Info_Eff = {'Intercept',1};
end_ineff = 1;
if(exist('BW','var'))
        
    bw_Eff = {};
    for i = 1:length(BW)
       bw_Eff = [bw_Eff; strrep(cellstr(num2str(nchoosek(1:length(BW),i))),'  ','_')];
    end

    for i = 1:size(bw_Eff,1)
        if((length(regexp(bw_Eff{i,1},'_'))+1) == 1)
            H = dummyvar(Data(:,BW(i)));
            Info_Eff{i+end_ineff,3} = H;
            if(exist('method','var'))
                switch method
                    case 'SR'
                        H = bsxfun(@minus, H(:,1:end-1), H(:,end));
                    case 'Treatment'
                        H(:,end) = [];
                end
            end
            Info_Eff{i+end_ineff,1} = bw_Eff{i,1};
            Info_Eff{i+end_ineff,2} = (size(X,2)+1):((size(X,2)) + size(H,2));
            X = [X H];
        else                                                                
            [sp, fp] = strtok(fliplr(bw_Eff{i,1}),'_'); fp = fp(2:end);
            fp = fliplr(fp);
            x_h1 = X(:,Info_Eff{strcmp(Info_Eff(:,1), fp),2});
            x_h2 = X(:,Info_Eff{strcmp(Info_Eff(:,1), sp),2});
            IA = zeros(size(x_h1,1), size(x_h1,2)*size(x_h2,2));
            for m = 1:size(x_h1,2)
                for n = 1:size(x_h2,2)
                    IA(:,(m-1)*size(x_h2,2)+n) = x_h1(:,m).*x_h2(:,n);
                end
            end
            Info_Eff{i+end_ineff,1} = bw_Eff{i,1};
            Info_Eff{i+end_ineff,2} = (size(X,2)+1):((size(X,2)) + size(IA,2));
            X = [X IA];
        end
    end
    
    Info_Eff{end+1,1} = 'S_Err';
    Info_Eff{end,3} = S_ind2;
    
else
    
    Info_Eff{1+end_ineff,1} = 'S_Err';
    Info_Eff{1+end_ineff,3} = S_ind2;
    
end

end_ineff = length(Info_Eff(:,1));

if(exist('WI','var'))
    
    wi_Eff = {};
    for i = 1:length(WI)
       wi_Eff = [wi_Eff; strrep(cellstr(num2str(nchoosek((1:length(WI))+sum(cellfun(@length, Info_Eff(:,1)) == 1),i))),'  ','_')];
    end
    
    for i = 1:length(wi_Eff)
     
        if((length(regexp(wi_Eff{i,1},'_'))+1) == 1)
            H = dummyvar(Data(:,WI(i)));
            Info_Eff{end_ineff+1,3} = H;
            if(exist('method','var'))
                switch method
                    case 'SR'
                        H = bsxfun(@minus, H(:,1:end-1), H(:,end));
                    case 'Treatment'
                        H(:,end) = [];
                end
            end
            Info_Eff{end_ineff+1,1} = wi_Eff{i,1};
            Info_Eff{end_ineff+1,2} = (size(X,2)+1):((size(X,2)) + size(H,2));
            X = [X H];
        else                                                                
            [sp, fp] = strtok(fliplr(wi_Eff{i,1}),'_'); fp = fp(2:end);
            fp = fliplr(fp); sp = fliplr(sp);
            x_h1 = X(:,Info_Eff{strcmp(Info_Eff(:,1), fp),2});
            x_h2 = X(:,Info_Eff{strcmp(Info_Eff(:,1), sp),2});
            IA = zeros(size(x_h1,1), size(x_h1,2)*size(x_h2,2));
            for m = 1:size(x_h1,2)
                for n = 1:size(x_h2,2)
                    IA(:,(m-1)*size(x_h2,2)+n) = x_h1(:,m).*x_h2(:,n);
                end
            end
            Info_Eff{end_ineff+1,1} = wi_Eff{i,1};
            Info_Eff{end_ineff+1,2} = (size(X,2)+1):((size(X,2)) + size(IA,2));
            X = [X IA];
            
            x_h1 = Info_Eff{strcmp(Info_Eff(:,1), fp),3};
            x_h2 = Info_Eff{strcmp(Info_Eff(:,1), sp),3};
            IA = zeros(size(x_h1,1), size(x_h1,2)*size(x_h2,2));
            for m = 1:size(x_h1,2)
                for n = 1:size(x_h2,2)
                    IA(:,(m-1)*size(x_h2,2)+n) = x_h1(:,m).*x_h2(:,n);
                end
            end
            Info_Eff{end_ineff+1,3} = IA;
        end
        
        if(exist('BW','var'))
            
            for j = 1:size(bw_Eff,1)
                x_h1 = X(:,Info_Eff{strcmp(Info_Eff(:,1), bw_Eff{j,1}),2});
                x_h2 = X(:,Info_Eff{strcmp(Info_Eff(:,1), wi_Eff{i,1}),2});
                IA = zeros(size(x_h1,1), size(x_h1,2)*size(x_h2,2));
                for m = 1:size(x_h1,2)
                    for n = 1:size(x_h2,2)
                        IA(:,(m-1)*size(x_h2,2)+n) = x_h1(:,m).*x_h2(:,n);
                    end
                end
                Info_Eff{1+end_ineff+j,1} = [bw_Eff{j,1},'_',wi_Eff{i,1}];
                Info_Eff{1+end_ineff+j,2} = (size(X,2)+1):((size(X,2)) + size(IA,2));
                X = [X IA];
            end
            
        end
        
        x_h1 = S_ind2;
        x_h2 = X(:,Info_Eff{strcmp(Info_Eff(:,1), wi_Eff{i,1}),2});
        IA = zeros(size(x_h1,1), size(x_h1,2)*size(x_h2,2));
        for m = 1:size(x_h1,2)
            for n = 1:size(x_h2,2)
                IA(:,(m-1)*size(x_h2,2)+n) = x_h1(:,m).*x_h2(:,n);
            end
        end
        Info_Eff{end+1,1} = ['S_',wi_Eff{i,1},'_Err'];
        Info_Eff{end,3} = IA;
        
        end_ineff = size(Info_Eff,1);
        
    end
    
    
end

end