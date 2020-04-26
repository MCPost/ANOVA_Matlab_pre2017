function [T, Info_Eff, X] = ANOVAperm(Y, Data, varargin)


%% List Variable Arguments Input
list_input = varargin;
list_input = reshape(list_input,2,[])';

if(sum(strcmp(list_input(:,1),'Sub')) == 1)
    if(isnumeric(list_input{strcmp(list_input,'Sub'),2}) && length(list_input{strcmp(list_input,'Sub'),2}) == 1)
        Sub = list_input{strcmp(list_input,'Sub'),2};
    else
        error('Error: Sub Index needs to be single integer!')
    end
else
    error('Error: Sub Index missing!');
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

if(~isempty(intersect(Sub,[BW WI])) || ~isempty(intersect(BW, WI)))
    error('Error: Sub, WI and BW have to be independent integers or vectors of integers!')
end

if(sum(strcmp(list_input(:,1), 'chan_label')) == 1)
    chan_label = list_input{strcmp(list_input(:,1), 'chan_label'),2};
    labels = {};
    for i = 1:size(chan_label,2)
        labels{i,1} = chan_label(i).label;
    end
else
    error('Error: Need chan_label');
end

if(sum(strcmp(list_input(:,1), 'type')) == 1)
    type = list_input{strcmp(list_input(:,1), 'type'),2};
else
    type = 3;
end

if(sum(strcmp(list_input(:,1), 'n_permutes')) == 1)
    n_permutes = list_input{strcmp(list_input(:,1), 'n_permutes'),2};
else
    n_permutes = 1000;
end

if(sum(strcmp(list_input(:,1), 'voxel_pval')) == 1)
    voxel_pval = list_input{strcmp(list_input(:,1), 'voxel_pval'),2};
else
    voxel_pval = 0.05;
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


%% Get Info
[X, Info_Eff, S_ind] = get_X_y_Info_S(Data, Sub, BW, WI, method);
N = size(S_ind,2);

sz_int = size(Y);
Err_list = [1; find(~cellfun(@isempty, regexp(Info_Eff(:,1),'Err')))];
if(~isempty(WI))
    
    for pt = 1:prod(sz_int(2:end))
        A = bsxfun(@times, Info_Eff{Err_list(end-1)+1,3}, squeeze(Y(:,pt)) - X*((X'*X)\X')*squeeze(Y(:,pt)));
        a = [];
        for j = 1:size(A,2)
            a(:,j) = A(Info_Eff{Err_list(end-1)+1,3}(:,j) ~= 0,j);
        end
        a = a - ones(size(a,1),size(a,1))'*a.*(1/size(a,1));
        Z(:,:,pt) = a'*a;
    end
    
    All_WI = unique(X(:,cell2mat(Info_Eff(Err_list(2:end-1)+1,2)')),'rows','stable');

end
clear a pt


%% Permutation Procedure
T = Info_Eff(:,1);
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

Hyp_Mat = {};

h = waitbar(0,sprintf('0 of %d Effects finished!',size(Info_Eff,1)),'name','Real F Values');
for i = 1:size(Info_Eff,1)
        
    if(sum(i == Err_list) == 0)
        
        % Real Data
        % Type of SS
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
        
        % Reduced Residuals
        HX1 = X1*((X1'*X1)\X1');
        for pt = 1:prod(sz_int(2:end))
            Y_rr(:,pt) = (eye(size(HX1,1)) - HX1)*squeeze(Y(:,pt));
        end
        Y_rr = reshape(Y_rr, [size(Y_rr,1) sz_int(2:end)]);
        
        Hyp_Mat{i,1} = Y_rr;
        
        % Effect SS
        X2rr = (eye(size(HX1,1)) - HX1)*X2;
        HX2rr = X2rr*((X2rr'*X2rr)\X2rr');
        
        for pt = 1:prod(sz_int(2:end))
            T{i,2}(pt) = squeeze(Y_rr(:,pt))'*HX2rr*squeeze(Y_rr(:,pt));
        end
        T{i,2} = reshape(T{i,2},sz_int(2:end));
        T{i,3} = size(X2rr,2);
        T{i,4} = T{i,2} ./ T{i,3};
        
        Hyp_Mat{i,2} = HX2rr;
        
        % Error SS  
        if(isempty(T{Err_list(find(Err_list > i,1,'first')),4})) %Err_list(find(Err_list > i,1,'first')-1)+1 == i
            Z2rr = orth((eye(size(HX1,1)) - HX1)*Info_Eff{Err_list(find(Err_list > i,1,'first')),3});
            HZ2rr = Z2rr*((Z2rr'*Z2rr)\Z2rr');

            for pt = 1:prod(sz_int(2:end))
                T{Err_list(find(Err_list > i,1,'first')),2}(pt) = squeeze(Y_rr(:,pt))'*(HZ2rr)*squeeze(Y_rr(:,pt));
            end
            T{Err_list(find(Err_list > i,1,'first')),2} = reshape(T{Err_list(find(Err_list > i,1,'first')),2},sz_int(2:end));
            T{Err_list(find(Err_list > i,1,'first')),3} = 32;%size(Z2rr,2);
            T{Err_list(find(Err_list > i,1,'first')),4} = T{Err_list(find(Err_list > i,1,'first')),2} ./ T{Err_list(find(Err_list > i,1,'first')),3};
        end
        
        Hyp_Mat{i,3} = HZ2rr;
        
        T{i,5} = T{i,4} ./ T{Err_list(find(Err_list > i,1,'first')),4};
        
        % Sphericity Correction
%         if(~isempty(Info_Eff{i,3}) && Err_list(2) < i)
%             
%             M = orth(All_WI(:,1:length(Info_Eff{i,2})))';
%             All_WI(:,1:length(Info_Eff{i,2})) = [];
%             
%             e_GG = zeros(1,prod(sz_int(2:end))); e_HF = zeros(1,prod(sz_int(2:end)));
%             for pt = 1:prod(sz_int(2:end))
%                 S = M*Z(:,:,pt)*M';
%                 k = size(S,1);
%                 
%                 e_GG(1,pt) = sum(diag(S))^2 / (k*sum(sum(S.^2)));
%                 e_HF(1,pt) = min(((N)*k*e_GG(1,pt) - 2) / (k*(N - T{Err_list(2)-1,3} - k*e_GG(1,pt))), 1);
%             end
%             
%             T(i:Err_list(Err_list > i)-1,7) = repmat({reshape(e_GG,[],sz_int(2:end))}, length(i:Err_list(Err_list > i)-1),1);
%             T(i:Err_list(Err_list > i)-1,9) = repmat({reshape(e_HF,[],sz_int(2:end))}, length(i:Err_list(Err_list > i)-1),1);
%             
%         end
    
    end
    waitbar(i/size(Info_Eff,1),h,sprintf('%d of %d Effects finished!',i,size(Info_Eff,1)),'name','Real F Values')
end
close(h)


%% Surrogate Data
df_eff = 1;
df_err = 32;

names_eff = transpose(2:size(T,1));
names_eff(~cellfun(@isempty, regexp(T(2:end,1),'_Err'))) = [];

h = waitbar(0,sprintf('0 of %d Permutation finished!',n_permutes),'name','Surrogate F-Values'); time = [];
for permi = 1:n_permutes
    tic
    for pt = 1:prod(sz_int(2:end))
        
        for eff = 1:length(names_eff)
            rand_vec = randperm(size(Hyp_Mat{names_eff(eff),1},1));
            T{names_eff(eff),6}(permi,pt) = (df_err.*(squeeze(Hyp_Mat{names_eff(eff),1}(rand_vec,pt))'*Hyp_Mat{names_eff(eff),2}*squeeze(Hyp_Mat{names_eff(eff),1}(rand_vec,pt))))./(df_eff.*(squeeze(Hyp_Mat{names_eff(eff),1}(rand_vec,pt))'*Hyp_Mat{names_eff(eff),3}*squeeze(Hyp_Mat{names_eff(eff),1}(rand_vec,pt))));
        end
        
    end
    time(permi) = toc;
    sec = median(time)*n_permutes - sum(time);
    waitbar(permi/n_permutes,h,sprintf('%d of %d Permutation finished!  (%dh %dmin)',permi,n_permutes,floor(sec/3600),round(60*(sec/3600 - floor(sec/3600)))),'name','Surrogate F-Values')
end
close(h)

prct = zeros(length(names_eff),length(labels));
for eff = 1:length(names_eff)
    T{names_eff(eff),6} = reshape(T{names_eff(eff),6}, [size(T{names_eff(eff),6},1) sz_int(2:end)]);
    help = permute(T{names_eff(eff),6},[2 1 3]);
    prct(eff,:) = prctile(help(:,:)',(1-voxel_pval)*100);
    %prct(eff,:) = finv((1 - voxel_pval)*ones(1,length(labels)), df_eff, df_err);
end

h = waitbar(0,'0 % of Cluster Correction finished!','name','Surrogate F-Values');
clust_pixel_info = [];
for permi = 1:70%n_permutes
    for eff = 1:length(names_eff)
        for chan = 1:length(labels)
            F_surrog = squeeze(T{names_eff(eff),6}(permi,chan,:));
            T{names_eff(eff),7}(permi,chan,1) = max(F_surrog(:));
            F_surrog(F_surrog < prct(eff,chan)) = 0;
            clustinfo = bwconncomp(F_surrog);
            T{names_eff(eff),7}(permi,chan,2) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]);
            T{names_eff(eff),7}(permi,chan,3) = max([ 0 cell2mat(cellfun(@(x) sum(F_surrog(x)),clustinfo.PixelIdxList,'UniformOutput',0))]);
        end
    end
    waitbar(permi/n_permutes,h,sprintf('%3.1f %% of Cluster Correction finished!',permi/n_permutes*100),'name','Surrogate F-Values')
end
close(h)

for eff = 1:length(names_eff)
    
    help = {};
    for chan = 1:size(chan_label,2)

        pdata = T{names_eff(eff),5}(chan,:);
        pdata(pdata < prct(eff,chan)) = 0;
        clustinfo = bwconncomp(pdata);
        help{chan,1} = clustinfo.PixelIdxList(cell2mat(cellfun(@(x) sum(pdata(x)),clustinfo.PixelIdxList,'UniformOutput',0)) > prctile(T{names_eff(eff),7}(chan,:,3),(1 - voxel_pval)*100));

    end

    help2 = help;
    cl = 2;
    for chan = 1:size(help,1)
        for clust1 = 1:size(help{chan},2)
            if(~isempty(help{chan,1}{clust1}))
                help2{chan,cl} = help{chan}{clust1}';
                help_inters = help2{chan,cl};
                for neighb = 1:size(help,1)
                    if(neighb ~= chan && ~isempty(help{neighb}))
                        for clust2 = 1:size(help{neighb},2)
                            if(length(intersect(help_inters,help{neighb}{clust2}')) > 2)
                                help_inters = intersect(help_inters, help{neighb}{clust2}');
                                help2{neighb,cl} = help{neighb}{clust2}';
                                help{neighb}{clust2} = [];
                            end
                        end
                    end
                end

                cl = cl + 1;
            end
        end
    end
    
    ind = [0,zeros(1,size(help2,2)-1)];
    for col = 2:size(help2,2)
        lab_ind = find(~cellfun(@isempty, help2(:,col)));
        for row = 1:length(lab_ind)
            if(isempty(intersect(labels(lab_ind)', chan_label(lab_ind(row)).neighblabel)))
                help2{lab_ind(row),col} = [];
            end
        end
        if(sum(~cellfun(@isempty, help2(:,col))) > 4)
            ind(1,col) = 0;
        else
            ind(1,col) = 1;
        end
    end

    help2(:,logical(ind)) = [];
    help2(:,1) = labels;
    
    T{names_eff(eff),8} = help2;
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, Info_Eff, S_ind] = get_X_y_Info_S(Data, Sub, BW, WI, method)

X = ones(size(Data,1),1);
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