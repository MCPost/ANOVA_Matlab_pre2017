%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% ANOVA implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Between-subjects Oneway ANOVA

% Equal Group Size
Data_BT_1F = [1:9; 1 1 1 2 2 2 3 3 3; 2 3 4 6 7 8 10 11 12]';

[T_BT_1F, Info_Eff, X] = ANOVA(Data_BT_1F, 'Sub',1, 'DV',3, 'BW',2, 'WI',[], 'BW_names',{'Group'}, 'effectsize','eta');
EMM_BT_1F = EstMargMeans(Data_BT_1F, T_BT_1F, Info_Eff, X, 'Sub',1, 'DV',3, 'BW',2);

% Unequal Group Size
keepout = [2 6];
Data_BT_1F_uneq = Data_BT_1F;
for i = 1:length(keepout)
    Data_BT_1F_uneq(Data_BT_1F_uneq(:,1) == keepout(i),:) = [];
end
clear keepout i

[T_BT_1F_uneq, Info_Eff, X] = ANOVA(Data_BT_1F_uneq, 'Sub',1, 'DV',3, 'BW',2, 'WI',[], 'BW_names',{'Group'}, 'effectsize','eta');
EMM_BT_1F_uneq = EstMargMeans(Data_BT_1F_uneq, T_BT_1F_uneq, Info_Eff, X, 'Sub',1, 'DV',3, 'BW',2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Between-subjects Two Factor 2x2 ANOVA

% Equal Group Size
Data_BT_2F = [1:20; kron([1 2],ones(1,10)); kron([1 2],ones(1,5)) kron([1 2],ones(1,5)); randi(10,1,20)]';

[T_BT_2F, Info_Eff, X] = ANOVA(Data_BT_2F, 'Sub',1, 'DV',4, 'BW',2:3, 'WI',[], 'BW_names',{'Sex','Group'}, 'effectsize','eta');
EMM_BT_2F = EstMargMeans(Data_BT_2F, T_BT_2F, Info_Eff, X, 'Sub',1, 'DV',4, 'BW',2:3);

% Unequal Group Size
keepout = [2 6 17];
Data_BT_2F_uneq = Data_BT_2F;
for i = 1:length(keepout)
    Data_BT_2F_uneq(Data_BT_2F_uneq(:,1) == keepout(i),:) = [];
end
clear keepout i

[T_BT_2F_uneq, Info_Eff, X] = ANOVA(Data_BT_2F_uneq, 'Sub',1, 'DV',4, 'BW',2:3, 'WI',[], 'BW_names',{'Sex','Group'}, 'effectsize','eta');
EMM_BT_2F_uneq = EstMargMeans(Data_BT_2F_uneq, T_BT_2F_uneq, Info_Eff, X, 'Sub',1, 'DV',4, 'BW',2:3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Between-subjects Three Factor 2x2x3 ANOVA

% Equal Group Size
Data_BT_3F = [1:60; kron([1 2],ones(1,30)); kron([1 2 1 2],ones(1,15)); kron([1 2 3 1 2 3 1 2 3 1 2 3],ones(1,5)); randi(10,1,60)]';

[T_BT_3F, Info_Eff, X] = ANOVA(Data_BT_3F, 'Sub',1, 'DV',5, 'BW',2:4, 'WI',[], 'BW_names',{'Sex','Group','Med'}, 'effectsize','eta');
EMM_BT_3F = EstMargMeans(Data_BT_3F, T_BT_3F, Info_Eff, X, 'Sub',1, 'DV',5, 'BW',2:4);

% Unequal Group Size
keepout = [2 6 9 17 28 33 42 57];
Data_BT_3F_uneq = Data_BT_3F;
for i = 1:length(keepout)
    Data_BT_3F_uneq(Data_BT_3F_uneq(:,1) == keepout(i),:) = [];
end
clear keepout i

[T_BT_3F_uneq, Info_Eff, X] = ANOVA(Data_BT_3F_uneq, 'Sub',1, 'DV',5, 'BW',2:4, 'WI',[],'BW_names',{'Sex','Group','Med'}, 'effectsize','eta');
EMM_BT_3F_uneq = EstMargMeans(Data_BT_3F_uneq, T_BT_3F_uneq, Info_Eff, X, 'Sub',1, 'DV',5, 'BW',2:4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Within-subjects Oneway ANOVA

% Equal Group Size
Data_WI_1F = [repmat(1:9,1,3); kron([1 2 3], ones(1,9)); randi(100,1,27)]';

[T_WI_1F, Info_Eff, X] = ANOVA(Data_WI_1F, 'Sub',1, 'DV',3, 'BW',[], 'WI',2, 'WI_names',{'Shape'}, 'effectsize','eta');
EMM_WI_1F = EstMargMeans(Data_WI_1F, T_WI_1F, Info_Eff, X, 'Sub',1, 'DV',3, 'WI',2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Within-subjects Two Factor 2x2 ANOVA

% Equal Group Size
Data_WI_2F = [repmat(1:3,1,9); kron([1 2 3], ones(1,9)); repmat(kron([1 2 3], ones(1,3)),1,3); randi(100,1,27)]';

[T_WI_2F, Info_Eff, X] = ANOVA(Data_WI_2F, 'Sub',1, 'DV',4, 'BW',[], 'WI',2:3, 'WI_names',{'Shape','Angle'}, 'effectsize','eta');
EMM_WI_2F = EstMargMeans(Data_WI_2F, T_WI_2F, Info_Eff, X, 'Sub',1, 'DV',4, 'WI',2:3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Within-subjects Three Factor 2x3x3 ANOVA

% Equal Group Size
Data_WI_3F = [kron([1 2 3 4], ones(1,18)); repmat(kron([1 2],ones(1,9)),1,4); repmat(kron([1 2 3], ones(1,3)),1,8); repmat([1 2 3],1,24); randi(100,1,72)]';

[T_WI_3F, Info_Eff, X] = ANOVA(Data_WI_3F, 'Sub',1, 'DV',5, 'BW',[], 'WI',2:4, 'WI_names',{'Treat','Shape','Angle'}, 'effectsize','eta');
EMM_WI_3F = EstMargMeans(Data_WI_3F, T_WI_3F, Info_Eff, X, 'Sub',1, 'DV',5, 'WI',2:4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mixed Two Factor 2x3 ANOVA

% Equal Group Size
Data_BT_WI_2F = [kron([1 2 3 4 5 6 7 8 9 10], ones(1,3)); kron([1 2],ones(1,15)); repmat([1 2 3],1,10); randi(100,1,30)]';

[T_BT_WI_2F, Info_Eff, X] = ANOVA(Data_BT_WI_2F, 'Sub',1, 'DV',4, 'BW',2, 'WI',3, 'BW_names', {'Sex'}, 'WI_names',{'Shape'}, 'effectsize','eta');
EMM_BT_WI_2F = EstMargMeans(Data_BT_WI_2F, T_BT_WI_2F, Info_Eff, X, 'Sub',1, 'DV',4, 'BW',2, 'WI',3);

% Unequal Group Size
keepout = [2 4 10];
Data_BT_WI_2F_uneq = Data_BT_WI_2F;
for i = 1:length(keepout)
    Data_BT_WI_2F_uneq(Data_BT_WI_2F_uneq(:,1) == keepout(i),:) = [];
end
clear keepout i

[T_BT_WI_2F_uneq, Info_Eff, X] = ANOVA(Data_BT_WI_2F_uneq, 'Sub',1, 'DV',4, 'BW',2, 'WI',3, 'BW_names', {'Sex'}, 'WI_names',{'Shape'}, 'effectsize','eta');
EMM_BT_WI_2F_uneq = EstMargMeans(Data_BT_WI_2F_uneq, T_BT_WI_2F_uneq, Info_Eff, X, 'Sub',1, 'DV',4, 'BW',2, 'WI',3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mixed Three Factor 2x3x3 ANOVA

% Equal Group Size
Data_BT_WI_3F = [kron([1 2 3 4 5 6 7 8 9 10], ones(1,9)); kron([1 2],ones(1,45)); repmat(kron([1 2 3], ones(1,3)),1,10); repmat([1 2 3],1,30); randi(100,1,90)]';

[T_BT_WI_3F, Info_Eff, X] = ANOVA(Data_BT_WI_3F, 'Sub',1, 'DV',5, 'BW',2, 'WI',3:4, 'BW_names', {'Sex'}, 'WI_names',{'Shape','Angle'}, 'effectsize','eta');
EMM_BT_WI_3F = EstMargMeans(Data_BT_WI_3F, T_BT_WI_3F, Info_Eff, X, 'Sub',1, 'DV',5, 'BW',2, 'WI',3:4);

% Unequal Group Size
keepout = [2 4 10];
Data_BT_WI_3F_uneq = Data_BT_WI_3F;
for i = 1:length(keepout)
    Data_BT_WI_3F_uneq(Data_BT_WI_3F_uneq(:,1) == keepout(i),:) = [];
end
clear keepout i

[T_BT_WI_3F_uneq, Info_Eff, X] = ANOVA(Data_BT_WI_3F_uneq, 'Sub',1, 'DV',5, 'BW',2, 'WI',3:4, 'BW_names', {'Sex'}, 'WI_names',{'Shape','Angle'}, 'effectsize','eta');
EMM_BT_WI_3F_uneq = EstMargMeans(Data_BT_WI_3F_uneq, T_BT_WI_3F_uneq, Info_Eff, X, 'Sub',1, 'DV',5, 'BW',2, 'WI',3:4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BIG Mixed Six Factor 2x3x3(BW) x 3x3x4(WI) ANOVA

% Equal Group Size
Data_BT_WI_6F = table2array(readtable('C:/Users/Christopher/Documents/SPSS/BigANOVA_long_1.csv','Delimiter',';','HeaderLines',0));

[T_BT_WI_6F, Info_Eff, X] = ANOVA(Data_BT_WI_6F, 'Sub',1, 'DV',8, 'BW',2:4, 'WI',5:7, 'BW_names', {'Sex','Group','Treat'}, 'WI_names',{'Month','Session','Timepoint'}, 'effectsize','eta');
EMM_BT_WI_6F = EstMargMeans(Data_BT_WI_6F, T_BT_WI_6F, Info_Eff, X, 'Sub',1, 'DV',8, 'BW',2:4, 'WI',5:7);

% Unequal Group Size
keepout = [10 17 24 29 45 60 63 78];
Data_BT_WI_6F_uneq = Data_BT_WI_6F;
for i = 1:length(keepout)
    Data_BT_WI_6F_uneq(Data_BT_WI_6F_uneq(:,1) == keepout(i),:) = [];
end

[T_BT_WI_6F_uneq, Info_Eff, X] = ANOVA(Data_BT_WI_6F_uneq, 'Sub',1, 'DV',8, 'BW',2:4, 'WI',5:7, 'BW_names', {'Sex','Group','Treat'}, 'WI_names',{'Month','Session','Timepoint'}, 'effectsize','eta');
EMM_BT_WI_6F_uneq = EstMargMeans(Data_BT_WI_6F_uneq, T_BT_WI_6F_uneq, Info_Eff, X, 'Sub',1, 'DV',8, 'BW',2:4, 'WI',5:7);

