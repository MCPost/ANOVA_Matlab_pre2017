%% Implement descriptives


Data = Data_BT_WI_3F_uneq;

[T, Info_Eff, X] = ANOVA(Data, 'Sub',1, 'AV',5, 'BW',2, 'WI',3:4, 'BW_names', {'Sex'}, 'WI_names',{'Shape','Angle'}, 'effectsize','eta');

a = accumarray(fliplr(Data(:,2:4)), Data(:,5), [], @(x) nanmean(x));
b = accumarray(fliplr(Data(:,2:4)), Data(:,5), [], @(x) nanstd(x,0,1));

a(:)
sqrt(T{12,4}/3)

