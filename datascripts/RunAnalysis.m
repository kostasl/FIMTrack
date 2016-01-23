
strDir = '/home/klagogia/Videos/FIMTrackVideos/OregonR/';

speed_thres = 1;

[WT_Dat,N,srcDir] = readDataForAnalysis(strDir);
velo_wildtype = selectFeature(WT_Dat, 'velo');

WT_list_velocities = double(velo_wildtype(:,2:end));
% list_velocities = list_velocities(tbl_velocities > speed_thres);
mean_velo = nanmean(WT_list_velocities(WT_list_velocities>speed_thres));

WT_OverallMeanSpeed = nanmean(mean_velo)
%mean_velo = mean_velo(~isnan(mean_velo))

boxplot(WT_list_velocities);
hist(WT_list_velocities,20);

strDir = '/home/klagogia/Videos/FIMTrackVideos/v1016 AbetaArc';
speed_thres = 1;
[ABv1016_Dat,N,srcDir] = readDataForAnalysis(strDir);
velo_ABv1016 = selectFeature(ABv1016_Dat, 'velo');

ABv1016_list_velocities = double(velo_ABv1016(:,2:end));
% list_velocities = list_velocities(tbl_velocities > speed_thres);
mean_velo = nanmean(ABv1016_list_velocities(ABv1016_list_velocities>speed_thres));

ABv1016_OverallMeanSpeed = nanmean(mean_velo)
%mean_velo = mean_velo(~isnan(mean_velo))
figure;
boxplot(ABv1016_list_velocities);


%%Save
save('ImportedData.mat');