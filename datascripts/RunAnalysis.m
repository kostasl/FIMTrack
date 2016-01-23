



%% Import Data Sets
strDir = '/home/klagogia/Videos/FIMTrackVideos/OregonR/';
[WT_Dat,N,srcDir] = readDataForAnalysis(strDir);
strDir = '/home/klagogia/Videos/FIMTrackVideos/v1016 AbetaArc';
[ABv1016_Dat,N,srcDir] = readDataForAnalysis(strDir);

%% EXTRACT DATA and ANALYSE for Mean Speeds
speed_thres = 0;
%Analyse WT
velo_wildtype = selectFeature(WT_Dat, 'velo');

WT_list_velocities = double(velo_wildtype(:,2:end));
% list_velocities = list_velocities(tbl_velocities > speed_thres);
mean_velo = nanmean(WT_list_velocities(WT_list_velocities>speed_thres));
median_velo = nanmedian(WT_list_velocities(WT_list_velocities>speed_thres));

WT_OverallMeanSpeed  = nanmean(mean_velo)
WT_OverallMedianSpeed = nanmedian(median_velo)
%mean_velo = mean_velo(~isnan(mean_velo))

%ABv1016 
velo_ABv1016 = selectFeature(ABv1016_Dat, 'velo');

ABv1016_list_velocities = double(velo_ABv1016(:,2:end));
% list_velocities = list_velocities(tbl_velocities > speed_thres);
mean_velo = nanmean(ABv1016_list_velocities(ABv1016_list_velocities>speed_thres));
median_velo = nanmedian(ABv1016_list_velocities(ABv1016_list_velocities>speed_thres));

ABv1016_OverallMeanSpeed = nanmean(mean_velo)
ABv1016_OverallMedianSpeed = nanmedian(median_velo)
%mean_velo = mean_velo(~isnan(mean_velo))

%% Produce Plots
boxplot(WT_list_velocities);
hist(WT_list_velocities,20);


figure;
boxplot(ABv1016_list_velocities);




%%Save
save('ImportedData.mat');