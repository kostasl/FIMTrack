



%% Import Data Sets
strDir = '/home/klagogia/Videos/FIMTrackVideos/OregonR/';
[WT_Dat,N,srcDir] = readDataForAnalysis(strDir);
strDir = '/home/klagogia/Videos/FIMTrackVideos/v1016 AbetaArc';
[ABv1016_Dat,N,srcDir] = readDataForAnalysis(strDir);

strDir = '/home/klagogia/Videos/FIMTrackVideos/v1190 APPTau';
[v1190_Dat,N,srcDir] = readDataForAnalysis(strDir);


strDir = '/home/klagogia/Videos/FIMTrackVideos/attP40xC155 (gen ctr BWD)';
[attP40_Dat,N,srcDir] = readDataForAnalysis(strDir);




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

%Translate Data column colNum to Data Source File
%colNum = 1;
%WT_list_velocities(colNum,:)
%velo_wildtype.Properties.VarNames{colNum+1} 

export(velo_wildtype,'File','WT_velocities.csv','WriteVarNames',true);

dlmwrite('WT_velocities.csv',[0 nanmean(WT_list_velocities)],'delimiter','\t','-append');
dlmwrite('WT_velocities.csv',[0 nanmedian(WT_list_velocities)],'delimiter','\t','-append');


%ABv1016 
velo_ABv1016 = selectFeature(ABv1016_Dat, 'velo');

ABv1016_list_velocities = double(velo_ABv1016(:,2:end));
% list_velocities = list_velocities(tbl_velocities > speed_thres);
mean_velo = nanmean(ABv1016_list_velocities(ABv1016_list_velocities>speed_thres));
median_velo = nanmedian(ABv1016_list_velocities(ABv1016_list_velocities>speed_thres));

ABv1016_OverallMeanSpeed = nanmean(mean_velo)
ABv1016_OverallMedianSpeed = nanmedian(median_velo)
%mean_velo = mean_velo(~isnan(mean_velo))

export(velo_ABv1016,'File','ABv1016_velocities.csv','WriteVarNames',true);

dlmwrite('ABv1016_velocities.csv',[0 nanmean(ABv1016_list_velocities)],'delimiter','\t','-append');
dlmwrite('ABv1016_velocities.csv',[0 nanmedian(ABv1016_list_velocities)],'delimiter','\t','-append');


%% v1190
velo_v1190 = selectFeature(v1190_Dat, 'velo');

v1190_list_velocities = double(velo_v1190(:,2:end));
% list_velocities = list_velocities(tbl_velocities > speed_thres);
mean_velo = nanmean(v1190_list_velocities(v1190_list_velocities>speed_thres));
median_velo = nanmedian(v1190_list_velocities(v1190_list_velocities>speed_thres));

v1190_OverallMeanSpeed = nanmean(mean_velo)
v1190_OverallMedianSpeed = nanmedian(median_velo)
%mean_velo = mean_velo(~isnan(mean_velo))

export(velo_v1190,'File','v1190_velocities.csv','WriteVarNames',true);

dlmwrite('v1190_velocities.csv',[0 nanmean(v1190_list_velocities)],'delimiter','\t','-append');
dlmwrite('v1190_velocities.csv',[0 nanmedian(v1190_list_velocities)],'delimiter','\t','-append');


%% attP40
velo_attP40 = selectFeature(attP40_Dat, 'velo');

attP40_list_velocities = double(velo_attP40(:,2:end));
% list_velocities = list_velocities(tbl_velocities > speed_thres);
mean_velo = nanmean(attP40_list_velocities(attP40_list_velocities>speed_thres));
median_velo = nanmedian(attP40_list_velocities(attP40_list_velocities>speed_thres));

attP40_OverallMeanSpeed = nanmean(mean_velo)
attP40_OverallMedianSpeed = nanmedian(median_velo)
%mean_velo = mean_velo(~isnan(mean_velo))
velo_attP40_exp = join(velo_attP40,mean_velo,median_velo)
export(velo_attP40,'File','attP40_velocities.csv','WriteVarNames',true);

dlmwrite('attP40_velocities.csv',[0 nanmean(attP40_list_velocities)],'delimiter','\t','-append');
dlmwrite('attP40_velocities.csv',[0 nanmedian(attP40_list_velocities)],'delimiter','\t','-append');

%% Produce Plots
%boxplot(WT_list_velocities);
%hist(WT_list_velocities,20);



%figure;
%boxplot(ABv1016_list_velocities);




%%Save
save('ImportedData.mat');