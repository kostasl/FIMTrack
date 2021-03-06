

%Directory must contain only results directories of genotypes - no other
%files
%Home Comp.
baseDir = '/media/kostasl/FlashDrive/FIMTracker/trackerOutput/';

% Office
baseDir = '/home/public/FIMTrackVideos/Tracked/'; 


s=dir(baseDir);

%strDir = '/media/kostasl/FlashDrive/FIMTracker/FIMtrackerOutput/attP40xC155 (gen ctr BWD)';
%[attP40_Dat,N,srcDir] = readDataForAnalysis(strDir);


minTrackletLength = 200;
speed_thres = 3; %Set to <0 to remove filter

%% Produce Plots
%boxplot(WT_list_velocities);
%hist(WT_list_velocities,20);


%figure;
%boxplot(ABv1016_list_velocities);
%Ignore  1st 2 contain .  & ..
for k=3:size(s)

    
    if s(k).isdir == 0 
        continue; %skip if this is not a directory
    end
    
    strDir = strcat(baseDir,'/' ,s(k).name,'/')
    %
    %import data
    [gen_Dat,N,srcDir] = readDataForAnalysis(strDir);

    %Dataset and array - use array for calculations / 
    gen_velodataset = selectFeature(gen_Dat, 'velo');
    gen_list_velocities = double(gen_velodataset(:,2:end));

    %Remove Short Tracklets from Dataset
    vNanCountInColumns = sum(~isnan(gen_list_velocities),1);
    colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
    gen_velodataset(:,colIndexToRemove+1) = [];
    %Update The Matrix Of Tracklets
    gen_list_velocities = double(gen_velodataset(:,2:end));

    if speed_thres >= 0 
        %Remove Low Speeds
        cntExcluded = sum(gen_list_velocities <= speed_thres);
        %Need the Col Row indexing so I can Translate to the dataset
        %Coordinates by adding col+1
        [cellRow,cellCol] = find(gen_list_velocities <= speed_thres);
        %Only way I could find now is just Iterate through Indexes set to NaN
        for i=1:size(cellRow) gen_velodataset{cellRow(i),cellCol(i)+1} = NaN; end
    end

    %Update The Matrix Of Tracklets
    gen_list_velocities = double(gen_velodataset(:,2:end));
    runVelo = gen_list_velocities;
    dataset = gen_velodataset;

    mean_velo = nanmean(runVelo);
    median_velo = nanmedian(runVelo);
    count_nonNan = sum(~isnan(runVelo));
    %mean_velo = mean_velo(~isnan(mean_velo))

    filename = sprintf('dataout/velocity/%s-gt%d.csv',s(k).name,speed_thres);
    export(dataset,'File',filename,'WriteVarNames',true);

    
    dlmwrite(filename,[0 mean_velo],'delimiter','\t','-append');
    dlmwrite(filename,[0 median_velo],'delimiter','\t','-append');
    dlmwrite(filename,[0 count_nonNan],'delimiter','\t','-append');
    dlmwrite(filename,[0 cntExcluded],'delimiter','\t','-append');
    
    sprintf(filename);
    nanmean(mean_velo)
   
end



%%Save
save('ImportedData.mat');

%%%%%%%%%%%%%%%%%%%%5
%% bending
% NOT Complete

%figure;
%boxplot(ABv1016_list_velocities);
%Ignore  1st 2 contain .  & ..

%Remove nending angles size less or equal to absolute value %Defines cast
%threshold
angle_thres = 18; 
for k=3:size(s)
    
    if s(k).isdir == 0 
        continue; %skip if this is not a directory
    end
    
    strDir = strcat(baseDir,'/' ,s(k).name,'/')
    %
    %import data
    [gen_Dat,N,srcDir] = readDataForAnalysis(strDir);

    %Dataset and array - use array for calculations / 
    gen_benddataset = selectFeature(gen_Dat, 'bending');
    vars = get(gen_benddataset,'VarNames');
    for i=2:numel(vars)  %Update Data set column by column
        gen_benddataset.(vars{i}) = abs(gen_benddataset.(vars{i}) - 180);
    end
    
    gen_list_bends = double(gen_benddataset(:,2:end));

    %Remove Short Tracklets from Dataset
    vNanCountInColumns = sum(~isnan(gen_list_bends),1);
    colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
    gen_benddataset(:,colIndexToRemove+1) = [];
    %Update The Matrix Of Tracklets
    gen_list_bends = double(gen_benddataset(:,2:end));

    %Remove Low Angles
    cntExcluded = sum(gen_list_bends <= angle_thres);
    [cellRow,cellCol] = find(abs(gen_list_bends) <= angle_thres);
    
    %Only way I could find now is just Iterate through Indexes set to NaN
    for i=1:size(cellRow) gen_benddataset{cellRow(i),cellCol(i)+1} = NaN; end

    %Update The Matrix Of Tracklets
    gen_list_bends = double(gen_benddataset(:,2:end));
    runVelo = gen_list_bends;
    dataset = gen_benddataset;

    mean_velo   = nanmean(runVelo);
    median_velo = nanmedian(runVelo);
    count_bends =  sum(~isnan(runVelo));
    
    %mean_velo = mean_velo(~isnan(mean_velo))

    filename = sprintf('dataout/bending/bend-%s-gt%d.csv',s(k).name,angle_thres);
    export(dataset,'File',filename,'WriteVarNames',true);

    
    dlmwrite(filename,['M' mean_velo],'delimiter','\t','-append');
    dlmwrite(filename,['m' median_velo],'delimiter','\t','-append');
    dlmwrite(filename,['C' count_bends],'delimiter','\t','-append');
    dlmwrite(filename,['X' cntExcluded],'delimiter','\t','-append');
    

    sprintf(filename);
    nanmean(mean_velo)
   
end




%% EXTRACT DATA and ANALYSE for Mean Speeds
% 
% % Import Data Sets
% strDir = strcat(baseDir,'/OregonR/');
% [WT_Dat,N,srcDir] = readDataForAnalysis(strDir);
% 
% assert(N > 0,'No data loaded! Couldnt find file');
% %Analyse WT
% WT_velodataset = selectFeature(WT_Dat, 'velo');
% 
% WT_list_velocities = double(WT_velodataset(:,2:end));
% 
% %Remove Short Tracklets from Dataset
% vNanCountInColumns = sum(~isnan(WT_list_velocities),1);
% colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
% WT_velodataset(:,colIndexToRemove+1) = [];
% %Update The Matrix Of Tracklets
% WT_list_velocities = double(WT_velodataset(:,2:end));
% 
% % list_velocities = list_velocities(tbl_velocities > speed_thres);
% mean_velo = nanmean(WT_list_velocities(WT_list_velocities>speed_thres));
% median_velo = nanmedian(WT_list_velocities(WT_list_velocities>speed_thres));
% 
% 
% WT_OverallMeanSpeed  = nanmean(mean_velo)
% WT_OverallMedianSpeed = nanmedian(median_velo)
% %mean_velo = mean_velo(~isnan(mean_velo))
% 
% %Translate Data column colNum to Data Source File
% %colNum = 1;
% %WT_list_velocities(colNum,:)
% %velo_wildtype.Properties.VarNames{colNum+1} 
% 
% export(WT_velodataset,'File','WT_velocities.csv','WriteVarNames',true);
% 
% dlmwrite('WT_velocities.csv',[0 nanmean(WT_list_velocities)],'delimiter','\t','-append');
% dlmwrite('WT_velocities.csv',[0 nanmedian(WT_list_velocities)],'delimiter','\t','-append');
% 
% 
% %% ABv1016
% 
% % Import Data Sets
% strDir = strcat(baseDir,'/v1016 AbetaArc/');
% [ABv1016_Dat,N,srcDir] = readDataForAnalysis(strDir);
% 
% 
% ABv1016_velodataset = selectFeature(ABv1016_Dat, 'velo');
% 
% ABv1016_list_velocities = double(ABv1016_velodataset(:,2:end));
% 
% %Remove Short Tracklets from Dataset
% vNanCountInColumns = sum(~isnan(ABv1016_list_velocities),1);
% colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
% ABv1016_velodataset(:,colIndexToRemove+1) = [];
% %Update The Matrix Of Tracklets
% ABv1016_list_velocities = double(ABv1016_velodataset(:,2:end));
% 
% 
% % list_velocities = list_velocities(tbl_velocities > speed_thres);
% mean_velo = nanmean(ABv1016_list_velocities(ABv1016_list_velocities>speed_thres));
% median_velo = nanmedian(ABv1016_list_velocities(ABv1016_list_velocities>speed_thres));
% 
% ABv1016_OverallMeanSpeed = nanmean(mean_velo)
% ABv1016_OverallMedianSpeed = nanmedian(median_velo)
% %mean_velo = mean_velo(~isnan(mean_velo))
% 
% export(ABv1016_velodataset,'File','ABv1016_velocities.csv','WriteVarNames',true);
% 
% dlmwrite('ABv1016_velocities.csv',[0 nanmean(ABv1016_list_velocities)],'delimiter','\t','-append');
% dlmwrite('ABv1016_velocities.csv',[0 nanmedian(ABv1016_list_velocities)],'delimiter','\t','-append');
% 
% 
% %% v1190
% 
% %import data
% strDir = strcat(baseDir,'/v1190 APPTau/');
% [v1190_Dat,N,srcDir] = readDataForAnalysis(strDir);
% 
% 
% v1190_velodataset = selectFeature(v1190_Dat, 'velo');
% 
% v1190_list_velocities = double(v1190_velodataset(:,2:end));
% 
% %Remove Short Tracklets from Dataset
% vNanCountInColumns = sum(~isnan(v1190_list_velocities),1);
% colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
% v1190_velodataset(:,colIndexToRemove+1) = [];
% %Update The Matrix Of Tracklets
% v1190_list_velocities = double(v1190_velodataset(:,2:end));
% 
% 
% 
% % list_velocities = list_velocities(tbl_velocities > speed_thres);
% mean_velo = nanmean(v1190_list_velocities(v1190_list_velocities>speed_thres));
% median_velo = nanmedian(v1190_list_velocities(v1190_list_velocities>speed_thres));
% 
% v1190_OverallMeanSpeed = nanmean(mean_velo)
% v1190_OverallMedianSpeed = nanmedian(median_velo)
% %mean_velo = mean_velo(~isnan(mean_velo))
% 
% export(v1190_velodataset,'File','v1190_velocities.csv','WriteVarNames',true);
% 
% dlmwrite('v1190_velocities.csv',[0 nanmean(v1190_list_velocities)],'delimiter','\t','-append');
% dlmwrite('v1190_velocities.csv',[0 nanmedian(v1190_list_velocities)],'delimiter','\t','-append');
% 
% %% BWD46 wt
% %import data
% strDir = strcat(baseDir,'/BWD46 Abeta(wt)/');
% [BWD46_Dat,N,srcDir] = readDataForAnalysis(strDir);
% 
% 
% BWD46_velodataset = selectFeature(BWD46_Dat, 'velo');
% 
% BWD46_list_velocities = double(BWD46_velodataset(:,2:end));
% 
% %Remove Short Tracklets from Dataset
% vNanCountInColumns = sum(~isnan(BWD46_list_velocities),1);
% colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
% BWD46_velodataset(:,colIndexToRemove+1) = [];
% %Update The Matrix Of Tracklets
% BWD46_list_velocities = double(BWD46_velodataset(:,2:end));
% 
% 
% 
% % list_velocities = list_velocities(tbl_velocities > speed_thres);
% mean_velo = nanmean(BWD46_list_velocities(BWD46_list_velocities>speed_thres));
% median_velo = nanmedian(BWD46_list_velocities(BWD46_list_velocities>speed_thres));
% 
% BWD46_OverallMeanSpeed = nanmean(mean_velo)
% BWD46_OverallMedianSpeed = nanmedian(median_velo)
% %mean_velo = mean_velo(~isnan(mean_velo))
% 
% export(BWD46_velodataset,'File','BWD46_velocities.csv','WriteVarNames',true);
% 
% dlmwrite('BWD46_velocities.csv',[0 nanmean(BWD46_list_velocities)],'delimiter','\t','-append');
% dlmwrite('BWD46_velocities.csv',[0 nanmedian(BWD46_list_velocities)],'delimiter','\t','-append');
% 
% 
% %% BWD47(Abeta_Ita)
% %import data
% strDir = strcat(baseDir,'/BWD47(Abeta_Ita)/');
% [BWD47_Dat,N,srcDir] = readDataForAnalysis(strDir);
% 
% 
% BWD47_velodataset = selectFeature(BWD47_Dat, 'velo');
% 
% BWD47_list_velocities = double(BWD47_velodataset(:,2:end));
% 
% %Remove Short Tracklets from Dataset
% vNanCountInColumns = sum(~isnan(BWD47_list_velocities),1);
% colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
% BWD47_velodataset(:,colIndexToRemove+1) = [];
% %Update The Matrix Of Tracklets
% BWD47_list_velocities = double(BWD47_velodataset(:,2:end));
% 
% 
% 
% % list_velocities = list_velocities(tbl_velocities > speed_thres);
% mean_velo = nanmean(BWD47_list_velocities(BWD47_list_velocities>speed_thres));
% median_velo = nanmedian(BWD47_list_velocities(BWD47_list_velocities>speed_thres));
% 
% BWD47_OverallMeanSpeed = nanmean(mean_velo)
% BWD47_OverallMedianSpeed = nanmedian(median_velo)
% %mean_velo = mean_velo(~isnan(mean_velo))
% 
% export(BWD47_velodataset,'File','BWD47_velocities.csv','WriteVarNames',true);
% 
% dlmwrite('BWD47_velocities.csv',[0 nanmean(BWD47_list_velocities)],'delimiter','\t','-append');
% dlmwrite('BWD47_velocities.csv',[0 nanmedian(BWD47_list_velocities)],'delimiter','\t','-append');
% 
% 
% 
% %% BWD48.2(AbetaArc)on attP2
% %import data
% strDir = strcat(baseDir,'/BWD48.2(AbetaArc)on attP2/');
% [BWD482_Dat,N,srcDir] = readDataForAnalysis(strDir);
% 
% 
% BWD482_velodataset = selectFeature(BWD482_Dat, 'velo');
% 
% BWD482_list_velocities = double(BWD482_velodataset(:,2:end));
% 
% %Remove Short Tracklets from Dataset
% vNanCountInColumns = sum(~isnan(BWD482_list_velocities),1);
% colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
% BWD482_velodataset(:,colIndexToRemove+1) = [];
% %Update The Matrix Of Tracklets
% BWD482_list_velocities = double(BWD482_velodataset(:,2:end));
% 
% 
% 
% % list_velocities = list_velocities(tbl_velocities > speed_thres);
% mean_velo = nanmean(BWD482_list_velocities(BWD482_list_velocities>speed_thres));
% median_velo = nanmedian(BWD482_list_velocities(BWD482_list_velocities>speed_thres));
% 
% BWD482_OverallMeanSpeed = nanmean(mean_velo)
% BWD482_OverallMedianSpeed = nanmedian(median_velo)
% %mean_velo = mean_velo(~isnan(mean_velo))
% 
% export(BWD482_velodataset,'File','BWD482_velocities.csv','WriteVarNames',true);
% 
% dlmwrite('BWD482_velocities.csv',[0 nanmean(BWD482_list_velocities)],'delimiter','\t','-append');
% dlmwrite('BWD482_velocities.csv',[0 nanmedian(BWD482_list_velocities)],'delimiter','\t','-append');
% 
% 
% 
% %%~~~~
% %% BWD48(AbetaArc)
% %import data
% strDir = strcat(baseDir,'/BWD48(AbetaArc)/');
% [BWD48_Dat,N,srcDir] = readDataForAnalysis(strDir);
% 
% 
% BWD48_velodataset = selectFeature(BWD48_Dat, 'velo');
% 
% BWD48_list_velocities = double(BWD48_velodataset(:,2:end));
% 
% %Remove Short Tracklets from Dataset
% vNanCountInColumns = sum(~isnan(BWD48_list_velocities),1);
% colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
% BWD48_velodataset(:,colIndexToRemove+1) = [];
% %Update The Matrix Of Tracklets
% BWD48_list_velocities = double(BWD48_velodataset(:,2:end));
% 
% 
% 
% % list_velocities = list_velocities(tbl_velocities > speed_thres);
% mean_velo = nanmean(BWD48_list_velocities(BWD48_list_velocities>speed_thres));
% median_velo = nanmedian(BWD48_list_velocities(BWD48_list_velocities>speed_thres));
% 
% BWD48_OverallMeanSpeed = nanmean(mean_velo)
% BWD48_OverallMedianSpeed = nanmedian(median_velo)
% %mean_velo = mean_velo(~isnan(mean_velo))
% 
% export(BWD48_velodataset,'File','BWD48_velocities.csv','WriteVarNames',true);
% 
% dlmwrite('BWD48_velocities.csv',[0 nanmean(BWD48_list_velocities)],'delimiter','\t','-append');
% dlmwrite('BWD48_velocities.csv',[0 nanmedian(BWD48_list_velocities)],'delimiter','\t','-append');
% 
% %% attP2xC155 (gen ctr BWD.2)
% %import data
% strDir = strcat(baseDir,'/attP2xC155 (gen ctr BWD.2)/');
% [attP2xC155_Dat,N,srcDir] = readDataForAnalysis(strDir);
% 
% 
% attP2xC155_velodataset = selectFeature(attP2xC155_Dat, 'velo');
% 
% attP2xC155_list_velocities = double(attP2xC155_velodataset(:,2:end));
% 
% %Remove Short Tracklets from Dataset
% vNanCountInColumns = sum(~isnan(attP2xC155_list_velocities),1);
% colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
% attP2xC155_velodataset(:,colIndexToRemove+1) = [];
% %Update The Matrix Of Tracklets
% attP2xC155_list_velocities = double(attP2xC155_velodataset(:,2:end));
% 
% 
% 
% % list_velocities = list_velocities(tbl_velocities > speed_thres);
% mean_velo = nanmean(attP2xC155_list_velocities(attP2xC155_list_velocities>speed_thres));
% median_velo = nanmedian(attP2xC155_list_velocities(attP2xC155_list_velocities>speed_thres));
% 
% attP2xC155_OverallMeanSpeed = nanmean(mean_velo)
% attP2xC155_OverallMedianSpeed = nanmedian(median_velo)
% %mean_velo = mean_velo(~isnan(mean_velo))
% 
% export(attP2xC155_velodataset,'File','attP2xC155_velocities.csv','WriteVarNames',true);
% 
% dlmwrite('attP2xC155_velocities.csv',[0 nanmean(attP2xC155_list_velocities)],'delimiter','\t','-append');
% dlmwrite('attP2xC155_velocities.csv',[0 nanmedian(attP2xC155_list_velocities)],'delimiter','\t','-append');
% 
% %% attP40xC155 (gen ctr BWD)
% %import data
% strDir = strcat(baseDir,'/attP40xC155 (gen ctr BWD)/');
% [attP40xC155_Dat,N,srcDir] = readDataForAnalysis(strDir);
% 
% %Dataset and array - use array for calculations / 
% attP40xC155_velodataset = selectFeature(attP40xC155_Dat, 'velo');
% attP40xC155_list_velocities = double(attP40xC155_velodataset(:,2:end));
% 
% %Remove Short Tracklets from Dataset
% vNanCountInColumns = sum(~isnan(attP40xC155_list_velocities),1);
% colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
% attP40xC155_velodataset(:,colIndexToRemove+1) = [];
% %Update The Matrix Of Tracklets
% attP40xC155_list_velocities = double(attP40xC155_velodataset(:,2:end));
% 
% %Remove Low Speeds
% [cellRow,cellCol] = find(attP40xC155_list_velocities <= speed_thres);
% %Only way I could find now is just Iterate through Indexes set to NaN
% for i=1:size(cellRow) attP40xC155_velodataset{cellRow(i),cellCol(i)+1} = NaN; end
% 
% %Update The Matrix Of Tracklets
% attP40xC155_list_velocities = double(attP40xC155_velodataset(:,2:end));
% runVelo = attP40xC155_list_velocities;
% dataset = attP40xC155_velodataset;
% 
% mean_velo = nanmean(runVelo);
% median_velo = nanmedian(runVelo);
% %mean_velo = mean_velo(~isnan(mean_velo))
% 
% filename = sprintf('attP40xC155_velocities-gt%d.csv',speed_thres);
% export(dataset,'File',filename,'WriteVarNames',true);
% 
% dlmwrite(filename,[0 mean_velo],'delimiter','\t','-append');
% dlmwrite(filename,[0 median_velo],'delimiter','\t','-append');
% 
% nanmean(mean_velo)
% 
% %%%
% 
% %% BWD46 wt TEST
% %import data
% % 
% % strDir = strcat(baseDir,'/BWD46 Abeta(wt)/FT_20160111_BWD46_Abeta(wt)_6/');
% % [BWD46T_Dat,N,srcDir] = readDataForAnalysis(strDir);
% % 
% % 
% % BWD46T_velodataset = selectFeature(BWD46T_Dat, 'velo');
% % 
% % BWD46T_list_velocities = double(BWD46T_velodataset(:,2:end));
% % 
% % %Remove Short Tracklets from Dataset
% % vNanCountInColumns = sum(~isnan(BWD46T_list_velocities),1);
% % colIndexToRemove = find(vNanCountInColumns < minTrackletLength);
% % 
% % BWD46T_velodataset(:,colIndexToRemove+1) = [];
% % %Update The Matrix Of Tracklets
% % BWD46T_list_velocities = double(BWD46T_velodataset(:,2:end));
% % 
% % 
% % % list_velocities = list_velocities(tbl_velocities > speed_thres);
% % mean_velo = nanmean(BWD46T_list_velocities(BWD46T_list_velocities>speed_thres));
% % median_velo = nanmedian(BWD46T_list_velocities(BWD46T_list_velocities>speed_thres));
% % 
% % BWD46T_OverallMeanSpeed = nanmean(mean_velo)
% % BWD46T_OverallMedianSpeed = nanmedian(median_velo)
% % %mean_velo = mean_velo(~isnan(mean_velo))
% % 
% % export(BWD46T_velodataset,'File','BWD46T_velocities.csv','WriteVarNames',true);
% % 
% % dlmwrite('BWD46T_velocities.csv',[0 nanmean(BWD46T_list_velocities)],'delimiter','\t','-append');
% % dlmwrite('BWD46T_velocities.csv',[0 nanmedian(BWD46T_list_velocities)],'delimiter','\t','-append');
% % 
% % 
% 
% 
% 
% 
% % %% attP40
% % velo_attP40 = selectFeature(attP40_Dat, 'velo');
% % 
% % attP40_list_velocities = double(velo_attP40(:,2:end));
% % % list_velocities = list_velocities(tbl_velocities > speed_thres);
% % mean_velo = nanmean(attP40_list_velocities(attP40_list_velocities>speed_thres));
% % median_velo = nanmedian(attP40_list_velocities(attP40_list_velocities>speed_thres));
% % 
% % attP40_OverallMeanSpeed = nanmean(mean_velo)
% % attP40_OverallMedianSpeed = nanmedian(median_velo)
% % %mean_velo = mean_velo(~isnan(mean_velo))
% % velo_attP40_exp = join(velo_attP40,mean_velo,median_velo)
% % export(velo_attP40,'File','attP40_velocities.csv','WriteVarNames',true);
% % 
% % dlmwrite('attP40_velocities.csv',[0 nanmean(attP40_list_velocities)],'delimiter','\t','-append');
% % dlmwrite('attP40_velocities.csv',[0 nanmedian(attP40_list_velocities)],'delimiter','\t','-append');

