function [ret_living, ret_dead] = analyzeCollisions_script(living, dead)

disp('************************************************')
disp('******** START COLLISION ANAYSIS SCRIPT ********')
disp('************************************************')
disp(['Start time: ' datestr(now)])

%% load tables based on paths

% all tables
% path_living =   '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/living/';
% path_dead =     '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/dead/';

path_living =   '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/alive/';
path_dead = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/dead/';

path_living = '/Volumes/Macintosh HD/Documents/Workspace/Matlab/Collision_Analysis/Tables/all csvs alive/';
path_dead = '/Volumes/Macintosh HD/Documents/Workspace/Matlab/Collision_Analysis/Tables/all csvs dead/';

% single tables
% path_living = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/living/NilsGFP4_overlay2/';
% path_living = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/living/NilsGFP4_overlay3/';
% path_living = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/living/NilsGFP4_overlay6/';
% path_living = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/living/NilsGFP4_overlay7/';
% path_living = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/living/NilsGFP5_overlay1/';
% path_living = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/living/NilsGFP5_overlay3/';
% path_living = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/living/NilsGFP5_overlay5/';

% path_dead = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/dead/NilsGFP4_overlay4_dead/';
% path_dead = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/dead/NilsGFP4_overlay5_dead/';
% path_dead = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/dead/NilsGFP4_overlay8_dead/';
% path_dead = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/dead/NilsGFP4_overlay9_dead/';
% path_dead = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/dead/NilsGFP5_overlay2_dead/';
% path_dead = '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/15_08_2014_Collision_Tables/dead/NilsGFP5_overlay4_dead/';

% path_living =   '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/alive/table_short/';
% path_living =   '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/alive/table2/';
% path_living =   '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/alive/table3/';
% path_living =   '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/alive/table4/';
% path_living =   '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/alive/table5/';
% path_living =   '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/alive/table6/';
% path_living =   '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/alive/table7/';

% path_dead =     '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/dead/table2/';
% path_dead =     '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/dead/table3/';
% path_dead =     '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/dead/table4/';
% path_dead =     '/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/Collision_Tables_Johannes/dead/table5_short/';


% living = readDataForAnalysis(path_living);
% dead = readDataForAnalysis(path_dead);




%% load tables based on mat files
% load('/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/New_Collision_Tables/living/living_gfp4_gfp5.mat');
% load('/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/New_Collision_Tables/dead/dead_gfp4_gfp5.mat');

 
% load('/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/New_Collision_Tables/living/living_gfp5.mat');
% load('/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/New_Collision_Tables/dead/dead_gfp5.mat');

% load('/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/New_Collision_Tables/living/living_gfp4.mat');
% load('/Users/bena/Documents/workspace/MATLAB/Collision_Analysis/New_Collision_Tables/dead/dead_gfp4.mat');

% LATEST:
% load('/Volumes/Macintosh HD/Documents/Workspace/Matlab/Collision_Analysis/Tables/all csvs alive/alive.mat')
% load('/Volumes/Macintosh HD/Documents/Workspace/Matlab/Collision_Analysis/Tables/all csvs dead/dead.mat')

%% print # of animals
disp('--------------------------------------')
disp(['#living animals: ' num2str(size(living,2))])
disp(['#dead animals: ' num2str(size(dead,2))])


%% set parameters
% min_time_before_after specifies the minimal time the larva must be visible before and
% after the collision
min_time_before_after = 20; %20
% min_time_coll specifies the miniamal time the collision must have
min_time_coll = 5;
% acc_gap_size specifies the acceptable gap size between two collision
% (if it is < min_gap_size, the collisions will be merged to one collision
acc_gap_size = 5; %10
% min_larva_area specifies the minimal larva area (during a collision)
min_larva_area = 800; %750 or 1100
% put everything in a parateters cell
parameters = {min_time_before_after, min_time_coll, acc_gap_size, min_larva_area};


%% extract valid, short and overlap collisions
[valid_collisions_wt_living, short_collisions_wt_living, overlaps_wt_living] = extractValidCollisions(living, parameters);
[valid_collisions_wt_dead, short_collisions_wt_dead, overlaps_wt_dead] = extractValidCollisions(dead, parameters);


%% print short collisions
disp('--------------------------------------')
disp(['#animals in short collisions (living): ' num2str(size(short_collisions_wt_living,2))])
disp(['#animals in short collisions (dead): ' num2str(size(short_collisions_wt_dead,2))])
tmp_living = cellfun(@(x) [size(x,2)], short_collisions_wt_living);
tmp_living = sum(tmp_living(2,:));
tmp_dead = cellfun(@(x) [size(x,2)], short_collisions_wt_dead);
tmp_dead = sum(tmp_dead(2,:));
disp(['#short collisions (living): ' num2str(tmp_living)])
disp(['#short collisions (dead): ' num2str(tmp_dead)])


%% print # of collisions
disp('--------------------------------------')
disp(['#animals in valid collisions (living): ' num2str(size(valid_collisions_wt_living,2))])
disp(['#animals in valid collisions (dead): ' num2str(size(valid_collisions_wt_dead,2))])
tmp_living = cellfun(@(x) [size(x,2)], valid_collisions_wt_living);
tmp_living = sum(tmp_living(2,:));
tmp_dead = cellfun(@(x) [size(x,2)], valid_collisions_wt_dead);
tmp_dead = sum(tmp_dead(2,:));
disp(['#valid collisions (living): ' num2str(tmp_living)])
disp(['#valid collisions (dead): ' num2str(tmp_dead)])


%% print # of overlaps
disp('--------------------------------------')
disp(['#animals in overlaps (living): ' num2str(size(overlaps_wt_living,2))])
disp(['#animals in overlaps (dead): ' num2str(size(overlaps_wt_dead,2))])
tmp_living = cellfun(@(x) [size(x,2)], overlaps_wt_living);
tmp_living = sum(tmp_living(2,:));
tmp_dead = cellfun(@(x) [size(x,2)], overlaps_wt_dead);
tmp_dead = sum(tmp_dead(2,:));
disp(['#overlap events (living): ' num2str(tmp_living)])
disp(['#overlap events (dead): ' num2str(tmp_dead)])


%% analyze valid collsisions
ret_living = analyzeValidCollisions(living,parameters,valid_collisions_wt_living);
ret_dead = analyzeValidCollisions(dead,parameters,valid_collisions_wt_dead);


%% generate plots

dead_velocities = selectFeature(dead, 'velosity');
living_velocities = selectFeature(living, 'velosity');

dead_isgfp = selectFeature(dead, 'is_gfp');
living_isgfp = selectFeature(living, 'is_gfp');

velocity_d = double(dead_velocities(:,2:end));
velocity_l = double(living_velocities(:,2:end));

isgfp_d = double(dead_isgfp(:,2:end));
isgfp_l = double(living_isgfp(:,2:end));

iswt_d = abs(isgfp_d -1);
iswt_l = abs(isgfp_l - 1);

iswt_l_nonan = iswt_l(~isnan(iswt_l));
iswt_d_nonan = iswt_d(~isnan(iswt_d));
velocity_d_nonan = velocity_d(~isnan(velocity_d));
velocity_l_nonan = velocity_l(~isnan(velocity_l));

velocity_d_wt = velocity_d_nonan(logical(iswt_d_nonan));
velocity_l_wt = velocity_l_nonan(logical(iswt_l_nonan));

velocity_d_wt_filtered = velocity_d_wt(velocity_d_wt<50);
velocity_l_wt_filtered = velocity_l_wt(velocity_l_wt<50);
velocity_l_wt_filtered = velocity_l_wt_filtered(velocity_l_wt_filtered>=0);
velocity_d_wt_filtered = velocity_d_wt_filtered(velocity_d_wt_filtered>=0);
v_d = velocity_d_wt_filtered;
v_l = velocity_l_wt_filtered(1:numel(velocity_d_wt_filtered));
v = [v_l, v_d];


 generateCollisionPlots(valid_collisions_wt_living,...
    valid_collisions_wt_dead,...
    short_collisions_wt_living,...
    short_collisions_wt_dead,...
    overlaps_wt_living,...
    overlaps_wt_dead,...
    ret_living,...
    ret_dead,v);


%% print finished time
disp('--------------------------------------')
disp(['Finished time: ' datestr(now)])
disp('************************************************')
end