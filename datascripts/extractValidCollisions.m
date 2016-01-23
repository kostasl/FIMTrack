function [valid_collisions_wt, short_collisions_wt, overlaps_wt] = extractValidCollisions(data, parameters)

collision_indicator = selectFeature(data,'in_collision');
gfp_indicator = selectFeature(data, 'is_gfp_larva');
area = selectFeature(data,'area');
is_pusher = selectFeature(data,'isPusher');

% head_x = selectFeature(data,'head_x');
% head_y = selectFeature(data,'head_y');
% sp_1_x = selectFeature(data,'spinepoint_1_x');
% sp_1_y = selectFeature(data,'spinepoint_1_y');
% sp_2_x = selectFeature(data,'spinepoint_2_x');
% sp_2_y = selectFeature(data,'spinepoint_2_y');
% sp_3_x = selectFeature(data,'spinepoint_3_x');
% sp_3_y = selectFeature(data,'spinepoint_3_y');
% tail_x = selectFeature(data,'tail_x');
% tail_y = selectFeature(data,'tail_y');

larva_names = data.Properties.VarNames(2:end);

min_time_BA     = parameters{1};
min_time_col    = parameters{2};
acc_gap_size    = parameters{3};
min_larva_area  = parameters{4};

%% return statements
% valid_collisions is a 2xN cell where N is the number of
% larvae which are involved in a valid collisions. In the first row the
% larvaIDs are given. The second row is a is a 3xM cell again where M is
% the number of collisions this particular larva is involved in
% EXAMPLE: valid_collision_living
%   'larva0'    'larva3'    'larva5'    'larva10'
%   [3x1]       [3x3]       [3x1]       [3x3]
% so that four larvae are involved in collisions.
% larva0 is involved in one collision, larva3 is involved in 3 collisions.
% EXAMPLE valid_collisions for larva0:
%   {[1;422];               [423;438];          [439;8140]}
%   valid frames before c.  collision frames    valid frames after coll.
valid_collisions_wt = cell(2,0);
valid_collisions_gfp = cell(2,0);

short_collisions_wt = cell(2,0);
short_collisions_gfp = cell(2,0);

overlaps_wt = cell(2,0);
overlaps_gfp = cell(2,0);

%% for all larvae
for k=1:numel(larva_names)
    
    %% get features from the larvae
    % get collision indicator for the larvae
    l_collision = collision_indicator.(larva_names{k});
    % is larva gfp larva?
    l_gfp_indicator = gfp_indicator.(larva_names{k});
    % get larval area
    l_area = area.(larva_names{k});
    % get pusher indicator
    l_is_pusher = is_pusher.(larva_names{k});
    
    % find first frame number in which the larva exists
    is_not_nan_idx = find(~isnan(l_collision)==1);
    
    % find first NaN
    start_to_exist = is_not_nan_idx(1);
    ends_to_exist = is_not_nan_idx(end);
    
    
    %% check for short_collisions
    % get time points (indices) in which the animal is in a collision
    in_collision_idx = find(l_collision == 1);
    
    % if there are collisions
    if(~isempty(in_collision_idx))
        short_collisions     = cell(1,0);
        % special case: there is only one short touch (for one frame):
        % numel(in_collision_idx)==1:
        if (numel(in_collision_idx) == 1)
            before = [start_to_exist-1, in_collision_idx-1];
            during = [in_collision_idx, in_collision_idx];
            after = [in_collision_idx+1, ends_to_exist+1];
            
            % get collision meta data:
            % do the larva survives the collision?
            if(after(1) <= size(l_collision, 1))
                survives_collision = ~isnan(l_collision(after(1)));
            else
                survives_collision = logical(0);
            end
            
            % distance to previous collision
            if (before(1) ~= 1) % if this is not the first frame
                dist_prev_collision = during(1) - before(1);
            else
                dist_prev_collision = 0;
            end
            
            % distance to next collision
            if (survives_collision)
                dist_next_collision = after(2) - during(2);
            else
                dist_next_collision = 0;
            end
            
            % get the is pusher flag for the first collision frame
            is_pusher_flag = l_is_pusher(during(1));
            % -----------------------------------------------------------
            
            if (collision_length < min_time_col &&...
                    dist_prev_collision > 2*min_time_col &&...
                    dist_next_collision > 2*min_time_col &&...
                    survives_collision == 1) %&& is_pusher_flag)
                short_collisions{1,end+1} = during;
            end
            
            % else: find collision indices
        else
            dIdx = diff(in_collision_idx)==1;
            
            % there might be two short collisions at t and t+2. In this
            % case dIdx = 0 (no consecutive index difference = 1) and f
            % would be empty... For now: ignore those cases.
            if (~(numel(dIdx) == 1 && dIdx == 0))
                
                f = find([0;dIdx] ~= [dIdx;0]);
                
                if(1 ~= in_collision_idx(f(1)))
                    %                 colIndex = [0; in_collision_idx(f)];
                    colIndex = [start_to_exist-1; in_collision_idx(f)];
                else
                    colIndex = in_collision_idx(f);
                    colIndex(0) = colIndex(0)-1;
                end
                
                if(size(l_collision, 1) ~= in_collision_idx(f(end)))
                    %                 colIndex = [colIndex; size(l_collision, 1)+1];
                    colIndex = [colIndex; ends_to_exist+1];
                else
                    colIndex(end) = colIndex(end)+1;
                end
                % analyze collision indices
                for i=1:2:size(colIndex,1)-3
                    % get frames before collision (e.g. [start, c1_start-1];
                    % important: frame 0 will be set to frame 1...
                    before  = [colIndex(i)+1; colIndex(i+1)-1];
                    
                    % get start and stop frame of the collision
                    % (e.g. [c1_start, c1_stop])
                    during  = [colIndex(i+1); colIndex(i+2)];
                    
                    % get frames after collision (e.g. [c1_end+1 c2_start-1])
                    after   = [colIndex(i+2)+1; colIndex(i+3)-1];
                    
                    % -----------------------------------------------------------
                    % get collision meta data:
                    % do the larva survives the collision?
                    if(after(1) <= size(l_collision, 1))
                        survives_collision = ~isnan(l_collision(after(1)));
                    else
                        survives_collision = logical(0);
                    end
                    
                    % length of the collision
                    collision_length = during(2) - during(1);
                    
                    % length of time before collision
                    before_length = before(2) - before(1);
                    
                    % length of time after collision
                    if (survives_collision)
                        after_length = after(2) - after(1);
                    else
                        after_length = 0;
                    end
                    
                    % distance to previous collision
                    if (before(1) ~= 1) % if this is not the first frame
                        dist_prev_collision = during(1) - before(1);
                    else
                        dist_prev_collision = 0;
                    end
                    
                    % distance to next collision
                    if (survives_collision)
                        dist_next_collision = after(2) - during(2);
                    else
                        dist_next_collision = 0;
                    end
                    
                    % get the is pusher flag for the first collision frame
                    is_pusher_flag = l_is_pusher(during(1));
                    % -----------------------------------------------------------
                    
                    if (collision_length < min_time_col &&...
                            dist_prev_collision > 2*min_time_col &&...
                            dist_next_collision > 2*min_time_col &&...
                            survives_collision == 1) % && is_pusher_flag
                        short_collisions{1,end+1} = during;
                    end
                end
            end
        end
        
        if(~isempty(short_collisions))
            % is the larva a gfp larva?
            is_gfp_larva = l_gfp_indicator(before(2));
            if (is_gfp_larva == 1)
                short_collisions_gfp{1,end+1}  = larva_names{k};
                short_collisions_gfp{2,end}    = short_collisions;
            else
                short_collisions_wt{1,end+1}  = larva_names{k};
                short_collisions_wt{2,end}    = short_collisions;
            end
        end
        
        
        
        %% check for valid collisions and overlaps
        % smooth collision indicator
        l_collision = fillCollisionGaps(l_collision,acc_gap_size);
        % get time points (indices) in which the animal is in a collision
        in_collision_idx = find(l_collision == 1);
        
        % if there are collisions
        if(~isempty(in_collision_idx) && numel(in_collision_idx) > 1)
            valid_collisions    = cell(3,0);
            overlaps            = cell(1,0);
            
            dIdx = diff(in_collision_idx)==1;
            
            % only several short collisions and no consecutive indices in 
            % which this animal is in a collision (ignore these cases)
            if (~(numel(dIdx) == 1 && dIdx == 0))
                
                f = find([0;dIdx] ~= [dIdx;0]);
                
                if(1 ~= in_collision_idx(f(1)))
                    %                 colIndex = [0; in_collision_idx(f)];
                    colIndex = [start_to_exist-1; in_collision_idx(f)];
                else
                    colIndex = in_collision_idx(f);
                    colIndex(0) = colIndex(0)-1;
                end
                
                if(size(l_collision, 1) ~= in_collision_idx(f(end)))
                    %                 colIndex = [colIndex; size(l_collision, 1)+1];
                    colIndex = [colIndex; ends_to_exist+1];
                else
                    colIndex(end) = colIndex(end)+1;
                end
                
                % analyze collision indices
                for i=1:2:size(colIndex,1)-3
                    % get frames before collision (e.g. [start, c1_start-1];
                    % important: frame 0 will be set to frame 1...
                    before  = [colIndex(i)+1; colIndex(i+1)-1];
                    
                    % get start and stop frame of the collision
                    % (e.g. [c1_start, c1_stop])
                    during  = [colIndex(i+1); colIndex(i+2)];
                    
                    % get frames after collision (e.g. [c1_end+1 c2_start-1])
                    after   = [colIndex(i+2)+1; colIndex(i+3)-1];
                    
                    % -----------------------------------------------------------
                    % get collision meta data:
                    % do the larva survives the collision?
                    survives_collision = ~isnan(l_collision(after(1)));
                    
                    % length of the collision
                    collision_length = during(2) - during(1);
                    
                    % length of time before collision
                    before_length = before(2) - before(1);
                    
                    % length of time after collision
                    if (survives_collision)
                        after_length = after(2) - after(1);
                    else
                        after_length = 0;
                    end
                    
                    % distance to previous collision
                    if (before(1) ~= 1) % if this is not the first frame
                        dist_prev_collision = during(1) - before(1);
                    else
                        dist_prev_collision = 0;
                    end
                    
                    % distance to next collision
                    if (survives_collision)
                        dist_next_collision = after(2) - during(2);
                    else
                        dist_next_collision = 0;
                    end
                    
                    % is the area below the area min threshold during
                    % collision?
                    below_min_size_threshold = l_area(during(1):during(2)) < min_larva_area;
                    valid_area_during_collision = sum(below_min_size_threshold) == 0;
                    
                    % get the is pusher flag for the first collision frame
                    is_pusher_flag = l_is_pusher(during(1));
                    % -----------------------------------------------------------
                    
                    if (collision_length >= min_time_col &&...
                            before_length >= min_time_BA &&...
                            after_length >= min_time_BA &&...
                            survives_collision == 1 &&...
                            valid_area_during_collision == 1) %&& is_pusher_flag
                        valid_collisions{1,end+1} = before;
                        valid_collisions{2,end}   = during;
                        valid_collisions{3,end}   = after;
                    elseif (~valid_area_during_collision)
                        overlaps{1,end+1} = during;
                    end
                end
                
                if(~isempty(valid_collisions))
                    % is the larva a gfp larva?
                    is_gfp_larva = l_gfp_indicator(before(2));
                    if (is_gfp_larva == 1)
                        valid_collisions_gfp{1,end+1}  = larva_names{k};
                        valid_collisions_gfp{2,end}    = valid_collisions;
                    else
                        valid_collisions_wt{1,end+1}  = larva_names{k};
                        valid_collisions_wt{2,end}    = valid_collisions;
                    end
                end
                
                if (~isempty(overlaps))
                    % is the larva a gfp larva?
                    is_gfp_larva = l_gfp_indicator(before(2));
                    if (is_gfp_larva == 1)
                        overlaps_gfp{1,end+1}  = larva_names{k};
                        overlaps_gfp{2,end}    = overlaps;
                    else
                        overlaps_wt{1,end+1}  = larva_names{k};
                        overlaps_wt{2,end}    = overlaps;
                    end
                end
            end
        end
    end
end
end