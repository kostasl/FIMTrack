function [Dataset, numberOfSubDirs, larvaPerSubDir] = readDataForAnalysis(rootPath)

Dataset = dataset();
larvaID = 0;
larvaPerSubDir = [];


%% get all CSV-Files
[CSVFiles numberOfSubDirs] = getAllFiles(rootPath);

%%
for i=1:size(CSVFiles,1)
    larvaID = 0;
    %% Initialize variables.
    filename = CSVFiles{i};
    fnameparts = strsplit(filename,filesep);
    
    delimiter = ',';
    startRow = 2;
    
    %% print file name to console
    filename
    
    %% Open the text file.
    fileID = fopen(filename,'r');
    lineIjustReadIn = fgetl(fileID);
    numColumns = sum(lineIjustReadIn == delimiter);
    larvaPerSubDir = [larvaPerSubDir numColumns];
    fseek(fileID, 0, -1);
    
    %% Format string for each line of text:
    formatSpec = ['%s' repmat('%f', 1, numColumns) '%[^\n\r]'];
    
    %% Read columns of data according to format string.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    
    %% Close the text file.
    fclose(fileID);
           
    %% Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post
    % processing code is included. To generate code which works for
    % unimportable data, select unimportable cells in a file and regenerate the
    % script.
    
    %cellLarva       = repmat({strcat(fnameparts(end-2),fnameparts(end-1),'/larva')}, 1, numColumns);
    %The following trick is required so isStrings works within the
    %setVarNames of the dataset
    tt1 = fnameparts(end-2);
    tt2 = fnameparts(end-1);
    cellLarva       = repmat({strcat(tt1{1},'_',tt2{1},'_larva')}, 1, numColumns);
    cellLarvaID     = cell(1, numColumns);
    
    for i=1:numColumns
        cellLarvaID{i} = larvaID;%Make new Column ID
        larvaID = larvaID + 1;
    end
    cellLarvaIDFinal    = cellfun(@(x,y) strcat(x,num2str(y)), cellLarva, cellLarvaID, 'uniformoutput', 0);
    cellLarvaIDFinal    = [cellstr('ID') cellLarvaIDFinal];
    
    %I fixed the string so its suitable for VarNames - Otherwise it causes
    %trouble
    cellLarvaIDFinal = strrep(cellLarvaIDFinal, 'output', '');
    cellLarvaIDFinal = strrep(cellLarvaIDFinal, '(', '');
    cellLarvaIDFinal = strrep(cellLarvaIDFinal, ')', '');
    cellLarvaIDFinal = strrep(cellLarvaIDFinal, '-', '_');
    
    %display(cellLarvaIDFinal);
    %% Create output variable
    Data = dataset(dataArray{1:end-1}, 'VarNames', cellLarvaIDFinal);
    if(isempty(Dataset))
        Dataset = Data;
    else
        % IMPORTANT: join(X,Y,'ID') needs an appropriate ID entry in X
        % (e.g. 'mom_x(1001)') so that X must be the dataset with more
        % entries (i.e. with more time steps). This is guaranteed by
        % the following if-condition. However, if join(Y,X,'ID') must
        % be used (else case), the larval variable names are not in
        % correct order anymore! This is fixed by the sort section
        % below...
        %%KL: Join two datasets , but restrict size to the smallest one in
        %%terms of number of frames - so join can work
        if (size(Data.ID,1) >= size(Dataset.ID,1))
            Dataset = join(Dataset, Data, 'ID');
        else
            Dataset = join(Data,Dataset,'ID');
        end
    end
end

%% remove animals containing only NaNs
% since all dataset entries (i.e. larvae) are cutted so that they have the
% same length, some animals from very long movies only exist in the frames
% which was cutted. The result is an animal in the dataset which consists
% of NaNs only (all meaningful entries are cut away). These animals need to
% be removed since there is no useful data in them.
% larva_names = Dataset.Properties.VarNames(2:end);
% collision_indicator = selectFeature(Dataset,'in_collision');
% larva_names_to_delete = cell(1,0);
% for k=1:numel(larva_names)
%     % get collision indicator for the larvae
%     l_collision = collision_indicator.(larva_names{k});
%     
%     % find first frame number in which the larva exists
%     is_not_nan_idx = find(~isnan(l_collision)==1);
%     
%     if(numel(is_not_nan_idx) == 0)
%        larva_names_to_delete{1,end+1}=larva_names{k};
%     end
% end
% 
% for k=1:numel(larva_names_to_delete)
%    Dataset.(larva_names_to_delete{k}) = []; 
% end
% NOTE ABOVE CODE REMOVED AS THERE IS NO COLLISSION INDICATOR in_collision

% %% Sort columns in talbe ('larva0', 'larva1', 'larva2', ...)
% 
% % larva names (e.g. 'larva101') might not be in order anymore
% % get the larva names ('ID', 'larva101', 'larva1334', ...)
% larva_names = Dataset.Properties.VarNames;
% % get the numbers from the larva names (101, 1334, ...)
% larva_numbers = cellfun(@str2num, cellfun(@(x) [regexprep(x,'larva','')], larva_names(2:end), 'UniformOutput', false));
% % sort the numbers and save the permutation (isort)
% [~,isort]=sort(larva_numbers);
% % add a 0 to the permutation and increment all values (because of the
% % 'ID' variable name which must be in the first column
% adjusted_isort = [0,isort]+1;
% % reorder the dataset based on the permutation
% Dataset = Dataset(:,adjusted_isort);
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans Data CSVFiles cellLarva cellLarvaID cellLarvaIDFinal larvaID;
end

