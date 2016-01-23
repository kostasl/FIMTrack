function [fileList numberOfSubDirs] = getAllFiles(dirName)
    
    numberOfSubDirs = 0;
    
    % Get the data for the current directory
    dirData = dir(dirName);
    
    % Find the index for directories
    dirIndex = [dirData.isdir];
    
    fileList = {dirData(~dirIndex).name}';  %' Get a list of the files
    idx = ~cellfun('isempty',regexp(fileList,'.csv'));
    fileList = fileList(idx);
    if ~isempty(fileList)
        fileList = cellfun(@(x) fullfile(dirName,x),...  % Prepend path to files
                       fileList,'UniformOutput', false);
    end
    
    % Get a list of the subdirectories
    subDirs = {dirData(dirIndex).name};
    
    % Find index of subdirectories that are not '.' or '..'
    validIndex = ~ismember(subDirs,{'.','..', 'Plots'});
    
    for iDir = find(validIndex)                         % Loop over valid subdirectories
        nextDir = fullfile(dirName,subDirs{iDir});      % Get the subdirectory path
        fileList = [fileList; getAllFiles(nextDir)];    % Recursively call getAllFiles
        numberOfSubDirs = numberOfSubDirs + 1;
    end
end

