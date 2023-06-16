
function updateTirFiles(dropDown)
    % Specify the folder path where the .tir files are located
    folderPath = '.\Simulation\Tires';

    % Get a list of all .tir files in the folder and its subfolders
    fileList = dir(fullfile(folderPath, '**/*.tir'));
    
    % Extract only the filenames without extensions
    fileNames = {fileList.name};
    
    dropDown.Items = fileNames;
end