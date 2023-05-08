clear all
close all

% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/WIP/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

% Copy subject summary from server
source_pth              = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/';
local_pth               = '/Users/fschneider/Documents/CPR_psychophysics/';
dest_dir                = '/Users/fschneider/Desktop/CPRsolo_Tobias/';
fname                   = 'Subjects_summary.xlsx';
[status,msg,msgID]      = copyfile([source_pth fname], local_pth);

% Import subject summary table
x                       = readtable([local_pth fname]);
sbj_lst                 = x.Abbreviation;
sbj_lst(cellfun(@isempty,sbj_lst)) = [];

% Preprocess data
for iSubj = 1:length(sbj_lst)
    
    if ~isfolder([local_pth sbj_lst{iSubj}])
        continue
    end
    
    disp(['Processing subject: ' sbj_lst{iSubj}])
                
    % Change into data source directory
    cd([local_pth sbj_lst{iSubj} '/summary/'])

    % Extract .mwk2 file names from source directory
    mat_files              = dir('*.mat');    
    
    for iFile = 1:length(mat_files)
        if contains(mat_files(iFile).name,'CPRsolo')
            [status,msg,msgID]	= copyfile([local_pth sbj_lst{iSubj} '/summary/' mat_files(iFile).name], dest_dir);
        end
    end
end