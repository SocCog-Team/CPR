addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/preprocessing
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab

cd /Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_social_context/
% cd /Users/fschneider/Desktop/social_context_study/

% Number of tested dyads
num_dyads = 5;

% Configuration file
cfg_pth = '/Users/fschneider/Documents/GitHub/CPR/Matlab/CFG/felix_context_study.cfg';

% Variables to import
var_import = {
    'INFO_', ...
    'TRIAL_', ...
    'IO_joystickDirection', ...
    'IO_joystickStrength',...
    'IO_fixation_flag',...
    'EYE_x_dva',...
    'EYE_y_dva',...
    '#stimDisplay'};

% Loop through each dyad folder...
for iDyad = 4:num_dyads
    % Extract all files
    file_list = dir(['Dyad' num2str(iDyad) '/mwk2/']);
    
    % Convert each CPR .mwk2 file to .h5 
    for iFile = 1:length(file_list)
        if contains(file_list(iFile).name, 'CPR')
            d = CPR_import_mwk2([file_list(iFile).folder '/' file_list(iFile).name], var_import, true, cfg_pth);
        end
    end
end