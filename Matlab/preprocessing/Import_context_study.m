addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/preprocessing
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab

% cd /Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_social_context/
cd /Users/fschneider/Desktop/social_context_study/

% Number of tested dyads
num_dyads = 18;

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
        'IO_joystickDirection2', ...
        'IO_joystickStrength2',...
        'IO_fixation2_flag',...
        'EYE_x2_dva',...
        'EYE_y2_dva',...
        '#stimDisplay'};

% Loop through each dyad folder...
for iDyad = 8:num_dyads
    % Extract all files
    file_list = dir(['Dyad' num2str(iDyad) '/mwk2/']);
    
    % Convert each CPR .mwk2 file to .h5 
    for iFile = 1:length(file_list)
        if contains(file_list(iFile).name, 'CPR')
            d = CPR_import_mwk2([file_list(iFile).folder '/' file_list(iFile).name], var_import, true, cfg_pth);
        end
    end
end

%% Manually

% file_path = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_social_context/Solo/AnE/20240705_ane_CPRsolo_block2_psy4.mwk2';
% file_path = '/Users/fschneider/Desktop/social_context_study/Solo/NiA/mwk2/20240815_nia_CPRsolo_block1_psy3.h5';
% d = CPR_import_mwk2(file_path, var_import, true, cfg_pth);

%% Solo
solo_dir            = '/Users/fschneider/Desktop/social_context_study/Solo/';
% solo_dir            = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_social_context/Solo/';
cfg_pth             = '/Users/fschneider/Documents/GitHub/CPR/Matlab/CFG/felix_context_study.cfg';
sbj_lst             = dir(solo_dir);
cnt                 = 0;

var_import = {
    'INFO_', ...
    'TRIAL_', ...
    'IO_joystickDirection', ...
    'IO_joystickStrength',...
    'IO_fixation_flag',...
    'EYE_x_dva',...
    'EYE_y_dva',...
    '#stimDisplay'};

for iSubj = 1:length(sbj_lst)
    if length(sbj_lst(iSubj).name) ~= 3
        continue
    end
    
    file_list = dir([solo_dir sbj_lst(iSubj).name '/mwk2/']);

    % Convert each CPR .mwk2 file to .h5 
    for iFile = 1:length(file_list)
        if contains(file_list(iFile).name, 'CPR')
            d = CPR_import_mwk2([file_list(iFile).folder '/' file_list(iFile).name], var_import, true, cfg_pth);
        end
    end
end