function [tbl,d,summ,crr] = CPR_psych_import(pth, fname, d)

% This function imports, organises and visualises CPR data.
%
% Input:  	.subj           String, Subject ID, 3 letter code
%           .pth            String, Data directory
%           .fname        	String, File identifier
%          	.d              Structure (optional), Contains trial
%                           information. Skips data import if included.
%
% Output:   .tbl            Cell, Contains state-wise data table with
%                           stimulus parameters as well as behavioural
%                           responses
%          	.d              Structure, Contains (raw) MWorks data
%          	.summ           Cell, Contains structure with summary
%                           information
%
% Example:	[tbl,d,summ] = CPR_psych_import('fxs', '/Users/fschneider/Documents', 'test.mwk2')
%
% Known bugs:
%
% Feature wish list:
%
% Version history
%   1.0     (fxs 2020-09-01) Initial version.
%   1.1     (fxs 2021-05-25) Added compatibility with dyadic setting.
%   1.1     (fxs 2021-05-25) Analysis of multiple targets per state
%                            possible.
%   1.2     (fxs 2021-11-24) Frame-wise import & analysis    

%% Add relevant directories

addpath /Users/fschneider/Documents/GitHub/CPR/Matlab
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

cd(pth)

%% Import data

% Define variables for import
if contains(fname,'dyadic')
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
elseif contains(fname,'agent')
    var_import = {
        'INFO_', ...
        'TRIAL_', ...
        'IO_joystickDirection', ...
        'IO_joystickStrength',...
        'IO_fixation_flag',...
        'EYE_x_dva',...
        'EYE_y_dva',...
        'AGNT_direction', ...
        'AGNT_strength',...
        '#stimDisplay'};
else
    var_import = {
        'INFO_', ...
        'TRIAL_', ...
        'IO_joystickDirection', ...
        'IO_joystickStrength',...
        'IO_fixation_flag',...
        'EYE_x_dva',...
        'EYE_y_dva',...
        '#stimDisplay'};
end

% Import H5 instead on MWK2
write = false; 

% Get file identifier/recording date
if iscell(fname)
    % Do something here
else
    fid = strsplit(fname,'_');
    
    if length(fid{2}) > 3
        sbj{1}                 	= fid{2}(1:3);
        sbj{2}                 	= fid{2}(4:end);
    else
        sbj{1}                 	= fid{2};
    end
    
    if contains(fname,'agent')
        sbj{2} = 'agnt';
        write = true;
    end
end

% Import .mwk2 data file
if nargin < 3
    d                	= CPR_import_mwk2(fname, var_import, write);
    
%     d                	= CPR_import_mwk2(fname, var_import, true);
%     d                   = CPR_data_correction(d, 'IO_joystickDirection', 'IO_joystickStrength');    % Correct for sample differences
%     d                   = CPR_data_correction(d, 'IO_joystickDirection2', 'IO_joystickStrength2');   
end

%% Extract data and organise in table

% Initialise variables
idx                         = [];
cyc                         = [];

% Create variable-specific indices
idx.cOn                     = d.event == 'TRIAL_start';
idx.cEnd                    = d.event == 'TRIAL_end';
idx.frame                   = d.event == 'STIM_displayUpdate';
idx.RDP_onset               = d.event == 'STIM_RDP_onset';
idx.RDP_dir                 = d.event == 'STIM_RDP_direction';
idx.RDP_coh                 = d.event == 'STIM_RDP_coherence';
idx.RDP_dot                	= d.event == 'STIM_RDP_dotPositions';
idx.trg_on                  = d.event == 'STIM_target_onset';
idx.JS_dir                  = d.event == 'IO_joystickDirection';
idx.JS_str                  = d.event == 'IO_joystickStrength';
idx.fixation              	= d.event == 'IO_fixation_flag';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.trg                     = d.event == 'TRIAL_reactionTrigger';
idx.eye_x_dva              	= d.event == 'EYE_x_dva';
idx.eye_y_dva             	= d.event == 'EYE_y_dva';

if strcmp(sbj{1}, 'cla') || strcmp(sbj{1}, 'nil')
    idx.reward              = d.event == 'INFO_Juice_ml';
    rew_str                 = 'ml';
else
    idx.reward              = d.event == 'INFO_Score';
    rew_str                 = '[a.u.]';
end

if contains(fname,'dyadic')
    idx.JS2_dir          	= d.event == 'IO_joystickDirection2';
    idx.JS2_str         	= d.event == 'IO_joystickStrength2';
    idx.fixation2       	= d.event == 'IO_fixation2_flag';
    idx.outcome2          	= d.event == 'TRIAL_outcome2';
    idx.eye_x2_dva       	= d.event == 'EYE_x2_dva';
    idx.eye_y2_dva       	= d.event == 'EYE_y2_dva';
    idx.reward2           	= d.event == 'INFO_Score2';
elseif contains(fname,'agent')
    idx.JS2_dir          	= d.event == 'AGNT_direction';
    idx.JS2_str         	= d.event == 'AGNT_strength';
    idx.fixation2       	= zeros(1,length(d.event));
    idx.outcome2          	= d.event == 'TRIAL_outcome2';
    idx.eye_x2_dva       	= zeros(1,length(d.event));
    idx.eye_y2_dva       	= zeros(1,length(d.event));
    idx.reward2           	= d.event == 'INFO_Score2';      
end

% Get trial timestamps
cyc.cOn                     = d.time(idx.cOn);
cyc.cEnd                    = d.time(idx.cEnd);

% Check if subject ID is dyad
if contains(fname,'dyadic') && length(fid{2}) < 6
    error('Dyatic experiment. Use concatenated 3-letter subject code: ''xxxyyy''')
end

% Construct data table(s)
tbl{1}                      = CPR_construct_table(sbj{1}, fname, d, cyc, idx, true);
    
if size(sbj,2) > 1 && isfield(idx,'JS2_dir')
    % Define fields that index variables of subject 2
    fields = {'JS2_dir',...
        'JS2_str',...
        'outcome2',...
        'fixation2',...
        'eye_x2_dva',...
        'eye_y2_dva',...
        'reward2'};
    
    % Create new structure. Use subj1 fieldnames but fill with subj2
    % variable index. Necessary to construct similar table for both subjects.
    idx2                    = rmfield(idx,fields);
    idx2.JS_dir          	= idx.JS2_dir;
    idx2.JS_str         	= idx.JS2_str;
    idx2.outcome          	= idx.outcome2 ;
    idx2.fixation       	= idx.fixation2;
    idx2.eye_x_dva       	= idx.eye_x2_dva;
    idx2.eye_y_dva       	= idx.eye_y2_dva;
    idx2.reward           	= idx.reward2;
    
    % Construct data table for subject2
    tbl{2}               	= CPR_construct_table(sbj{2}, fname, d, cyc, idx2, true);
end

% Plot performance summary
% [summ,crr]                  = CPR_plot_performance_summary(d,idx,tbl,rew_str,sbj,fname); % buggy
% [~]                         = CPR_target_reward(tbl{1});
summ = []; crr = [];
end
