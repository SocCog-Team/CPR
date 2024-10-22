function out = CPR_main_online_v2(data,varargin)
% This function interfaces with MWorks and provides online stimulus  
% information for the CPR paradigm. MWorks needs to call this function in
% the Matlab Interface Window. Don't forget to select the variables needed 
% for the creation of the param structure.
%
% Input:  	.data           MWorks online input
%
% Output:   .out            Logical, required output 
%
% Example:	d = CPR_import_mwk2('fname', [], false)
%
% Known bugs:
%
% Feature wish list:
%
% Version history
%   1.0     (fxs 2024-04-18) Initial version.

tic
setup = 'MAC40656';

if strcmp(setup, 'MAC40656')
    pth = '/Users/fschneider/Desktop/';
    addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/PSY_Setup/random_walk/
    addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
    addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
else
    pth = '/Volumes/CPR/';
    addpath('/Volumes/CPR/');
    addpath /Users/cnl/Desktop/Felix/online/random_walk
    addpath /Users/cnl/Desktop/Felix/online/CircStat2012a
    addpath /Users/cnl/ownCloud/Shared/MatLab4MWorks
end

% Import data
if isfield(data, 'time') && isfield(data, 'value') && isfield(data, 'event')
    d                           = data;
else
    d                           = MW_readOnline(data,'~checkTrialBorder','debugLevel',0);
end

% Set parameters
param.trial                     = cell2mat(d.value(d.event == 'TRIAL_start'));% Trial number
param.coh_duration_ms       	= 10000;                % Duration of coherence block
param.walk_duration_ms      	= 60000;                % Duration of random walk
param.Fs                        = 1000 / 120;           % Screen sampling rate [Hz]
param.snr_list                  = [.2 .4 .6 .8];        % Stimulus coherence
param.jump_probability          = 0.0025;               % Probability of stimulus direction jump
param.min_jump_interval_ms      = param.Fs * 100;       % Min interval between direction jumps
param.feedback_probability    	= 0.0025;               % Probability of reward
param.min_feedback_interval_ms  = param.Fs * 100;       % Min interval between reward administration

STIM                            = CPR_create_random_walk(param);            % Draw RDP stimulus parameters
[~]                             = CPR_write_txt(STIM,pth);                  % Write parameters to .txt files
out                             = true;
toc

end