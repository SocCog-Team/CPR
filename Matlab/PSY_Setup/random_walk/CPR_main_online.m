function out = CPR_main_online(data,varargin)
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
setup = 'PSY4';

if strcmp(setup, 'MAC40656')
    pth = '/Users/fschneider/Desktop/';
    addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/PSY_Setup/random_walk
    addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
    addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
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

% Extract parameter
param.trial                     = cell2mat(d.value(d.event == 'TRIAL_start'));                  % Trial number
param.walk_duration_ms      	= cell2mat(d.value(d.event == 'TMP_cycle_duration_ms'));    	% Duration of random walk (stimulus cycle)
param.coh_duration_ms       	= cell2mat(d.value(d.event == 'TMP_coh_block_duration_ms'));  	% Duration of coherence block
param.polar_step_size       	= cell2mat(d.value(d.event == 'TMP_polar_step_size'));          % Step size of random walk in polar space
param.reward_probability       	= cell2mat(d.value(d.event == 'TMP_feedback_probability')); 	% Probability of reward target appearance
param.feedback_interval_ms      = cell2mat(d.value(d.event == 'TMP_feedback_ITI_ms'));      	% Interval between feedback presentation
param.Fs                        = 1000 / 120;                                                   % Screen sampling rate
param.pth                       = pth;
tmp_snr                         = d.value(d.event == 'TMP_snr_list');                           % Stimulus coherence list
param.snr_list                  = cellfun(@double, tmp_snr{1});

STIM                            = CPR_create_random_walk(param);            % Draw RDP stimulus parameters
[~]                             = CPR_write_txt(STIM,pth);                  % Write parameters to .txt files
out                             = true;
toc

end