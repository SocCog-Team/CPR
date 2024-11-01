function out = CPR_main_online_v3(data,varargin)
% This function interfaces with MWorks and provides online stimulus  
% information for the CPR paradigm. MWorks needs to call this function in
% the Matlab Interface Window. Don't forget to select the variables needed 
% for the creation of the param structure.
%
% Input:  	.data           MWorks online input
%
% Output:   .out            Logical, required output 
%
% Example:	d = CPR_main_online(data)
%
% Known bugs:
%
% Feature wish list:
%
% Version history
%   1.0     (fxs 2024-04-18) Initial version.
%   1.1     (fxs 2024-10-30) Random walk implemented. Variables renamed.

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

% Import MWorks data 
if isfield(data, 'time') && isfield(data, 'value') && isfield(data, 'event')
    d                           = data;
else
    d                           = MW_readOnline(data,'~checkTrialBorder','debugLevel',0);
end

% Extract parameter
param.trial                     = cell2mat(d.value(d.event == 'TRIAL_start'));                  % Trial number
param.cycle_duration_ms      	= cell2mat(d.value(d.event == 'TMP_cycle_duration_ms'));    	% Duration of random walk (stimulus cycle)
param.coh_duration_ms       	= cell2mat(d.value(d.event == 'TMP_coh_block_duration_ms'));  	% Duration of coherence block
param.feedback_probability     	= cell2mat(d.value(d.event == 'TMP_feedback_probability')); 	% Probability of reward target appearance
param.min_feedback_interval_ms 	= cell2mat(d.value(d.event == 'TMP_feedback_ITI_ms'));      	% Interval between consecutive feedback presentations
param.Fs                        = 1000 / 120;                                                   % Screen sampling rate
param.pth                       = pth;
tmp_snr                         = d.value(d.event == 'TMP_snr_list');                           % Stimulus coherence list
param.snr_list                  = cellfun(@double, tmp_snr{1});

% Create and write new stimulus cycle
STIM                            = CPR_create_random_walk_v3(param);        	% Draw RDP stimulus parameters
[~]                             = CPR_write_txt(STIM,pth);                  % Write parameters to .txt files
out                             = true;
toc
end