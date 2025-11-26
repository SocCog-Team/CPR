function out = CPR_main_online_v5(data,varargin)
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
%   1.2     (fxs 2024-10-30) Created random walk agent.
%   1.3     (fxs 2025-09-02) Added replay option.

tic
setup = 'MAC40656';

if strcmp(setup, 'MAC40656')
    pth = '/Volumes/cnl/Desktop/Felix/CPR/';
    addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/PSY_Setup/random_walk/
    addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
    addpath(genpath('/Users/fschneider/Documents/GitLab/matlab4mworks/'));

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
param.agent_flag                = cell2mat(d.value(d.event == 'TMP_show_agent'));
param.replay_flag               = cell2mat(d.value(d.event == 'TMP_show_replay'));
param.Fs                        = 1000 / 120;                                                   % Screen sampling rate
param.pth                       = pth;
tmp_snr                         = d.value(d.event == 'TMP_snr_list');                           % Stimulus coherence list
param.snr_list                  = cellfun(@double, tmp_snr{1});

% Create and write new stimulus cycle
STIM                            = CPR_create_random_walk_v3(param);        	% Draw RDP stimulus parameters

% Prepare agent response
if param.agent_flag == true
    if param.replay_flag == true
        solo_summary_file           = '/Users/cnl/Desktop/Felix/online/random_walk/replay_data/20250812_fih_psy4.mat'; % better player
        [STIM, AGNT]             	= CPR_create_replay_agnt(solo_summary_file, STIM, param.trial);
    else
        %         AGNT.dir_sigma          	= 15;                                       % Direction sigma [deg]
        %         AGNT.str_sigma            	= .025;                                     % Eccentricity sigma [%]
        %         AGNT.lag                  	= 50;                                       % Delay to reference point [samples]
        %         AGNT.win                 	= 50;                                       % Smoothing window size [samples]
        %         AGNT.smooth_kernel        	= 'gaussian';                               % Smoothing kernel [samples]
        %         AGNT                        = CPR_create_agent_random_walk(STIM,AGNT);
        
        onnx_path               = '/Users/fschneider/Documents/GitHub/CPR/Matlab/PSY_Setup/random_walk/lstm_player.onnx';
        coherence               = repmat(.8, length(STIM.RDP_direction_deg),1);
        acc_mean                = .5;
        tlt_mean                = .5;
        noise                   = 0;
        use_ts                  = false;
        [AGNT.dir_smooth,  AGNT.str_smooth] = pred_player(onnx_path, STIM.RDP_direction_deg, coherence, acc_mean, tlt_mean, noise, use_ts);
        AGNT
    end
    [~]                       	= CPR_write_txt(STIM,AGNT, pth);            % Write parameters to .txt files
else
    [~]                       	= CPR_write_txt(STIM,[],pth);               % Write parameters to .txt files
end

out                             = true;
toc
end