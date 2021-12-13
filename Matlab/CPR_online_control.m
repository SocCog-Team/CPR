function out = CPR_online_control(data,varargin)
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
%   1.0     (fxs 2020-11-01) Initial version.

tic
setup                           = 'MAC40656';

if strcmp(setup, 'MAC40656')
    pth                         = '/Users/fschneider/Desktop/';
    addpath '/Users/fschneider/Documents/GitHub/CPR/Matlab'
    addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
    addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
else
    pth                         = '/Volumes/CPR/';
    addpath('/Volumes/CPR/');
    addpath /Users/cnl/Desktop/Felix/online
    addpath /Users/cnl/Desktop/Felix/online/CircStat2012a
    addpath /Users/cnl/ownCloud/Shared/MatLab4MWorks
end

% Import data

d                               = MW_readOnline(data,'~checkTrialBorder','debugLevel',0);

% Extract parameter
param.trial                     = cell2mat(d.value(d.event == 'TRIAL_start'));
param.agent_flag                = cell2mat(d.value(d.event == 'TMP_show_agent'));
param.NoStates                  = cell2mat(d.value(d.event == 'TMP_NoStates'));
param.NoCoherenceStates         = cell2mat(d.value(d.event == 'TMP_NoCoherenceStates'));
param.state_min_ms              = cell2mat(d.value(d.event == 'TMP_state_min_ms'));
param.state_max_ms              = cell2mat(d.value(d.event == 'TMP_state_max_ms'));
tmp_snr                         = d.value(d.event == 'TMP_snr_list');
tmp_ddir                        = d.value(d.event == 'TMP_directionChange_list');
param.snr_list                  = cellfun(@double, tmp_snr{1});
param.directionChange_list      = cellfun(@double, tmp_ddir{1});
param.target_ITI_ms             = cell2mat(d.value(d.event == 'TMP_target_ITI_ms'));
param.target_blocked_ms         = cell2mat(d.value(d.event == 'TMP_target_ban_duration_ms'));
param.target_prob               = .005;
param.Fs                        = 100;
param.pth                       = pth;

% Draw RDP stimulus parameters
STATES                          = CPR_stimulus_control(param);

% Prepare agent response
if param.agent_flag == true
    AGNT.dir_sigma          	= 30;                                       % Direction sigma [deg]
    AGNT.str_sigma            	= .1;                                       % Strength sigma [%]
    AGNT.lag                  	= 50;                                       % Delay to reference point [samples]
    AGNT.win                 	= 50;                                       % Smoothing window size [samples]
    AGNT.smooth_kernel        	= 'gaussian';                               % Smoothing kernel [samples]
    AGNT                      	= CPR_prepare_agent(STATES,AGNT);           % Compute agent
    [~]                         = CPR_update_params_txt(STATES, AGNT, pth);	% Write parameters to .txt files
else
    [~]                         = CPR_update_params_txt(STATES,[],pth);   	% Write parameters to .txt files
end

out                             = true;
toc

end