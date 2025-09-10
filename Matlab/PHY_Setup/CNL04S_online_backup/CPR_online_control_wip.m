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
% Example:	out = CPR_online_control(data)
%
% Known bugs:
%
% Feature wish list:
%
% Versionhistory
%   1.0     (fxs 2020-11-01) Initial version.
%   1.1     (fxs 2023-04-16) Receptive field processing
%   1.2     (fxs 2025-02-14) Receptive field optimization

tic
setup                           = 'PHYSIO4';

global RF_map RF_tun f ff f_tun ax_lst nRep stim_pos

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
if isfield(data, 'time') && isfield(data, 'value') && isfield(data, 'event')
    d                           = data;
else
    d                         	= MW_readOnline(data,'~checkTrialBorder','debugLevel',0);
    d.time                    	= double(d.time);
end

% Check experimental task
x                               = d.value(d.event == "INFO_task");
INFO_task                       = x{1};

% Receptive field mapping
if  strcmp(INFO_task, "RF_mapping")
    if isempty(f) || ~ishandle(f)
        f                           = figure('units','normalized','position',[0 0 1 1]);
        
        for iChannel = 1:32
            ax_lst(iChannel)     	= subplot(4,8,iChannel);
        end
        
        ff                           = figure('units','normalized','position',[0 0 1 1]);
        for iChannel = 33:64
            ax_lst(iChannel)     	= subplot(4,8,iChannel-32);
        end

    end
    
    param.trial               	= cell2mat(d.value(d.event == 'TRIAL_start'));
    if param.trial <= 1
        RF_map                  = [];
        nRep                    = [];
        stim_pos                = [];
        [RF_map,nRep]           = RF_estimation_online_copy(d,f,ff,ax_lst);
    else
        % [RF_map,nRep]        	= RF_estimation_online(d,f,ax_lst,RF_map,nRep);
        [RF_map,nRep,stim_pos]        	= RF_estimation_online_copy(d,f,ff,ax_lst,RF_map,nRep,stim_pos);
    end
    
% Receptive field tuning
elseif  strcmp(INFO_task, "RF_tuning")
    
    if isempty(f_tun) | ~ishandle(f_tun)
        f_tun                       = figure('units','normalized','position',[0 0 1 1]);
    end
    
    param.trial               	= cell2mat(d.value(d.event == 'TRIAL_start'));
    if param.trial <= 1
        RF_tun                  = [];
        RF_tun               	= RF_tuning_online(d,f_tun);
    else
        RF_tun                	= RF_tuning_online(d,f_tun,RF_tun);
    end
    
% Continuous perceptual report
elseif strcmp(INFO_task, "CPR_solo_stepfunction_neutral") || strcmp(INFO_task, "CPR_dyadic_stepfunction_neutral") || strcmp(INFO_task, "CPR_dyadic_computer_stepfunction_neutral")
    
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
        AGNT.dir_sigma          	= 45;                                       % Direction sigma [deg]
        AGNT.str_sigma            	= .2;                                       % Eccentricity sigma [%]
        AGNT.str                    = [];                                       % Eccentricity [norm]
        AGNT.lag                  	= 50;                                       % Delay to reference point [samples]
        AGNT.win                 	= 50;                                       % Smoothing window size [samples]
        AGNT.smooth_kernel        	= 'gaussian';                               % Smoothing kernel [samples]
        AGNT                      	= CPR_prepare_agent(STATES,AGNT);           % Compute agent
        [~]                         = CPR_update_params_txt(STATES, AGNT, pth);	% Write parameters to .txt files
    else
        [~]                         = CPR_update_params_txt(STATES,[],pth);   	% Write parameters to .txt files
    end
end

out                             = true;
toc

end