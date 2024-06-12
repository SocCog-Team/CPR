close all
clear all

addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
    
file_pth                    = '/Users/fschneider/Desktop/Social_context_pilot/';
file_name                   = '20240514_CPRneutral_hekfxs_test.h5';
d                        	= MW_readH5([file_pth file_name]); % ...load .h5 file

%%

file_pth                    = '/Users/fschneider/Desktop/Social_context_pilot/';
file_name                   = '20240514_CPRsolo_hek_test.h5';
d                        	= MW_readH5([file_pth file_name]); % ...load .h5 file
myHandle = figure;
myTL = MW_timelineOLD_RAB([], 'setHandle', myHandle);
myTL = MW_timelineOLD_RAB(myTL, 'addData', d, 'plot','trialNumber',2);

 %%

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
idx.JS2_dir                 = d.event == 'IO_joystickDirection2';
idx.JS2_str                 = d.event == 'IO_joystickStrength2';
idx.fixation              	= d.event == 'IO_fixation_flag';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.trg                     = d.event == 'TRIAL_reactionTrigger';

% Get trial timestamps
cyc.cOn                     = d.time(idx.cOn);
cyc.cEnd                    = d.time(idx.cEnd);

figure;hold on
title('Target timestamps')

for iCyc = 1:length(cyc.cEnd)-1
    % Trial index
    cycIdx                  = [];
    cycIdx                  = d.time >= cyc.cOn(iCyc) & d.time <= cyc.cEnd(iCyc);
    
    % Extract dots
    dot_position           	= d.value(idx.RDP_dot);
    dot_position_ts       	= d.time(idx.RDP_dot);
    frame_ts                = d.time(idx.frame);
    
    % Stimulus data
    tmp_dir_ts              = getTrialData(d.time, cycIdx, idx.RDP_dir);                    % RDP_direction timestamps
    tmp_dir                 = getTrialData(d.value, cycIdx, idx.RDP_dir);                   % RDP_direction
    tmp_coh                 = getTrialData(d.value, cycIdx, idx.RDP_coh);                   % RDP coherence
    tmp_coh_ts            	= getTrialData(d.time, cycIdx, idx.RDP_coh);                    % RDP coherence timestamps
    
    tmp_target_val         	= getTrialData(d.value, cycIdx, idx.trg);                       % Target flag value
    tmp_target_ts         	= getTrialData(d.time, cycIdx, idx.trg);                        % Target timestamps
    tmp_targetOn_ts         = tmp_target_ts(tmp_target_val == 1);                          	% Target onset timestamps
    nTrg(iCyc)              = length(tmp_targetOn_ts);
    
    tmp_outcome          	= getTrialData(d.value, cycIdx, idx.outcome);                 	% Target outcome player1
    tmp_outcome_ts      	= getTrialData(d.time, cycIdx, idx.outcome);                 	% Target outcome timestamp

    % Joystick data
    tmp_joystick_dir     	= getTrialData(d.value, cycIdx, idx.JS_dir);                 	% JS player1
    tmp_joystick_str    	= getTrialData(d.value, cycIdx, idx.JS_str);                 
    tmp_joystick_dir_ts     = getTrialData(d.time, cycIdx, idx.JS_dir);                 	
 
    tmp_joystick2_dir     	= getTrialData(d.value, cycIdx, idx.JS2_dir);                 	% JS player2
    tmp_joystick2_str    	= getTrialData(d.value, cycIdx, idx.JS2_str);                 
    tmp_joystick2_dir_ts    = getTrialData(d.time, cycIdx, idx.JS2_dir);                 	
 
    scatter(tmp_targetOn_ts, repmat(iCyc,1,length(tmp_targetOn_ts)),'x')
end

%%% Check stuff
duration_cycle = (cyc.cEnd - cyc.cOn) ./ 1e6 % [s]
median_dot_ts = median(diff(dot_position_ts))
median_frame_ts = median(diff(frame_ts))
median_js_ts = median(diff(tmp_joystick_dir_ts))
num_targets = [min(nTrg) max(nTrg), mean(nTrg)]

%%% Plot %%%
figure
plot((tmp_joystick_dir_ts-tmp_joystick_dir_ts(1))/1e3,tmp_joystick_dir);
title('Joystick direction')

figure
cts = [tmp_coh_ts tmp_joystick_dir_ts(end)]./1e3;
plot((tmp_joystick_dir_ts-tmp_joystick_dir_ts(1))/1e3,tmp_joystick_str); hold on
stairs(cts-cts(1), [tmp_coh tmp_coh(end)])
title('Joystick eccentricity & RDP coherence')

figure;hold on
plot(frame_ts,[nan diff(frame_ts)])
scatter(cyc.cOn,ones(1,length(cyc.cOn)), 'rx')
title('Frame timestamps & Cycle onset')

figure;hold on
plot(dot_position_ts,[nan diff(dot_position_ts)])
scatter(cyc.cOn,ones(1,length(cyc.cOn)), 'rx')
title('Dot timestamps & Cycle onset')

figure
histogram(diff(tmp_dir_ts),10)
title('RDP ts difference')

figure
plot(diff(tmp_dir_ts))
title('RDP ts difference')
