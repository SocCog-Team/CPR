clear all
close all

addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitLab/matlab4mworks/
addpath /Users/fschneider/Documents/GitLab/matlab4mworks/code_maybe_needed/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/mat_to_summary/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/

pth = '/Volumes/T7_Shield/';
fname = '20260126_mah_CPRsolo_block2_psy3_sao.mwk2';
cfg_pth = '/Users/fschneider/Documents/GitHub/CPR/Matlab/CFG/felix_random_walk.cfg';
import_flag = true;

% Define relevant variables
var_lst = {
    'INFO_', ...
    'IO_joystickDirection', ...
    'IO_joystickStrength',...
    'IO_joystickDirection2', ...
    'IO_joystickStrength2',...
    'STIM',...
    'TRIAL_',...
    '#stimDisplay'};

if import_flag
    % Import .MWK2 and write as .H5
    d = MW_readData([pth fname], 'include', var_lst, '~typeOutcomeCheck');
    MW_writeH5(d, [pth fname(1:end-5) '.h5'], 'replace', 'privateCFG', cfg_pth)
else
    % Import .H5
    d = MW_readData([pth fname(1:end-5) '.h5'], 'include', var_lst, '~typeOutcomeCheck');
end

%%
% Index variables of interest
idx                         = [];
idx.cOn                     = d.event == 'TRIAL_start';
idx.cOn                     = d.event == 'TRIAL_start';
idx.cEnd                    = d.event == 'TRIAL_end';
idx.frame                   = d.event == 'STIM_displayUpdate';
idx.RDP_onset               = d.event == 'STIM_RDP_onset';
idx.RDP_dir                 = d.event == 'STIM_RDP_direction';
idx.RDP_coh                 = d.event == 'STIM_RDP_coherence';
idx.RDP_dot                	= d.event == 'STIM_RDP_dotPositions';
idx.JS_dir                  = d.event == 'IO_joystickDirection';
idx.JS_str                  = d.event == 'IO_joystickStrength';
idx.JS2_dir                 = d.event == 'IO_joystickDirection2';
idx.JS2_str                 = d.event == 'IO_joystickStrength2';
idx.fixation              	= d.event == 'IO_fixation_flag';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.outcome2                = d.event == 'TRIAL_outcome2';
idx.feedback                = d.event == 'TRIAL_reactionTrigger';
idx.cum_score            	= d.event == 'INFO_Score';
idx.cum_score2              = d.event == 'INFO_Score2';
idx.feedback_score          = d.event == 'INFO_FeedbackScore';
idx.feedback_score2         = d.event == 'INFO_FeedbackScore2';
idx.bonus1                  = d.event == 'INFO_bonus_ply1_cents';
idx.bonus2                  = d.event == 'INFO_bonus_ply2_cents';
idx.performance             = d.event == 'INFO_performance_percent';
idx.performance2            = d.event == 'INFO_performance_percent2';
idx.task                    = d.event == 'INFO_task';

out                         = response_readout(d, idx);
save([pth 'replay_data_' fname(1:end-9) '.mat'], 'out', '-v7.3')

%%%%%%%%%%%%%%%%%
%%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%

function out = response_readout(d, idx)

% Get cycle timestamps
cyc.cOn                     = d.time(idx.cOn);
cyc.cEnd                    = d.time(idx.cEnd);

last_entry = 0;
last_rew = 0;

for iCyc = 1:length(cyc.cEnd)
    % Trial index
    cycIdx                  = [];
    cycIdx                  = d.time >= cyc.cOn(iCyc) & d.time <= cyc.cEnd(iCyc);

    % Task info
    task                    =  getTrialData(d.value, cycIdx, idx.task);

    % Trial flag
    if contains(task, 'CPR_solo')
        out.raw.solo_flag(iCyc) = true;
    elseif contains(task, 'CPR_dyadic')
        out.raw.solo_flag(iCyc) = false;
    else
        out.raw.solo_flag(iCyc) = 999;
        continue
    end

    % Extract frame times
    frame_ts             	= getTrialData(d.time, cycIdx, idx.frame);                      % Frame timestamps

    % Extract coherence blocks
    tmp_rdp_coh             = getTrialData(d.value, cycIdx, idx.RDP_coh);                   % RDP coherence
    tmp_rdp_coh_ts          = getTrialData(d.time, cycIdx, idx.RDP_coh);                    % RDP coherence timestamps

    % Extract stimulus direction
    tmp_rdp_dir            	= getTrialData(d.value, cycIdx, idx.RDP_dir);                   % RDP_direction
    tmp_rdp_dir_ts        	= getTrialData(d.time, cycIdx, idx.RDP_dir);                    % RDP_direction timestamps

    % Joystick data
    tmp_js_dir_ts        	= getTrialData(d.time, cycIdx, idx.JS_dir);                     % JS player1
    tmp_js_dir            	= getTrialData(d.value, cycIdx, idx.JS_dir);
    tmp_js_tlt           	= getTrialData(d.value, cycIdx, idx.JS_str);
    
    % Target onset
    clear tmp_fbk_ts tmp_fbk_val feedback outcome outcome2 outcome_ts outcome2_ts
    tmp_fbk_ts              = getTrialData(d.time, cycIdx, idx.feedback);
    tmp_fbk_val             = getTrialData(d.value, cycIdx, idx.feedback);

    if ~isnan(tmp_fbk_val)
        feedback                = tmp_fbk_ts(logical(tmp_fbk_val));
        outcome                 = getTrialData(d.value, cycIdx, idx.outcome);
        outcome_ts              = getTrialData(d.time, cycIdx, idx.outcome);

        if length(feedback) ~= length(outcome) && ~strcmp(outcome{end},'FixationBreak')
            warning(['Cyc:' num2str(iCyc) '| nFeedback ~= nOutcome'])
        end

    else
        feedback                  = [];
        outcome                 = [];
        outcome2                = [];
    end

               
    % Initialise variables
    clear rdp_dir js_dir js2_dir js_err js2_err js_tlt js2_tlt js_acc js2_acc js_dff rdp_dff
    rdp_dir                 = nan(1,length(frame_ts));
    rdp_coh                 = nan(1,length(frame_ts));
    js_dir                  = nan(1,length(frame_ts));
    js_tlt                  = nan(1,length(frame_ts));

    % Create frame-wise data
    for iFrame = 1:length(frame_ts)
        if sum(tmp_rdp_dir_ts < frame_ts(iFrame)) ~=0
            rdp_dir(iFrame) = tmp_rdp_dir(find(tmp_rdp_dir_ts < frame_ts(iFrame),1,'last'));
            rdp_coh(iFrame) = tmp_rdp_coh(find(tmp_rdp_coh_ts < frame_ts(iFrame),1,'last'));
            js_dir(iFrame)  = tmp_js_dir(find(tmp_js_dir_ts < frame_ts(iFrame),1,'last'));
            js_tlt(iFrame)  = tmp_js_tlt(find(tmp_js_dir_ts < frame_ts(iFrame),1,'last'));
        end
    end

    % Calculate accuracy
    js_err                  = rad2deg(circ_dist(deg2rad(js_dir),deg2rad(rdp_dir)));         % Get circular distance to RDP direction
    js_acc                	= abs(1 - abs(js_err / 180));                                   % Calculate accuracy

    % Sample-by-sample difference (Derivative)
    js_dff                  = [0 rad2deg(circ_dist(deg2rad(js_dir(1:end-1)),deg2rad(js_dir(2:end))))];
    rdp_dff                 = [0 rad2deg(circ_dist(deg2rad(rdp_dir(1:end-1)),deg2rad(rdp_dir(2:end))))];
    
    ex                    	= isnan(js_dff) | isnan(rdp_dff);
    js_dff(ex)            	= 0;
    rdp_dff(ex)            	= 0;

    % Save raw signals
    out.raw.ts{iCyc}     	= frame_ts;
    out.raw.rdp_dir{iCyc} 	= rdp_dir;
    out.raw.rdp_coh{iCyc} 	= rdp_coh;
    out.raw.rdp_dff{iCyc} 	= rdp_dff;
    out.raw.js_dir{iCyc}  	= js_dir;
    out.raw.js_dff{iCyc}  	= js_dff;
    out.raw.js_err{iCyc}   	= abs(js_err);
    out.raw.js_acc{iCyc}   	= js_acc;
    out.raw.js_tlt{iCyc}   	= js_tlt;
    out.raw.fbk_on_ts{iCyc} = feedback;
    if ~isempty(feedback)
        out.raw.fbk_on_smpl{iCyc} = ceil((double(feedback - cyc.cOn(iCyc)) ./1e3) ./ (1000/120)); %subtract cycleOn, convert to ms, divide by ms/frame, round up
    else
        out.raw.fbk_on_smpl{iCyc} = [];
    end
    out.raw.outcome{iCyc}   = outcome;

end
end
