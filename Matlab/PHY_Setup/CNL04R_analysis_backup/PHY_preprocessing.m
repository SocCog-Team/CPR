function [out] = PHY_preprocessing(fname,source_dir,dest_dir,cfg_pth,import_flag)
%
% This function calculates and visualises reaction times
%
% Input:        .fname          String, File name
%               .source_dir     String, Path to dara souce directory
%               .dest_dir       String, Path to file destination directory
%               .cfg_pth        String, Path to configuration .cfg file
%               .import_flag    Logical, Indicates if imported from .mwk2
% Output:       .out            Structure, Contains information about
%                               experiment, stimulus, and neural response
%
% Felix Schneider, CNL
%
% Version history
%   1.0     (fxs 2025-09-05) Initial version.
%
%%% TO DO %%%
% (1) Add LFP/MUAe signal
%

%% Adjust paths
addpath(genpath('/Users/cnl/Documents/GitLab/matlab4mworks/'));

mkw2Filename    = [source_dir 'mwk2/' fname '.mwk2'];
h5Filename      = [source_dir 'h5/' fname '.h5'];
pl2Filename     = [source_dir 'pl2/' fname '.pl2'];

%% Import MWorks data

close all
clear d

var_import  = {
    'INFO_', ...
    'TRIAL_', ...
    'STIM_', ...
    'IO_joystickDirection', ...
    'IO_joystickStrength',...
    'IO_fixation_flag',...
    'IO_sync_16bit',...
    'IO_syncWord',...
    'IO_rewardA_ml',...
    'EYE_x_dva',...
    'EYE_y_dva',...
    '#stimDisplayUpdate'};

if import_flag
    % Import .MWK2 and write as .H5
    d = MW_readData(mkw2Filename, 'include', var_import, '~typeOutcomeCheck');
    MW_writeH5(d, [source_dir 'h5/' fname '.h5'], 'replace', 'privateCFG', cfg_pth)
else
    % Import .H5
    d = MW_readData(h5Filename, 'include', var_import, '~typeOutcomeCheck');
end

%% Sync MWorks and Plexon files

if isfile([dest_dir 'syncParam_' fname '.mat'])
    load([dest_dir '/syncParam_' fname]) % Load synchronisation parameters
else
    [sync_gain, sync_offset] = MW_getSyncParam(mkw2Filename, pl2Filename,'precision',4000, 'syncVarName','IO_syncWord');
    save([dest_dir '/syncParam_' fname],'sync_gain','sync_offset','-v7.3')
end

sync_gain = (sync_gain);
sync_offset = double(sync_offset);
%% Process data

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
idx.ttype                   = d.event == 'TRIAL_type';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.trg                     = d.event == 'TRIAL_reactionTrigger';
idx.eye_x_dva              	= d.event == 'EYE_x_dva';
idx.eye_y_dva             	= d.event == 'EYE_y_dva';
idx.reward                  = d.event == 'INFO_Juice_ml';
idx.pump                    = d.event == 'IO_rewardA_ml';
idx.score                   = d.event == 'INFO_score';
idx.task                    = d.event == 'INFO_task';
rew_str                     = 'ml';

% Get trial timestamps
cyc.cOn                     = double(d.time(idx.cOn));
cyc.cEnd                    = double(d.time(idx.cEnd));
ccnt                        = 0;

% Experiment information
exp_info = split(fname,'_');
exp.date = exp_info{1};
exp.monkey = exp_info{2};
exp.block =  exp_info{4};
exp.setup =  exp_info{5};
exp.rec_num = exp_info{6};
exp.experimenter = exp_info{7};

% Stimulus cycle information
for iCyc = 1:length(cyc.cEnd)
    disp(['Processing Cycle: '  num2str(iCyc)])

    % Trial index
    cycIdx                  = [];
    cycIdx                  = d.time >= cyc.cOn(iCyc) & d.time <= cyc.cEnd(iCyc);
    task_str{iCyc}        	= d.value(cycIdx & idx.task);

    if isempty(task_str)
        continue
    end

    if strcmp(task_str{iCyc}, 'CPR_solo_stepfunction_neutral')
        % Cycle counter
        ccnt                    = ccnt+1;

        % Stimulus cycle parameters
        x_pos                   = cell2mat(d.value(cycIdx & d.event == 'STIM_RDP_posX'));
        y_pos                   = cell2mat(d.value(cycIdx & d.event == 'STIM_RDP_posY'));
        cyc.rdp_center_xy{ccnt} = [x_pos y_pos];

        cyc.rdp_dir{ccnt}       = cell2mat(d.value(cycIdx & idx.RDP_dir));
        cyc.rdp_coh{ccnt}       = cell2mat(d.value(cycIdx & idx.RDP_coh));

        cyc.rdp_dir_ts{ccnt}    = double(d.time(cycIdx & idx.RDP_dir));
        cyc.rdp_coh_ts{ccnt}    = double(d.time(cycIdx & idx.RDP_coh));

        cyc.js_dir{ccnt}        = cell2mat(d.value(cycIdx & idx.JS_dir));
        cyc.js_tlt{ccnt}        = cell2mat(d.value(cycIdx & idx.JS_str));
        cyc.js_ts{ccnt}         = double(d.time(cycIdx & idx.JS_str));

        tmp_trg_val             = cell2mat(d.value(cycIdx & idx.trg_on));
        tmp_trg_ts              = double(d.time(cycIdx & idx.trg_on));
        cyc.feedback_ts{ccnt}   = tmp_trg_ts(tmp_trg_val==1);
        cyc.outcome{ccnt}       = d.value(cycIdx & idx.outcome);
        cyc.reward_ind{ccnt}    = cell2mat(d.value(cycIdx & idx.pump));
        cyc.reward_cum{ccnt}    = cell2mat(d.value(cycIdx & idx.reward));
        cyc.score_cum{ccnt}     = cell2mat(d.value(cycIdx & idx.score));
    end
end

%% Import spiking data
spk_files               = dir([dest_dir 'dataspikes*.mat']); % Get all spike files

for iChan = 1%length(spk_files)
    clear spks units
    disp(['Processing : '  spk_files(iChan).name])
    chan_info           = split(spk_files(iChan).name,'_');
    chan_str            = [chan_info{2} '_' chan_info{3}(1:3)];

    % Import (sorted) spike file
    load([spk_files(iChan).folder '/'  spk_files(iChan).name])

    % Correkt timestamps
    spks(:,1)           = cluster_class(:,1); % ID
    spks(:,2)           = double((cluster_class(:,2).*1e3) .* sync_gain) + sync_offset; % timestamp
    spks                = double(spks);

    %%% RF analysis %%%
    %%% 200ms after trial onset without stimulus! %%%
    clear RF
    [brain.RF.(chan_str).nSpikes, brain.RF.stim_id, brain.RF.stim_pos, brain.RF.trl_num,brain.RF.raw.(chan_str)] = RF_mapping(d,idx,spks);

    %%% CPR analysis %%%
    units               = unique(spks(:,1));

    for iUnit = 1:length(units)
        disp(['Processing Unit: '  num2str(iUnit)])
        unit_str        = ['unit' num2str(units(iUnit))];

        ccnt            = 0; % Reset cycle counter
        for iCyc = 1:length(cyc.cEnd)
            if strcmp(task_str{iCyc}, 'CPR_solo_stepfunction_neutral')
                ccnt            = ccnt+1; % stimulus cycle count

                % Cycle-wise spiking for each cluster
                cyc.cpr_cyle(ccnt,:) = [cyc.cOn(iCyc) cyc.cEnd(iCyc)];
                offset          = 300e3; 
                unitIdx         = spks(:,1) == units(iUnit); 
                blIdx           = spks(:,2) >= cyc.cOn(iCyc)-offset & spks(:,2) < cyc.cOn(iCyc);
                cycleIdx        = spks(:,2) >= cyc.cOn(iCyc) & spks(:,2) <= cyc.cEnd(iCyc);
                brain.CPR.spks.cyc.(chan_str).(unit_str){ccnt} = spks(unitIdx & cycleIdx,2);
                brain.CPR.spks.bl.(chan_str).(unit_str){ccnt} = spks(unitIdx & blIdx,2);

                % Other signals - add later
                brain.CPR.lfp   = [];
                brain.CPR.muae  = [];
            end
        end
    end
end

%% Format output
out.exp         = exp;      % Experimental information
out.cyc         = cyc;      % Stimulus cycle information
out.brain       = brain;    % Neural data

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nSpikes, stim_id, stim_pos, trl_num,raw] = RF_mapping(d,idx,spks)

idx.task                    = d.event == 'INFO_task';
idx.tStart                  = d.event == 'TRIAL_start';
idx.tEnd                    = d.event == 'TRIAL_end';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.type                    = d.event == 'TRIAL_type';

trl.start                   = d.time(idx.tStart);
trl.end                     = d.time(idx.tEnd);

tmp_task                    = d.value(idx.task);
tmp_task_ts                 = d.time(idx.task);
for iTrl = 1:length(trl.start)
    trl.task{iTrl}          = tmp_task{find(tmp_task_ts > trl.start(iTrl),1,'first')};
end

units                       = unique(spks(:,1));

x_position                  = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
y_position                  = fliplr([-12 -9 -6 -3 0 3 6 9 12]);
win_offset                  = [50e3 100e3];
stim_id                     = reshape(1:(length(x_position)*length(y_position)),[length(y_position),length(x_position)]);
stim_cnt                    = 0;

% Trial loop
for iTrl = 1:length(trl.start)

    if ~strcmp(trl.task{iTrl}, 'RF_mapping')
        continue
    end

    % Trial index
    trlIdx                      = [];
    trlIdx                      = d.time >= trl.start(iTrl) & d.time <= trl.end(iTrl);
    spks_trlIdx                 = spks(:,2) >= trl.start(iTrl) & spks(:,2) <= trl.end(iTrl);
    raw.spks_id{iTrl}           = spks(spks_trlIdx,1);
    raw.spks_ts{iTrl}           = spks(spks_trlIdx,2) - double(trl.start(iTrl));

    % Trial outcome
    outcome                     = getTrialData(d.value, trlIdx, idx.outcome);
    ttype                       = getTrialData(d.value, trlIdx, idx.type);
    trl_ID                      = getTrialData(d.value, trlIdx, idx.tStart);
    rdp_time                    = getTrialData(d.time, trlIdx, idx.type);

    if ~iscell(ttype)
        continue
    end

    if strcmp(outcome, 'fixation break')
        nStim                   = length(ttype) - 1; % Exclude last stimulus
    else
        nStim                   = length(ttype);
    end

    % Baseline spiking
    blIdx                     	= [];
    blIdx                       = spks(:,2) >= trl.start(iTrl) & spks(:,2) <= rdp_time(1);
    bl_spike_id                 = spks(blIdx,1);
    bl_spike_time               = spks(blIdx,2);

    % Stimulus position loop
    for iStim = 2:nStim

        % Stimulus index - offset by certain lag
        sIdx                    = [];
        %%% last stim missing %%%
        sIdx                    = spks(:,2) >= rdp_time(iStim-1)+win_offset(1) & spks(:,2) <= rdp_time(iStim)+win_offset(2);
        % Stimulus position
        pos                     = split(ttype{iStim-1},'_');
        rdp_x                   = str2num(pos{1}(2:end));
        rdp_y                   = str2num(pos{2}(2:end));
        stim_cnt                = stim_cnt+1;
        stim_pos(stim_cnt)      = stim_id(rdp_y == y_position, rdp_x == x_position);
        trl_num(stim_cnt)       = trl_ID;

        % Spike count
        spike_id                = spks(sIdx,1);
        spike_time              = spks(sIdx,2);

        sid{stim_cnt}           = spike_id;
        sts{stim_cnt}           = spike_time;

        % Unit loop
        for iUnit = 1:length(units)
            bl_spikes          	= length(bl_spike_time(bl_spike_id == units(iUnit)));
            nSpikes(iUnit,stim_cnt) = length(spike_time(spike_id == units(iUnit))) - bl_spikes;
        end
    end
end
end
