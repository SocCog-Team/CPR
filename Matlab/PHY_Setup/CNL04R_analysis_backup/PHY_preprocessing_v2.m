function [out] = PHY_preprocessing(fname, source_dir, dest_dir, cfg_pth, import_flag)
%
% Preprocessing pipeline: imports MWorks and Plexon data, computes
% stimulus-aligned spike counts, and characterises receptive fields.
%
% INPUT
%   fname        String   file name (used to construct all paths)
%   source_dir   String   path to data source directory
%   dest_dir     String   path to file destination directory
%   cfg_pth      String   path to configuration .cfg file
%   import_flag  Logical  true = import from .mwk2, false = load .h5
%
% OUTPUT
%   out   struct   experiment, stimulus, joystick, and neural data
%
% Felix Schneider, CNL
%
% Version history
%   1.0  (fxs 2025-09-05)  Initial version.
%   1.1  (fxs 2026-05-21)  Integrated RF_characterisation.
%
%%% TO DO %%%
% (1) Add LFP/MUAe signal


% =========================================================================
%% Adjust paths
% =========================================================================

addpath(genpath('/Users/cnl/Documents/GitLab/matlab4mworks/'));

mkw2Filename = [source_dir 'mwk2/' fname '.mwk2'];
h5Filename   = [source_dir 'h5/'   fname '.h5'];
pl2Filename  = [source_dir 'pl2/'  fname '.pl2'];


% =========================================================================
%% Import MWorks data
% =========================================================================

close all
clear d

var_import = {
    'INFO_', ...
    'TRIAL_', ...
    'STIM_', ...
    'IO_joystickDirection', ...
    'IO_joystickStrength', ...
    'IO_fixation_flag', ...
    'IO_sync_16bit', ...
    'IO_syncWord', ...
    'IO_rewardA_ml', ...
    'EYE_x_dva', ...
    'EYE_y_dva', ...
    '#stimDisplayUpdate'};

if import_flag
    d = MW_readData(mkw2Filename, 'include', var_import, '~typeOutcomeCheck');
    MW_writeH5(d, h5Filename, 'replace', 'privateCFG', cfg_pth);
else
    d = MW_readData(h5Filename, 'include', var_import, '~typeOutcomeCheck');
end


% =========================================================================
%% Sync MWorks and Plexon files
% =========================================================================

if isfile([dest_dir 'syncParam_' fname '.mat'])
    load([dest_dir '/syncParam_' fname])
else
    [sync_gain, sync_offset] = MW_getSyncParam(mkw2Filename, pl2Filename, ...
        'precision', 4000, 'syncVarName', 'IO_sync_16bit');
    save([dest_dir '/syncParam_' fname], 'sync_gain', 'sync_offset', '-v7.3')
end

sync_gain   = double(sync_gain);
sync_offset = double(sync_offset);


% =========================================================================
%% Process behavioural data
% =========================================================================

idx  = [];
stim = [];

% Event indices
idx.cOn      = d.event == 'TRIAL_start';
idx.cEnd     = d.event == 'TRIAL_end';
idx.frame    = d.event == 'STIM_displayUpdate';
idx.RDP_onset= d.event == 'STIM_RDP_onset';
idx.RDP_dir  = d.event == 'STIM_RDP_direction';
idx.RDP_coh  = d.event == 'STIM_RDP_coherence';
idx.RDP_dot  = d.event == 'STIM_RDP_dotPositions';
idx.trg_on   = d.event == 'STIM_target_onset';
idx.JS_dir   = d.event == 'IO_joystickDirection';
idx.JS_str   = d.event == 'IO_joystickStrength';
idx.JS2_dir  = d.event == 'IO_joystickDirection2';
idx.JS2_str  = d.event == 'IO_joystickStrength2';
idx.fixation = d.event == 'IO_fixation_flag';
idx.ttype    = d.event == 'TRIAL_type';
idx.outcome  = d.event == 'TRIAL_outcome';
idx.trg      = d.event == 'TRIAL_reactionTrigger';
idx.eye_x    = d.event == 'EYE_x_dva';
idx.eye_y    = d.event == 'EYE_y_dva';
idx.reward   = d.event == 'INFO_Juice_ml';
idx.pump     = d.event == 'IO_rewardA';
idx.score    = d.event == 'INFO_score';
idx.task     = d.event == 'INFO_task';

% Trial timestamps
stim.cOn  = double(d.time(idx.cOn));
stim.cEnd = double(d.time(idx.cEnd));
ccnt      = 0;

% Experiment metadata from filename
exp_info         = split(fname, '_');
exp.date         = exp_info{1};
exp.monkey       = exp_info{2};
exp.block        = exp_info{4};
exp.setup        = exp_info{5};
exp.rec_num      = exp_info{6};
exp.experimenter = exp_info{7};

% Stimulus cycle loop — CPR trials only
for iCyc = 1:numel(stim.cEnd)
    disp(['Processing cycle: ' num2str(iCyc)])

    cycIdx         = d.time >= stim.cOn(iCyc) & d.time <= stim.cEnd(iCyc);
    stim.task{iCyc}= d.value(cycIdx & idx.task);

    if isempty(stim.task{iCyc}) || ~contains(stim.task{iCyc}, 'CPR')
        continue
    end

    ccnt = ccnt + 1;

    stim.cpr_solo{ccnt}  = contains(stim.task{iCyc}, 'solo');
    stim.cpr_dyad{ccnt}  = contains(stim.task{iCyc}, 'dyad');
    stim.cpr_catch{ccnt} = contains(stim.task{iCyc}, 'Catch');

    x_pos = cell2mat(d.value(cycIdx & d.event == 'STIM_RDP_posX'));
    y_pos = cell2mat(d.value(cycIdx & d.event == 'STIM_RDP_posY'));
    stim.rdp_center_xy{ccnt} = [x_pos y_pos];

    stim.rdp_dir{ccnt}    = cell2mat(d.value(cycIdx & idx.RDP_dir));
    stim.rdp_coh{ccnt}    = cell2mat(d.value(cycIdx & idx.RDP_coh));
    stim.frme_ts{ccnt}    = double(d.time(cycIdx & idx.frame));
    stim.rdp_dir_ts{ccnt} = double(d.time(cycIdx & idx.RDP_dir));
    stim.rdp_coh_ts{ccnt} = double(d.time(cycIdx & idx.RDP_coh));

    tmp_trg_val             = cell2mat(d.value(cycIdx & idx.trg_on));
    tmp_trg_ts              = double(d.time(cycIdx & idx.trg_on));
    stim.feedback_ts{ccnt}  = tmp_trg_ts(tmp_trg_val == 1);
    stim.outcome{ccnt}      = d.value(cycIdx & idx.outcome);
    stim.reward_ind{ccnt}   = cell2mat(d.value(cycIdx & idx.pump));
    stim.reward_cum{ccnt}   = cellfun(@double, d.value(cycIdx & idx.reward));
    stim.score_cum{ccnt}    = cell2mat(d.value(cycIdx & idx.score));

    joy.js_monk_dir{ccnt} = cell2mat(d.value(cycIdx & idx.JS_dir));
    joy.js_monk_tlt{ccnt} = cell2mat(d.value(cycIdx & idx.JS_str));
    joy.js_monk_ts{ccnt}  = double(d.time(cycIdx & idx.JS_str));
    joy.js_hum_dir{ccnt}  = cell2mat(d.value(cycIdx & idx.JS2_dir));
    joy.js_hum_tlt{ccnt}  = cell2mat(d.value(cycIdx & idx.JS2_str));
    joy.js_hum_ts{ccnt}   = double(d.time(cycIdx & idx.JS2_str));
end


% =========================================================================
%% Import spiking data and run RF characterisation
% =========================================================================

spk_files = dir([dest_dir 'dataspikes*negthr.mat']);

for iChan = 1:numel(spk_files)

    clear spks units
    disp(['Processing: ' spk_files(iChan).name])

    % Build the channel string (e.g. 'ch001_neg') from the filename
    chan_info = split(spk_files(iChan).name, '_');
    chan_str  = [chan_info{2} '_' chan_info{3}(1:3)];

    % Load sorted spike data and apply synchronisation correction
    load([spk_files(iChan).folder '/' spk_files(iChan).name])
    spks(:,1) = cluster_class(:,1);                                        % Phy unit ID
    spks(:,2) = double((cluster_class(:,2).*1e3) .* sync_gain) + sync_offset; % µs timestamp
    spks      = double(spks);

    % -----------------------------------------------------------------
    % RF mapping — build the baseline-corrected spike-count map for all
    % units on this channel in one pass.  stim_id and stim_pos are shared
    % across channels and are written on the first channel; subsequent
    % channels overwrite with identical values.
    % -----------------------------------------------------------------
    [brain.RF.(chan_str).nSpikes, ...
     brain.RF.stim_id, ...
     brain.RF.stim_pos, ...
     brain.RF.trl_num, ...
     brain.RF.raw.(chan_str)] = RF_mapping(d, idx, spks);

    % -----------------------------------------------------------------
    % RF characterisation — run per unit.
    % RF_characterisation combines:
    %   Part A: fixation position [0,0] onset response vs baseline
    %           (Wilcoxon signed-rank test)
    %   Part B: 2-D Gaussian RF fit, area, distance to fixation, and
    %           cursor overlap (2-sigma ellipse ∩ cursor aperture)
    %
    % Results are stored under brain.RF.(chan_str).(unit_str) where
    % unit_str is the zero-padded Phy unit ID, e.g. 'unit002'.
    % -----------------------------------------------------------------
    units = unique(spks(:,1));

    for iUnit = 1:numel(units)
        unit_str = sprintf('unit%03d', units(iUnit));
        disp(['  RF characterisation — unit ' unit_str])

        brain.RF.(chan_str).(unit_str) = RF_characterisation( ...
            d, idx, spks, ...
            brain.RF.stim_id, ...
            brain.RF.stim_pos, ...
            brain.RF.(chan_str).nSpikes, ...
            iUnit);
    end

    % -----------------------------------------------------------------
    % CPR analysis — collect spike trains aligned to cycle onset for
    % every unit, including a 300 ms pre-trial baseline period.
    % -----------------------------------------------------------------
    ccnt = 0;

    for iUnit = 1:numel(units)
        disp(['  CPR — unit ' num2str(iUnit)])
        unit_str = ['unit' num2str(units(iUnit))];
        ccnt     = 0;

        for iCyc = 1:numel(stim.cEnd)
            if contains(stim.task{iCyc}, 'CPR') || contains(stim.task{iCyc}, 'Catch')
                ccnt = ccnt + 1;

                stim.cpr_cyle(ccnt,:) = [stim.cOn(iCyc) stim.cEnd(iCyc)];
                offset    = 300e3;   % µs — pre-stimulus baseline window
                unitIdx   = spks(:,1) == units(iUnit);
                cycleIdx  = spks(:,2) >= stim.cOn(iCyc) - offset & ...
                             spks(:,2) <= stim.cEnd(iCyc);

                % Spike timestamps relative to cycle onset
                brain.CPR.spks.(chan_str).(unit_str){ccnt} = ...
                    spks(unitIdx & cycleIdx, 2) - stim.cOn(iCyc);
            end
        end
    end

end % iChan


% =========================================================================
%% Format output
% =========================================================================

out.exp   = exp;
out.stim  = stim;
out.joy   = joy;
out.brain = brain;

end % PHY_preprocessing


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nSpikes, stim_id, stim_pos, trl_num, raw] = RF_mapping(d, idx, spks)
% RF_MAPPING  Build a baseline-corrected spike-count map across the RF grid.
%
% For each RF_mapping trial, spikes are counted in a response window
% (RDP onset +50 ms to +150 ms) and baseline-corrected by subtracting the
% spike count in the pre-stimulus window (trial start → first RDP onset).
%
% All timestamps are in MICROSECONDS.

idx.task    = d.event == 'INFO_task';
idx.tStart  = d.event == 'TRIAL_start';
idx.tEnd    = d.event == 'TRIAL_end';
idx.outcome = d.event == 'TRIAL_outcome';
idx.type    = d.event == 'TRIAL_type';

trl.start = d.time(idx.tStart);
trl.end   = d.time(idx.tEnd);

tmp_task    = d.value(idx.task);
tmp_task_ts = d.time(idx.task);
for iTrl = 1:numel(trl.start)
    trl.task{iTrl} = tmp_task{ find(tmp_task_ts > trl.start(iTrl), 1, 'first') };
end

units      = unique(spks(:,1));
x_position = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
y_position = fliplr([-12 -9 -6 -3 0 3 6 9 12]);
win_offset = [50e3 150e3];   % µs

stim_id  = reshape(1:(numel(x_position)*numel(y_position)), ...
                   [numel(y_position), numel(x_position)]);
stim_cnt = 0;

for iTrl = 1:numel(trl.start)

    if ~strcmp(trl.task{iTrl}, 'RF_mapping')
        continue
    end

    trlIdx      = d.time >= trl.start(iTrl) & d.time <= trl.end(iTrl);
    spks_trlIdx = spks(:,2) >= trl.start(iTrl) & spks(:,2) <= trl.end(iTrl);

    % Store raw spike times aligned to trial onset for this trial
    raw.spks_id{iTrl} = spks(spks_trlIdx, 1);
    raw.spks_ts{iTrl} = spks(spks_trlIdx, 2) - double(trl.start(iTrl));

    outcome  = getTrialData(d.value, trlIdx, idx.outcome);
    ttype    = getTrialData(d.value, trlIdx, idx.type);
    trl_ID   = getTrialData(d.value, trlIdx, idx.tStart);
    rdp_time = getTrialData(d.time,  trlIdx, idx.type);

    if ~iscell(ttype)
        continue
    end

    if strcmp(outcome, 'fixation break')
        nStim = numel(ttype) - 1;
    else
        nStim = numel(ttype);
    end

    % Baseline: trial start → first RDP onset (no stimulus on screen)
    blIdx        = spks(:,2) >= trl.start(iTrl) & spks(:,2) <= rdp_time(1);
    bl_spike_id  = spks(blIdx, 1);
    bl_spike_time= spks(blIdx, 2);

    for iStim = 1:nStim

        % Response window anchored to this stimulus onset
        sIdx       = spks(:,2) >= rdp_time(iStim) + win_offset(1) & ...
                     spks(:,2) <= rdp_time(iStim) + win_offset(2);
        spike_id   = spks(sIdx, 1);
        spike_time = spks(sIdx, 2);

        % Parse stimulus position from the ttype string (e.g. 'X-6_Y3')
        pos   = split(ttype{iStim}, '_');
        rdp_x = str2double(pos{1}(2:end));
        rdp_y = str2double(pos{2}(2:end));

        stim_cnt            = stim_cnt + 1;
        stim_pos(stim_cnt)  = stim_id(rdp_y == y_position, rdp_x == x_position);
        trl_num(stim_cnt)   = trl_ID;

        % Baseline-corrected spike count per unit
        for iUnit = 1:numel(units)
            bl_spikes              = sum(bl_spike_id == units(iUnit));
            nSpikes(iUnit,stim_cnt)= sum(spike_id    == units(iUnit)) - bl_spikes;
        end
    end
end

end % RF_mapping


function out = RF_characterisation(d, idx, spks, stim_id, stim_pos, nSpikes, iUnit)
% RF_CHARACTERISATION  Characterise the receptive field of a single unit
%                      and assess its functional relationship to fixation.
%
% Integrates two complementary analyses into one preprocessing step:
%
%   PART A — Fixation position response  (position [0,0])
%     Spike counts are extracted from a stimulus-free baseline window
%     (trial start → first RDP onset) and a response window (stim onset
%     +50 ms to +150 ms) for every [0,0] presentation.  A Wilcoxon
%     signed-rank test assesses whether the response differs significantly
%     from baseline across presentations.
%
%   PART B — Receptive field mapping
%     A baseline-corrected spike-count map is built from nSpikes.  A 2-D
%     Gaussian is fitted to estimate RF centre and size.  Spatial overlap
%     between the 2-sigma RF ellipse (~86% of Gaussian volume) and the
%     cursor aperture is quantified numerically on a 0.05 dva grid.
%
% All timestamps in spks are in MICROSECONDS.
%
% INPUT
%   d         struct          event/time/value data struct
%   idx       struct          pre-built event index struct
%   spks      [N×2]           col 1 = Phy unit ID, col 2 = timestamp µs
%   stim_id   [nRows×nCols]   grid position → stimulus ID mapping
%   stim_pos  [1×nTrials]     stimulus ID shown on each trial
%   nSpikes   [nUnits×nTrials] BL-corrected spike counts from RF_mapping
%   iUnit     scalar integer  unit index — row into nSpikes, also used to
%                             look up the Phy ID from unique(spks(:,1))
%
% OUTPUT
%   out   struct
%           out.params   — all analysis parameters
%           out.fixpos   — Part A results
%           out.rf       — Part B results


% =========================================================================
% 1.  SHARED PARAMETERS
% =========================================================================

x_dva        = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
y_dva        = fliplr([-12 -9 -6 -3 0 3 6 9 12]);
grid_spacing = 3;
n_stim       = stim_id(end, end);

fix_dva       = [0, 0];
fixwin_radius = 1.50;   % dva
cursor_radius = 1.75;   % dva
cursor_area   = pi * cursor_radius^2;

% iUnit indexes the sorted list of Phy unit IDs on this channel.
% This must match the ordering used in RF_mapping (both use unique(spks(:,1))).
all_unit_ids = unique(spks(:,1));
unit_id      = all_unit_ids(iUnit);

% Part A
target_pos = 'X0_Y0';
win_offset = [50e3, 150e3];   % µs — must match RF_mapping win_offset
alpha      = 0.05;

% Part B
min_spike_threshold = 5;     % peak BL-corrected count below this → skip fit
overlap_grid_step   = 0.05;  % dva
rf_sd_boundary      = 2;     % number of SDs defining the RF boundary
rf_boundary_thr     = rf_sd_boundary^2;   % threshold on normalised squared distance


% =========================================================================
% 2.  INITIALISE OUTPUT WITH SAFE DEFAULTS
% =========================================================================

out.unit_idx = iUnit;
out.unit_id  = unit_id;

out.fixpos.n_presentations  = 0;
out.fixpos.baseline_counts  = [];
out.fixpos.response_counts  = [];
out.fixpos.baseline_Hz      = [];
out.fixpos.response_Hz      = [];
out.fixpos.mean_baseline_Hz = NaN;
out.fixpos.sem_baseline_Hz  = NaN;
out.fixpos.mean_response_Hz = NaN;
out.fixpos.sem_response_Hz  = NaN;
out.fixpos.delta_Hz         = NaN;
out.fixpos.wilcoxon_W       = NaN;
out.fixpos.p_value          = NaN;
out.fixpos.significant      = false;
out.fixpos.trial_numbers    = [];
out.fixpos.stim_onset_us    = [];
out.fixpos.outcomes         = {};
out.fixpos.baseline_dur_us  = [];
out.fixpos.response_dur_us  = [];

out.rf.responsive            = false;
out.rf.fit_ok                = false;
out.rf.spike_map             = [];
out.rf.RF_x_dva              = NaN;
out.rf.RF_y_dva              = NaN;
out.rf.RF_sigma_x_dva        = NaN;
out.rf.RF_sigma_y_dva        = NaN;
out.rf.RF_area_dva2          = NaN;
out.rf.dist_to_fix_dva       = NaN;
out.rf.mean_act_in_cursor    = NaN;
out.rf.mean_spk_inside_RF    = NaN;
out.rf.mean_spk_outside_RF   = NaN;
out.rf.overlap_flag          = false;
out.rf.overlap_area_dva2     = NaN;
out.rf.overlap_pct_of_RF     = NaN;
out.rf.overlap_pct_of_cursor = NaN;

out.params.unit_id           = unit_id;
out.params.iUnit             = iUnit;
out.params.x_dva             = x_dva;
out.params.y_dva             = y_dva;
out.params.grid_spacing      = grid_spacing;
out.params.fix_dva           = fix_dva;
out.params.fixwin_radius     = fixwin_radius;
out.params.cursor_radius     = cursor_radius;
out.params.target_pos        = target_pos;
out.params.win_offset_us     = win_offset;
out.params.alpha             = alpha;
out.params.min_spike_thr     = min_spike_threshold;
out.params.overlap_grid_step = overlap_grid_step;
out.params.rf_sd_boundary    = rf_sd_boundary;


% =========================================================================
% 3.  EVENT INDICES
% =========================================================================

idx.task    = d.event == 'INFO_task';
idx.tStart  = d.event == 'TRIAL_start';
idx.tEnd    = d.event == 'TRIAL_end';
idx.outcome = d.event == 'TRIAL_outcome';
idx.type    = d.event == 'TRIAL_type';

trl.start = d.time(idx.tStart);
trl.end   = d.time(idx.tEnd);

tmp_task    = d.value(idx.task);
tmp_task_ts = d.time(idx.task);
for iTrl = 1:numel(trl.start)
    trl.task{iTrl} = tmp_task{ find(tmp_task_ts > trl.start(iTrl), 1, 'first') };
end


% =========================================================================
% PART A — FIXATION POSITION RESPONSE
% =========================================================================

max_pres         = 5000;
bl_counts_raw    = nan(1, max_pres);
resp_counts_raw  = nan(1, max_pres);
pres_trial_num   = nan(1, max_pres);
pres_stim_onset  = nan(1, max_pres);
pres_outcome     = cell(1, max_pres);
pres_bl_dur_us   = nan(1, max_pres);
pres_resp_dur_us = nan(1, max_pres);
pres_cnt         = 0;

for iTrl = 1:numel(trl.start)

    if ~strcmp(trl.task{iTrl}, 'RF_mapping')
        continue
    end

    trlIdx   = d.time >= trl.start(iTrl) & d.time <= trl.end(iTrl);
    outcome  = getTrialData(d.value, trlIdx, idx.outcome);
    ttype    = getTrialData(d.value, trlIdx, idx.type);
    trl_ID   = getTrialData(d.value, trlIdx, idx.tStart);
    rdp_time = getTrialData(d.time,  trlIdx, idx.type);

    if ~iscell(ttype)
        continue
    end

    if strcmp(outcome, 'fixation break')
        n_stim_trl = numel(ttype) - 1;
    else
        n_stim_trl = numel(ttype);
    end

    % Baseline: trial start → first RDP onset (stimulus-free period)
    bl_t_start = trl.start(iTrl);
    bl_t_end   = rdp_time(1);
    bl_dur_us  = bl_t_end - bl_t_start;

    bl_mask    = spks(:,2) >= bl_t_start & spks(:,2) < bl_t_end;
    bl_spk_ids = spks(bl_mask, 1);

    for iStim = 1:n_stim_trl

        if ~strcmp(ttype{iStim}, target_pos)
            continue
        end

        resp_t_start = rdp_time(iStim) + win_offset(1);
        resp_t_end   = rdp_time(iStim) + win_offset(2);
        resp_dur_us  = resp_t_end - resp_t_start;

        resp_mask    = spks(:,2) >= resp_t_start & spks(:,2) <= resp_t_end;
        resp_spk_ids = spks(resp_mask, 1);

        pres_cnt = pres_cnt + 1;

        bl_counts_raw(pres_cnt)    = sum(bl_spk_ids  == unit_id);
        resp_counts_raw(pres_cnt)  = sum(resp_spk_ids == unit_id);
        pres_trial_num(pres_cnt)   = trl_ID;
        pres_stim_onset(pres_cnt)  = rdp_time(iStim);
        pres_outcome{pres_cnt}     = outcome;
        pres_bl_dur_us(pres_cnt)   = bl_dur_us;
        pres_resp_dur_us(pres_cnt) = resp_dur_us;

    end
end

if pres_cnt > 0

    bl_counts   = double(bl_counts_raw(1:pres_cnt));
    resp_counts = double(resp_counts_raw(1:pres_cnt));
    bl_dur_s    = double(pres_bl_dur_us(1:pres_cnt))   / 1e6;
    resp_dur_s  = double(pres_resp_dur_us(1:pres_cnt)) / 1e6;

    baseline_Hz = bl_counts   ./ bl_dur_s;
    response_Hz = resp_counts ./ resp_dur_s;

    if any(bl_counts ~= resp_counts)
        [pval, hsig, wstats] = signrank(bl_counts', resp_counts', 'alpha', alpha);
        wstat = wstats.signedrank;
    else
        pval  = 1;
        hsig  = false;
        wstat = 0;
    end

    out.fixpos.n_presentations  = pres_cnt;
    out.fixpos.baseline_counts  = bl_counts;
    out.fixpos.response_counts  = resp_counts;
    out.fixpos.baseline_Hz      = baseline_Hz;
    out.fixpos.response_Hz      = response_Hz;
    out.fixpos.mean_baseline_Hz = mean(baseline_Hz);
    out.fixpos.sem_baseline_Hz  = std(baseline_Hz) / sqrt(pres_cnt);
    out.fixpos.mean_response_Hz = mean(response_Hz);
    out.fixpos.sem_response_Hz  = std(response_Hz) / sqrt(pres_cnt);
    out.fixpos.delta_Hz         = mean(response_Hz) - mean(baseline_Hz);
    out.fixpos.wilcoxon_W       = wstat;
    out.fixpos.p_value          = pval;
    out.fixpos.significant      = logical(hsig);
    out.fixpos.trial_numbers    = pres_trial_num(1:pres_cnt);
    out.fixpos.stim_onset_us    = pres_stim_onset(1:pres_cnt);
    out.fixpos.outcomes         = pres_outcome(1:pres_cnt);
    out.fixpos.baseline_dur_us  = pres_bl_dur_us(1:pres_cnt);
    out.fixpos.response_dur_us  = pres_resp_dur_us(1:pres_cnt);

end


% =========================================================================
% PART B — RECEPTIVE FIELD MAPPING AND CURSOR OVERLAP
% =========================================================================

[X_dva_grid, Y_dva_grid] = meshgrid(x_dva, y_dva);
dist_from_fix_dva         = sqrt(X_dva_grid.^2 + Y_dva_grid.^2);
cursor_mask               = dist_from_fix_dva <= cursor_radius;
xy_flat                   = [X_dva_grid(:), Y_dva_grid(:)];

% Build spike map from the already-computed nSpikes matrix.
% nSpikes rows correspond to unique(spks(:,1)) in ascending order, which
% is the same ordering used here, so iUnit indexes the correct row.
spike_map = nan(size(stim_id));
for iPos = 1:n_stim
    spike_map(stim_id == iPos) = sum(nSpikes(iUnit, stim_pos == iPos));
end
out.rf.spike_map          = spike_map;
out.rf.mean_act_in_cursor = mean(spike_map(cursor_mask), 'omitnan');

% Quality gate — skip units with a flat response map
if max(spike_map(:), [], 'omitnan') < min_spike_threshold
    return
end
out.rf.responsive = true;

% 2-D Gaussian fit
% p = [amplitude, x0, y0, sigma_x, sigma_y, offset]
gauss2d_fn = @(p, xy) ...
    p(1) .* exp( -( (xy(:,1)-p(2)).^2 ./ (2.*p(4).^2) + ...
                    (xy(:,2)-p(3)).^2 ./ (2.*p(5).^2) ) ) + p(6);

p0  = [5,   0,   0,   6,   6,   0];
lb  = [0, -24, -12, 0.5, 0.5, -Inf];
ub  = [Inf,  3,  12,  20,  20,  Inf];
opt = optimoptions('lsqcurvefit', 'Display', 'off');

z_flat = spike_map(:);
valid  = ~isnan(z_flat);

if sum(valid) >= 6
    try
        p_fit = lsqcurvefit(gauss2d_fn, p0, xy_flat(valid,:), z_flat(valid), lb, ub, opt);
        out.rf.RF_x_dva       = p_fit(2);
        out.rf.RF_y_dva       = p_fit(3);
        out.rf.RF_sigma_x_dva = p_fit(4);
        out.rf.RF_sigma_y_dva = p_fit(5);
        out.rf.fit_ok         = true;
    catch ME
        warning('RF_characterisation: unit %d fit failed (%s).', iUnit, ME.message);
    end
end

if out.rf.fit_ok
    x0 = out.rf.RF_x_dva;
    y0 = out.rf.RF_y_dva;
    sx = out.rf.RF_sigma_x_dva;
    sy = out.rf.RF_sigma_y_dva;

    % RF area at the 2-sigma boundary: π·(2σx)·(2σy)
    out.rf.RF_area_dva2    = pi * (rf_sd_boundary*sx) * (rf_sd_boundary*sy);
    out.rf.dist_to_fix_dva = sqrt(x0^2 + y0^2);

    % Mean spike count inside vs outside the 2-sigma RF ellipse
    ellipse_dist           = ((X_dva_grid-x0)./sx).^2 + ((Y_dva_grid-y0)./sy).^2;
    inside_rf              = ellipse_dist <= rf_boundary_thr;
    out.rf.mean_spk_inside_RF  = mean(spike_map(inside_rf),  'omitnan');
    out.rf.mean_spk_outside_RF = mean(spike_map(~inside_rf), 'omitnan');

    % Numerical overlap: 2-sigma RF ellipse ∩ cursor aperture (0.05 dva grid)
    x_lo = min(x0 - rf_sd_boundary*sx, -cursor_radius) - overlap_grid_step;
    x_hi = max(x0 + rf_sd_boundary*sx,  cursor_radius) + overlap_grid_step;
    y_lo = min(y0 - rf_sd_boundary*sy, -cursor_radius) - overlap_grid_step;
    y_hi = max(y0 + rf_sd_boundary*sy,  cursor_radius) + overlap_grid_step;

    [Xg, Yg]  = meshgrid(x_lo:overlap_grid_step:x_hi, y_lo:overlap_grid_step:y_hi);
    cell_area = overlap_grid_step^2;

    in_rf     = ((Xg-x0)./sx).^2 + ((Yg-y0)./sy).^2 <= rf_boundary_thr;
    in_cursor = Xg.^2 + Yg.^2 <= cursor_radius^2;
    in_both   = in_rf & in_cursor;

    out.rf.overlap_flag          = any(in_both(:));
    out.rf.overlap_area_dva2     = sum(in_both(:)) * cell_area;
    out.rf.overlap_pct_of_RF     = 100 * out.rf.overlap_area_dva2 / out.rf.RF_area_dva2;
    out.rf.overlap_pct_of_cursor = 100 * out.rf.overlap_area_dva2 / cursor_area;
end

end % RF_characterisation