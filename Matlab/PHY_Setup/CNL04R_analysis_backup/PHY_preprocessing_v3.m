function [out] = PHY_preprocessing_v3(fname, source_dir, dest_dir, cfg_pth, import_flag, signal_source, mua_opts)
%
% Preprocessing pipeline: imports MWorks and Plexon data, computes
% stimulus-aligned spike counts, and characterises receptive fields.
% Optionally also derives multi-unit activity (MUAe + MUAt) per channel.
%
% INPUT
%   fname         String   file name (used to construct all paths)
%   source_dir    String   path to data source directory
%   dest_dir      String   path to file destination directory
%   cfg_pth       String   path to configuration .cfg file (h5 variable list)
%   import_flag   Logical  true = import from .mwk2 (and rewrite the .h5),
%                          false = load the .h5. The .h5 keeps only variables
%                          listed in MW_readFile.cfg / cfg_pth and stores floats
%                          as single; h5 files written before 2026-09-11 lack
%                          the CTRL_ variables -> re-import them once (true).
%   signal_source (opt)    cellstr/char, subset of {'sorted','muat','muae'};
%                          which neural signals to produce. Default {'sorted'}.
%                          {'behaviour'} = no neural signal: no sync, no PL2.
%   mua_opts      (opt)    cell of name-value pairs forwarded to fn_compute_MUA,
%                          e.g. {'thr_factor',3.5,'env_lp_Hz',250}. Default {}.
%
% OUTPUT
%   out   struct   experiment, stimulus, joystick, and neural data
%                    .stim.cpr_cyle [nCyc x 2]            CPR cycle onset/end (µs)
%                    .stim.catch_trl [nCatch x 2]         catch trials with cursor: start/end (µs)
%                    .stim.catch_cursor/_ncur/_type       cursor onset (µs), arcs shown, TRIAL_type
%                    .brain.CPR.spks.(chXXX_neg).(unitN)  sorted units
%                    .brain.CPR.spks.(chXXX_mua).muat0    MUAt hash   (if 'muat')
%                    .brain.CPR.muae.(chXXX_mua)          MUAe env    (if 'muae'):
%                       .env / .t_us / .baseline          per CPR cycle; baseline = mean env
%                                                         [onset-300, onset) ms (raw, mV)
%                       .catch_env/_t_us/_baseline        the same per catch trial (cursor onset)
%                       .baseline_ref                     pooled running reference per cycle
%                       .catch_baseline_ref, .ref_k       ... and per catch trial (fn_muae_reference)
%
% Requires for MUA:  fn_compute_MUA.m, fn_muae_reference.m, *_ch*_wb.mat (broadband;
%                    or datafilt_ch*.mat), syncParam_<fname>.mat.  Caches per
%                    channel -> mua_ch###.mat.
%
% Felix Schneider, CNL
%
% Version history
%   1.0  (fxs 2025-09-05)  Initial version.
%   1.1  (fxs 2026-05-21)  Integrated RF_characterisation.
%   2.0  (fxs 2026-07-03)  Added MUAe/MUAt multi-unit path (Stark & Abeles 2007).
%   2.1  (fxs 2026-07-03)  MUA from broadband (*_wb) by default + per-channel cache.
%   2.2  (fxs 2026-09-10)  Partner outcome/reward, reward scale per cycle,
%                          behaviour-only mode ({'behaviour'}), mwk2 name
%                          fallback, time-based per-target outcome / juice /
%                          arc overlap (target_events; replaces index-paired
%                          INFO_Juice_ml, kept as *_counter).
%   2.3  (fxs 2026-09-11)  Speed: per-variable event series + binary search in
%                          the cycle loop (no whole-table masks per cycle),
%                          MUAe cycle slices by binary search, RF-mapping
%                          trials parsed once per session (~25x faster per
%                          MUAe channel; output identical to 2.2). cpr_cyle
%                          built once in the cycle loop. Joystick direction
%                          taken at the strength samples' indices (a pair split
%                          at a cycle edge misaligned direction by one sample;
%                          rec077 crashed in sort_states). IO_rewardA also
%                          matched under its h5 name IO_rewardA_ml; a failed
%                          .h5 write no longer stops the run.
%   2.4  (fxs 2026-09-15)  Catch trials (INFO_task 'Catch_trial'): stim.catch_*
%                          table (cursor onset = first STIM_catcharc(2)_onset
%                          == 1 in the trial; trials without cursor skipped)
%                          and, per MUAe channel, the catch envelope with its
%                          pre-cursor baseline window. Pooled running MUAe
%                          reference per CPR cycle and catch trial
%                          (fn_muae_reference: mean of the 10 nearest valid
%                          baseline windows before + 10 after, CPR and catch
%                          pooled, own window left out) -> baseline_ref,
%                          catch_baseline_ref, ref_k. Existing fields unchanged.


% =========================================================================
%% Optional-argument defaults
% =========================================================================

if nargin < 6 || isempty(signal_source); signal_source = {'sorted'}; end
if nargin < 7 || isempty(mua_opts);      mua_opts      = {};         end
if ischar(signal_source) || isstring(signal_source)
    signal_source = cellstr(signal_source);
end
want_sorted = any(strcmpi(signal_source, 'sorted'));
want_muat   = any(strcmpi(signal_source, 'muat'));
want_muae   = any(strcmpi(signal_source, 'muae'));
% signal_source = {'behaviour'} (no neural signal) -> BEHAVIOUR-ONLY: the
% MWorks-Plexon sync and every PL2 read are skipped, so sessions with an
% unusable PL2 still yield the behavioural summary (exp/stim/joy; brain empty).
want_neural = want_sorted || want_muat || want_muae;


% =========================================================================
%% Adjust paths
% =========================================================================

% matlab4mworks: the driver (PHY_main_analysis_v4) puts it on the path; add the
% laptop copy only when it is missing (genpath on every call is slow).
if ~exist('MW_readData', 'file')
    addpath(genpath('/Users/fschneider/Documents/GitLab/matlab4mworks/'));
end

mkw2Filename = [source_dir 'mwk2/' fname '.mwk2'];
h5Filename   = [source_dir 'h5/'   fname '.h5'];
pl2Filename  = [source_dir 'pl2/'  fname '.pl2'];
if ~isfile(mkw2Filename)
    % Tolerate stray characters before the extension (rec081 was saved as
    % '..._rec081_ann .mwk2') instead of renaming raw data.
    cand = dir([source_dir 'mwk2/' fname '*.mwk2']);
    if numel(cand) == 1
        warning('PHY_preprocessing_v3:mwk2name', 'Using MWorks file "%s".', cand.name);
        mkw2Filename = fullfile(cand.folder, cand.name);
    end
end
if ~isfolder(dest_dir); mkdir(dest_dir); end   % ensure the session output folder exists


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
    'IO_rewardA', ...
    'AGNT_', ...
    'CTRL_reward_target_', ...   % reward scale min/max (+ per-hit target volume)
    'CTRL_reward2_target_ml', ...% partner's reward-equivalent computed at its hit
    'CTRL_reward_power', ...     % hit reward = min + (max-min)*strength^power
    'CTRL_arc', ...              % CTRL_arc_flag / CTRL_arc2_flag: target inside the arc (hit rule)
    'EYE_x_dva', ...
    'EYE_y_dva', ...
    '#stimDisplayUpdate'};

% The .h5 is a cache for later import_flag=false runs (seconds instead of
% ~20 min). MW_writeH5 writes only variables defined in MW_readFile.cfg or
% cfg_pth (floats as single) and renames IO_rewardA -> IO_rewardA_ml, so every
% variable used below must be listed in cfg_pth (CTRL_/AGNT_arc added
% 2026-09-11). A failed write must not stop the run: d comes from the .mwk2.
if import_flag
    d = MW_readData(mkw2Filename, 'include', var_import, '~typeOutcomeCheck');
    try
        MW_writeH5(d, h5Filename, 'replace', 'privateCFG', cfg_pth);
    catch ME
        warning('PHY_preprocessing_v3:h5write', ...
            'Writing %s failed (%s) - continuing without the .h5 cache.', h5Filename, ME.message);
        if isfile(h5Filename); movefile(h5Filename, [h5Filename '.failed']); end   % never leave a partial cache
    end
else
    d = MW_readData(h5Filename, 'include', var_import, '~typeOutcomeCheck');
end


% =========================================================================
%% Sync MWorks and Plexon files
% =========================================================================

if ~want_neural
    sync_gain = NaN;   sync_offset = NaN;          % behaviour-only: no PL2 needed
elseif isfile([dest_dir 'syncParam_' fname '.mat'])
    load([dest_dir 'syncParam_' fname '.mat'], 'sync_gain', 'sync_offset')
else
    [sync_gain, sync_offset] = MW_getSyncParam(mkw2Filename, pl2Filename, ...
        'precision', 4000, 'syncVarName', 'IO_sync_16bit');
    save([dest_dir 'syncParam_' fname '.mat'], 'sync_gain', 'sync_offset', '-v7.3')
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
idx.AGNT_dir = d.event == 'AGNT_direction';    % computer-agent partner: direction (arc centre)
idx.AGNT_str = d.event == 'AGNT_strength';      %                         strength 0-1 (=1-arc_width/180)
idx.fixation = d.event == 'IO_fixation_flag';
idx.ttype    = d.event == 'TRIAL_type';
idx.outcome  = d.event == 'TRIAL_outcome';
idx.trg      = d.event == 'TRIAL_reactionTrigger';
idx.eye_x    = d.event == 'EYE_x_dva';
idx.eye_y    = d.event == 'EYE_y_dva';
idx.reward   = d.event == 'INFO_Juice_ml';
idx.pump     = d.event == 'IO_rewardA' | ...     % reward-board command (mL): raw name (.mwk2)
               d.event == 'IO_rewardA_ml';       % ... renamed by MW_writeH5 (.h5)
idx.score    = d.event == 'INFO_score';
idx.task     = d.event == 'INFO_task';
% --- Player 2 (human / computer-agent partner) + reward scale -------------
% TRIAL_outcome2 / CTRL_reward2_target_ml give the partner's per-target outcome
% and reward-EQUIVALENT (mL; the agent/human receive no juice) — extracted per
% target by time in the EV block below (INFO_Juice2_ml is kept only as a raw log:
% it does not log one entry per target in every task version).
% CTRL_reward_target_min/max_ml are the reward scale; they were changed between
% experiments, so they are stored per cycle to allow normalising reward to
% [min max] downstream (and to detect any within-session change).
idx.outcome2 = d.event == 'TRIAL_outcome2';
idx.reward2  = d.event == 'INFO_Juice2_ml';
idx.rew_min  = d.event == 'CTRL_reward_target_min_ml';
idx.rew_max  = d.event == 'CTRL_reward_target_max_ml';

% --- Per-target outcome, juice and cursor-target overlap: TIME-BASED ----------
% (2026-09-10) The cumulative counters INFO_Juice_ml / INFO_Juice2_ml do NOT log
% one entry per target in every task version: PHY_CPR_20260601 'a' (rec088) logs
% two per dyadic target, PHY_CPR_20260601 'b' (rec089-094) none for a player's
% miss when only the other player hit. Pairing the k-th counter value with the
% k-th outcome (the old reward_ind) therefore misassigned 13-30% of the
% per-target rewards in rec088-094 (session_audit.csv, target_truth). Validated
% on every target of 19 sessions against the task's own 'TARGET REWARD [ml]'
% reports in the raw .mwk2, the per-target truth is:
%   onset   = INFO_TargetCounter increment (state 'Target presentation'); the
%             displayed onset (STIM_target_onset) follows 8-12 ms later
%   outcome = first TRIAL_outcome / TRIAL_outcome2 'hit'|'miss' after the onset
%   juice   = IO_rewardA pulses (the reward-board command, mL) after a monkey
%             hit and before the next target whose value equals the
%             CTRL_reward_target_ml computed at that hit. 2 pulses = the task
%             re-sends the reward when a fixation break follows the hit (trial
%             aborted with CTRL_hit_flag still set); a hit whose command was 0
%             gets 0. Cycle-bonus / catch / fixation pulses never match.
%   partner = CTRL_reward2_target_ml at the partner's hit (reward-equivalent;
%             no juice is given to the partner)
%   overlap = the task's own arc flag (CTRL_arc_flag; partner CTRL_arc2_flag,
%             AGNT_arc_flag for the agent) at the onset or turning 1 within the
%             50 ms target window. Solo: overlap == hit. Dyad: whoever covers the
%             target FIRST ends it, so a 'miss' can hide a later overlap.
idx.tcount   = d.event == 'INFO_TargetCounter';
idx.crew     = d.event == 'CTRL_reward_target_ml';
idx.crew2    = d.event == 'CTRL_reward2_target_ml';
idx.arcflag  = d.event == 'CTRL_arc_flag';
idx.arcflag2 = d.event == 'CTRL_arc2_flag';
idx.agntarc  = d.event == 'AGNT_arc_flag';
EV       = struct();
EV.tc    = ev_series(d, idx.tcount,  'num');
EV.out   = ev_series(d, idx.outcome, 'str');
EV.out2  = ev_series(d, idx.outcome2, 'str');
EV.crew  = ev_series(d, idx.crew,    'num');
EV.crew2 = ev_series(d, idx.crew2,   'num');
EV.pump  = ev_series(d, idx.pump,    'num');
EV.arc   = ev_series(d, idx.arcflag, 'num');
EV.arc2  = ev_series(d, idx.arcflag2,'num');
EV.aarc  = ev_series(d, idx.agntarc, 'num');
if isempty(EV.tc.t)
    EV.on_t = zeros(0, 1);
else
    tc_inc  = [true; diff(EV.tc.v) > 0];        % target-counter increments = onsets
    EV.on_t = EV.tc.t(tc_inc);
end
if isempty(EV.pump.t) || isempty(EV.crew.t)
    warning('PHY_preprocessing_v3:targets', ...
        'IO_rewardA / CTRL_reward_target_ml not imported - per-target juice will be NaN (re-import with import_flag=true).');
end

% --- Per-variable event series for the cycle loop ---------------------------
% Each variable is extracted from the event table ONCE; a cycle's events are
% then one contiguous index range found by binary search (win_slice), instead
% of ~20 whole-table masks (10^7 events) per cycle. Values keep d.value's
% orientation, so every stim/joy field is identical to d.value(cycIdx & idx.X).
V.task     = var_series(d, idx.task);
V.posX     = var_series(d, d.event == 'STIM_RDP_posX');
V.posY     = var_series(d, d.event == 'STIM_RDP_posY');
V.rdp_dir  = var_series(d, idx.RDP_dir);
V.rdp_coh  = var_series(d, idx.RDP_coh);
V.frame    = var_series(d, idx.frame);
V.trg_on   = var_series(d, idx.trg_on);
V.outcome  = var_series(d, idx.outcome);
V.outcome2 = var_series(d, idx.outcome2);
V.reward   = var_series(d, idx.reward);
V.reward2  = var_series(d, idx.reward2);
V.score    = var_series(d, idx.score);
V.rew_min  = var_series(d, idx.rew_min);
V.rew_max  = var_series(d, idx.rew_max);
V.js_dir   = var_series(d, idx.JS_dir);
V.js_str   = var_series(d, idx.JS_str);
V.js2_dir  = var_series(d, idx.JS2_dir);
V.js2_str  = var_series(d, idx.JS2_str);
V.agnt_dir = var_series(d, idx.AGNT_dir);
V.agnt_str = var_series(d, idx.AGNT_str);
% Joystick samples come in pairs: IO_joystickStrength(2), then its Direction
% 0-0.4 ms later. Windowing both by time can split a pair at a cycle boundary
% (direction one sample short/long -> misaligned; rec077 crashed in
% sort_states). When the pairs are intact, direction is therefore taken at the
% strength samples' indices; else it is windowed by time (with a warning).
js_paired  = is_paired(V.js_dir,  V.js_str);
js2_paired = is_paired(V.js2_dir, V.js2_str);
if ~(js_paired && js2_paired)
    warning('PHY_preprocessing_v3:joystick', ...
        'Joystick direction/strength not logged in pairs - direction windowed by time (may be 1 sample off at cycle edges).');
end

% Trial timestamps
stim.cOn      = double(d.time(idx.cOn));
stim.cEnd     = double(d.time(idx.cEnd));
stim.cpr_cyle = zeros(0, 2);   % CPR cycle table [onset end] (µs), one row per ccnt
ccnt          = 0;

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
    c_on  = stim.cOn(iCyc);
    c_end = stim.cEnd(iCyc);

    stim.task{iCyc} = win_slice(V.task, c_on, c_end);

    if isempty(stim.task{iCyc}) || ~contains(stim.task{iCyc}, 'CPR')
        continue
    end

    ccnt = ccnt + 1;
    stim.cpr_cyle(ccnt, :) = [c_on c_end];

    stim.cpr_solo{ccnt}  = contains(stim.task{iCyc}, 'solo');
    stim.cpr_dyad{ccnt}  = contains(stim.task{iCyc}, 'dyad');
    stim.cpr_catch{ccnt} = contains(stim.task{iCyc}, 'Catch');
    stim.cpr_agent{ccnt} = contains(stim.task{iCyc}, 'computer');   % computer-agent dyad partner

    x_pos = cell2mat(win_slice(V.posX, c_on, c_end));
    y_pos = cell2mat(win_slice(V.posY, c_on, c_end));
    stim.rdp_center_xy{ccnt} = [x_pos y_pos];

    [v, t] = win_slice(V.rdp_dir, c_on, c_end);
    stim.rdp_dir{ccnt}    = cell2mat(v);
    stim.rdp_dir_ts{ccnt} = t;
    [v, t] = win_slice(V.rdp_coh, c_on, c_end);
    stim.rdp_coh{ccnt}    = cell2mat(v);
    stim.rdp_coh_ts{ccnt} = t;
    [~, t] = win_slice(V.frame, c_on, c_end);
    stim.frme_ts{ccnt}    = t;

    [v, tmp_trg_ts]         = win_slice(V.trg_on, c_on, c_end);
    tmp_trg_val             = cell2mat(v);
    stim.feedback_ts{ccnt}  = tmp_trg_ts(tmp_trg_val == 1);
    stim.outcome{ccnt}      = win_slice(V.outcome, c_on, c_end);
    stim.reward_cum{ccnt}   = cellfun(@double, win_slice(V.reward, c_on, c_end));
    stim.score_cum{ccnt}    = cell2mat(win_slice(V.score, c_on, c_end));
    stim.outcome2{ccnt}     = win_slice(V.outcome2, c_on, c_end);   % partner per-target outcome

    % Reward scale in force for this cycle (last value at or before cycle end;
    % falls back to the first logged value). Used to normalise reward to
    % [min max] downstream, since the scale changed between experiments.
    stim.reward_min_ml{ccnt} = last_val_before(V.rew_min, c_end);
    stim.reward_max_ml{ccnt} = last_val_before(V.rew_max, c_end);

    % --- Per-target reward from the cumulative counter (diagnostic only) --
    % The old reconstruction, kept as reward_ind_counter / reward2_ind_counter;
    % reward_ind / reward2_ind are replaced by the time-based extraction below.
    %  * reward_cyc_start = the last INFO_Juice_ml value BEFORE cycle onset (cOn).
    %    Catch-trial rewards and the cycle-completion bonus are dispensed in the
    %    inter-cycle gap, so they sit inside this baseline and are excluded.
    %  * a HIT earns the rise in the counter since the previous hit; a MISS earns
    %    0 and does NOT advance the baseline. Pairs the k-th counter value with
    %    the k-th outcome, which fails from PHY_CPR_20260601 on (see EV block).
    k = bs_lt(V.reward.t, c_on);                    % counter values before cycle onset ...
    if k == 0; base = 0; else; base = double(V.reward.v{k}); end   % ... keep the last
    stim.reward_cyc_start{ccnt} = base;

    rc = stim.reward_cum{ccnt}(:);   oc = stim.outcome{ccnt};   last_hit = base;
    ri = zeros(numel(rc), 1);
    for e = 1:numel(rc)
        if e <= numel(oc) && strcmp(oc{e}, 'hit') && isfinite(rc(e))
            ri(e)    = max(rc(e) - last_hit, 0);   % target juice = rise since last hit
            last_hit = rc(e);
        end                                        % miss / non-hit -> 0, baseline unchanged
    end
    stim.reward_ind{ccnt} = ri;

    % --- Partner (player 2) per-target reward-equivalent (mL), counter -----
    % Same reconstruction from the cumulative partner reward INFO_Juice2_ml
    % (+= CTRL_reward2_target_ml per target) aligned to TRIAL_outcome2; kept as
    % reward2_ind_counter. The computer agent receives no juice, but the reward
    % its performance earned is logged, so this is comparable across partners.
    stim.reward2_cum{ccnt} = cellfun(@double, win_slice(V.reward2, c_on, c_end));
    k = bs_lt(V.reward2.t, c_on);
    if k == 0; base2 = 0; else; base2 = double(V.reward2.v{k}); end
    stim.reward2_cyc_start{ccnt} = base2;

    rc2 = stim.reward2_cum{ccnt}(:);   oc2 = stim.outcome2{ccnt};   last_hit2 = base2;
    ri2 = zeros(numel(rc2), 1);
    for e = 1:numel(rc2)
        if e <= numel(oc2) && strcmp(oc2{e}, 'hit') && isfinite(rc2(e))
            ri2(e)    = max(rc2(e) - last_hit2, 0);
            last_hit2 = rc2(e);
        end
    end
    stim.reward2_ind{ccnt} = ri2;

    % --- Per-target truth (time-based; see the EV block above) ---------------
    % Replaces the index-paired counter reconstruction above, which is kept as
    % *_counter for diagnostics only. All arrays are aligned to feedback_ts.
    stim.reward_ind_counter{ccnt}  = stim.reward_ind{ccnt};
    stim.reward2_ind_counter{ccnt} = stim.reward2_ind{ccnt};
    tg = target_events(EV, stim.feedback_ts{ccnt});
    stim.outcome{ccnt}         = tg.outcome;    % 'hit' | 'miss' | '' per displayed target
    stim.outcome2{ccnt}        = tg.outcome2;
    stim.reward_ind{ccnt}      = tg.juice;      % juice COMMANDED for this target (mL; miss 0)
    stim.reward_npulse{ccnt}   = tg.npulse;     % reward pulses for this target (2 = re-reward)
    stim.reward_ctrl{ccnt}     = tg.ctrl;       % CTRL_reward_target_ml computed at the hit
    stim.reward2_ind{ccnt}     = tg.rew2;       % partner reward-equivalent (mL; miss 0)
    stim.arc_hit{ccnt}         = tg.arc;        % monkey arc covered the target in the window
    stim.arc_ms{ccnt}          = tg.arc_ms;     % ... first time it did (ms from logic onset)
    stim.arc2_hit{ccnt}        = tg.arc2;       % same for the partner
    stim.arc2_ms{ccnt}         = tg.arc2_ms;
    stim.target_logic_ts{ccnt} = tg.t_on;       % logic onset (MWorks us)

    [v, t, i0, i1] = win_slice(V.js_str, c_on, c_end);
    joy.js_monk_tlt{ccnt} = cell2mat(v);
    joy.js_monk_ts{ccnt}  = t;
    if js_paired
        joy.js_monk_dir{ccnt} = cell2mat(V.js_dir.v(i0:i1));   % direction of each strength sample
    else
        joy.js_monk_dir{ccnt} = cell2mat(win_slice(V.js_dir, c_on, c_end));
    end
    % 2nd-player trajectory. For a HUMAN partner this is IO_joystickDirection2 /
    % Strength2. For a COMPUTER-AGENT partner (stim.cpr_agent) the partner is the
    % pre-generated agent (AGNT_direction / AGNT_strength, ~245k events each); map
    % to the js_hum convention as direction = AGNT_direction and tilt =
    % AGNT_strength (0-1, the agent's confidence; = 1 - AGNT_arc_width/180).
    % Resample step-hold onto the 2nd-joystick ~100 Hz grid so downstream
    % sample->ms timing holds. (Requires the AGNT_ vars in the import list; the
    % agent sessions must be re-imported with import_flag=true.)
    [v, ts, i0, i1] = win_slice(V.js2_str, c_on, c_end);
    if stim.cpr_agent{ccnt}
        [ad, at] = win_slice(V.agnt_dir, c_on, c_end);   ad = cell2mat(ad);
        [as, st] = win_slice(V.agnt_str, c_on, c_end);   as = cell2mat(as);
        try
            [atu, ia] = unique(at(:));   [stu, iss] = unique(st(:));
            joy.js_hum_dir{ccnt} = interp1(atu, ad(ia),  ts(:), 'previous', 'extrap');
            joy.js_hum_tlt{ccnt} = interp1(stu, as(iss), ts(:), 'previous', 'extrap');
            joy.js_hum_ts{ccnt}  = ts(:);
        catch
            joy.js_hum_dir{ccnt} = ad(:);   joy.js_hum_tlt{ccnt} = as(:);   joy.js_hum_ts{ccnt} = at(:);
        end
    else
        if js2_paired
            joy.js_hum_dir{ccnt} = cell2mat(V.js2_dir.v(i0:i1));
        else
            joy.js_hum_dir{ccnt} = cell2mat(win_slice(V.js2_dir, c_on, c_end));
        end
        joy.js_hum_tlt{ccnt}  = cell2mat(v);
        joy.js_hum_ts{ccnt}   = ts;
    end
end

% --- Catch trials (INFO_task 'Catch_trial', from rec079) ----------------------
% No random-dot pattern, only the cursor arc(s) centred on fixation: one arc
% (TRIAL_type Catch0_*: STIM_catcharc or STIM_catcharc2) or both (Catch1_*).
% Cursor onset = first STIM_catcharc(2)_onset == 1 within the trial; trials
% aborted before it are skipped (as calc_muae_catch). The 300 ms before cursor
% onset are a baseline window like the 300 ms before a CPR cycle onset; both
% feed the pooled MUAe reference (MUAe block, fn_muae_reference).
EVc   = {ev_series(d, d.event == 'STIM_catcharc_onset',  'num'), ...
         ev_series(d, d.event == 'STIM_catcharc2_onset', 'num')};
EVtyp = ev_series(d, idx.ttype, 'str');
stim.catch_trl    = zeros(0, 2);   % [TRIAL_start TRIAL_end] (µs), one row per catch trial with cursor
stim.catch_cursor = zeros(0, 1);   % cursor onset (µs)
stim.catch_ncur   = zeros(0, 1);   % arcs shown (1 | 2)
stim.catch_type   = cell(0, 1);    % TRIAL_type logged between trial start and cursor onset ('' if none)
for iCat = 1:numel(stim.cEnd)
    if ~any(strcmp(stim.task{iCat}, 'Catch_trial')); continue; end
    trl_on  = stim.cOn(iCat);
    trl_end = stim.cEnd(iCat);
    t_arc   = [first_one(EVc{1}, trl_on, trl_end), first_one(EVc{2}, trl_on, trl_end)];
    if all(isnan(t_arc)); continue; end              % aborted before cursor onset
    t_cur   = min(t_arc);                            % min ignores NaN
    k_typ   = bs_last(EVtyp.t, t_cur);               % last TRIAL_type at/before cursor onset
    stim.catch_trl(end+1, :)    = [trl_on trl_end];
    stim.catch_cursor(end+1, 1) = t_cur;
    stim.catch_ncur(end+1, 1)   = nnz(~isnan(t_arc));
    stim.catch_type{end+1, 1}   = '';
    if k_typ >= 1 && EVtyp.t(k_typ) >= trl_on        % ... logged within this trial
        stim.catch_type{end, 1} = EVtyp.v{k_typ};
    end
end
fprintf('  Behaviour: %d CPR cycles, %d catch trials with cursor (of %d trials).\n', ...
        ccnt, size(stim.catch_trl, 1), numel(stim.cEnd));
clear V


% =========================================================================
%% Import spiking data and run RF characterisation
% =========================================================================

if want_sorted   % skip entirely when 'sorted' is not selected (e.g. MUAe-only)

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
        unit_str = sprintf('unit%d', units(iUnit));
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
    % every unit, including a 300 ms pre-trial baseline period. Cycles =
    % rows of stim.cpr_cyle (CPR trials only, index-aligned with rdp_dir).
    % -----------------------------------------------------------------
    offset = 300e3;   % µs — pre-stimulus baseline window
    for iUnit = 1:numel(units)
        disp(['  CPR — unit ' num2str(iUnit)])
        unit_str = ['unit' num2str(units(iUnit))];
        unitIdx  = spks(:,1) == units(iUnit);

        for iCyc = 1:size(stim.cpr_cyle, 1)
            c_on     = stim.cpr_cyle(iCyc, 1);
            cycleIdx = spks(:,2) >= c_on - offset & spks(:,2) <= stim.cpr_cyle(iCyc, 2);

            % Spike timestamps relative to cycle onset
            brain.CPR.spks.(chan_str).(unit_str){iCyc} = spks(unitIdx & cycleIdx, 2) - c_on;
        end
    end

end % iChan

end % if want_sorted


% =========================================================================
%% Multi-unit activity — MUAe (brain.CPR.muae) + MUAt (spks pseudo-units)
% =========================================================================
% Computed on the CONTINUOUS signal of each channel (broadband from the PL2,
% band-passed in fn_compute_MUA; or datafilt, see mua_source), then sliced into
% cycles — the same "process continuous, then bin" logic as the sorted spikes,
% so MUAe/MUAt share the sorted-spike/behaviour time base. See fn_compute_MUA
% (Stark & Abeles 2007). Drift-robust: no unit isolation required.
% The envelope time base assumes one continuous PL2 recording from t = 0
% (env_t = (0:M-1)/env_fs, synced with the MWorks affine).
% Channels: taken from the local caches (mua_ch###.mat) when any exist, else
% from the PL2's enabled channels — an interrupted extraction therefore needs
% its missing caches recomputed (delete the session's mua_ch*.mat to force).

if want_muat || want_muae

    if ~exist('brain', 'var'); brain = struct(); end

    n_cyc       = size(stim.cpr_cyle, 1);   % CPR cycles (built in the cycle loop)
    n_catch     = size(stim.catch_trl, 1);  % catch trials with cursor (behaviour block)
    mua_offset  = 300e3;   % µs — pre-stimulus baseline window (matches spikes)
    muae_ref_k  = 10;      % pooled MUAe reference: valid baseline windows on each side (fn_muae_reference)
    muae_rf_thr = 5;       % min peak % modulation to attempt an MUAe RF fit
    rf_P        = [];      % RF-mapping presentation table: parsed on the first
                           % channel by fn_RF_responses_muae, reused on the rest

    % --- MUA input source & caching ------------------------------------
    %  'wb'       broadband — extracted on demand from the PL2 (PLX_writeWideband
    %             logic in fn_ensure_wideband) if the *_ch##_wb.mat file is
    %             missing; fn_compute_MUA band-passes it (rebandpass=true).
    %             No datafilt needed for a MUAe-only workflow.
    %  'datafilt' pre-filtered datafilt_ch*.mat (333-5000 Hz); use when comparing
    %             MUAt against the sorter, which ran on datafilt.
    % Caching, cheapest source first, per channel:
    %   mua_ch###.mat  final MUA   -> skip wb read + filtering
    %   *_ch##_wb.mat  broadband   -> skip PL2 extraction
    mua_source    = 'wb';
    use_mua_cache = true;

    switch lower(mua_source)
        case 'wb'
            compute_args = [{'rebandpass', true}, mua_opts];
            pl2_fqn      = [source_dir 'pl2/' fname '.pl2'];
            pl2idx       = [];
            % Prefer channels already available locally (MUA cache or wb), so a
            % cached run needs neither the PL2 nor the Plexon SDK. Only open the
            % PL2 to enumerate channels when nothing local exists.
            dd = [dir([dest_dir 'mua_ch*.mat']); dir([dest_dir '*_ch*_wb.mat'])];
            if ~isempty(dd)
                chan_nums = [];
                for k = 1:numel(dd)
                    t = regexp(dd(k).name, 'ch(\d+)', 'tokens', 'once');
                    if ~isempty(t); chan_nums(end+1) = str2double(t{1}); end %#ok<AGROW>
                end
            elseif isfile(pl2_fqn)
                fn_addpath_plexon();   % adds whichever SDK location exists (no warning)
                pl2idx    = PL2GetFileIndex(pl2_fqn);
                nA        = min(64, numel(pl2idx.AnalogChannels));
                chan_nums = find(cellfun(@(c) c.Enabled, pl2idx.AnalogChannels(1:nA)));
            else
                chan_nums = [];
            end
        case 'datafilt'
            compute_args = mua_opts;
            pl2_fqn = ''; pl2idx = [];
            chan_nums = [];
            dd = dir([dest_dir 'datafilt_ch*.mat']);
            for k = 1:numel(dd)
                t = regexp(dd(k).name, 'ch(\d+)', 'tokens', 'once');
                if ~isempty(t); chan_nums(end+1) = str2double(t{1}); end %#ok<AGROW>
            end
        otherwise
            error('PHY_preprocessing:badSource', 'mua_source must be ''wb'' or ''datafilt''.');
    end
    chan_nums = unique(chan_nums);
    if isempty(chan_nums)
        warning('PHY_preprocessing:noMUAinput', ...
            'signal_source requests MUA but no %s / PL2 input found for %s — skipping MUA.', mua_source, fname);
    end

    % Recompute the cache if any of these change:
    cache_key = struct('source', lower(mua_source), 'do_muat', want_muat, ...
                       'compute_args', {compute_args});

    for ci = 1:numel(chan_nums)
        chan_num     = chan_nums(ci);
        chan_str_mua = sprintf('ch%03d_mua', chan_num);   % 9-char key (fits uid(1:9))
        disp(['  MUA — ch' num2str(chan_num, '%03d')])

        % --- MUAe + MUAt on the full recording (Option A), cached --------
        cache_file = fullfile(dest_dir, sprintf('mua_ch%03d.mat', chan_num));
        m = [];
        if use_mua_cache && isfile(cache_file)
            C = load(cache_file);   % loads m (+ cache_key if present); no warning
            if isfield(C, 'm')
                % Accept a matching cache_key, OR a file that carries none — e.g.
                % MUAe extracted offline by run_muae_extract — provided it holds
                % what this run needs (MUAt crossings, if requested).
                key_ok  = ~isfield(C, 'cache_key') || isequaln(C.cache_key, cache_key);
                muat_ok = ~want_muat || (isfield(C.m, 'params') && ...
                          isfield(C.m.params, 'do_muat') && C.m.params.do_muat);
                if key_ok && muat_ok
                    m = C.m;   % reuse -> skip read + filtering (and PL2)
                end
            end
            clear C
        end
        if isempty(m)
            switch lower(mua_source)
                case 'wb'
                    raw = fn_ensure_wideband(dest_dir, exp_info, chan_num, pl2_fqn, pl2idx, false);   % persist_wb=false (mua cache is the real cache; wb would be 142 GB/session)
                    if isempty(raw); continue; end   % channel not recorded
                    sig = raw.ad_raw_mv;   fs = raw.Fs;   clear raw
                case 'datafilt'
                    S   = load(fullfile(dest_dir, sprintf('datafilt_ch%03d.mat', chan_num)), 'data', 'par');
                    sig = S.data;          fs = S.par.sr;   clear S
            end
            m     = fn_compute_MUA(sig, fs, 'do_muat', want_muat, compute_args{:});
            m.env = single(m.env);   % ample for an envelope; halves cache + memory
            clear sig
            if use_mua_cache
                save(cache_file, 'm', 'cache_key', '-v7.3');
            end
        end

        % --- local data clock -> MWorks µs (same affine as the spikes) ---
        % env sample-time vector is uniform and not stored -> reconstruct it.
        env_t_s    = (0:numel(m.env)-1) / m.env_fs;
        env_t_us   = (env_t_s      * 1e6) * sync_gain + sync_offset;
        muat_ts_us = (m.muat_t_s   * 1e6) * sync_gain + sync_offset;

        % --- MUAt: bin crossings into cycles as pseudo-unit 'muat0' -------
        if want_muat
            for iCyc = 1:n_cyc
                c_on  = stim.cpr_cyle(iCyc, 1);
                c_end = stim.cpr_cyle(iCyc, 2);
                sel   = muat_ts_us >= (c_on - mua_offset) & muat_ts_us <= c_end;
                brain.CPR.spks.(chan_str_mua).muat0{iCyc} = muat_ts_us(sel).' - c_on;
            end
        end

        % --- MUAe: per-cycle / per-catch-trial envelope, baselines, reference
        if want_muae
            brain.CPR.muae.(chan_str_mua).env_fs   = m.env_fs;
            brain.CPR.muae.(chan_str_mua).chan_num = chan_num;
            brain.CPR.muae.(chan_str_mua).env      = cell(1, n_cyc);
            brain.CPR.muae.(chan_str_mua).t_us     = cell(1, n_cyc);
            brain.CPR.muae.(chan_str_mua).baseline = nan(1, n_cyc);

            % env_t_us increases monotonically, so each window is one index
            % range (binary search): baseline [c_on-offset, c_on), segment
            % [c_on-offset, c_end] — the same samples as the logical masks.
            for iCyc = 1:n_cyc
                c_on  = stim.cpr_cyle(iCyc, 1);
                c_end = stim.cpr_cyle(iCyc, 2);

                i0 = bs_lt(env_t_us, c_on - mua_offset) + 1;   % first sample >= c_on - offset
                ib = bs_lt(env_t_us, c_on);                    % last sample  <  c_on
                i1 = bs_last(env_t_us, c_end);                 % last sample  <= c_end

                brain.CPR.muae.(chan_str_mua).baseline(iCyc) = mean(m.env(i0:ib), 'omitnan');
                brain.CPR.muae.(chan_str_mua).env{iCyc}      = single(m.env(i0:i1));
                brain.CPR.muae.(chan_str_mua).t_us{iCyc}     = env_t_us(i0:i1) - c_on;
            end

            % Catch trials: the same windows around cursor onset — baseline
            % [cursor-offset, cursor), segment [cursor-offset, trial end].
            brain.CPR.muae.(chan_str_mua).catch_env      = cell(1, n_catch);
            brain.CPR.muae.(chan_str_mua).catch_t_us     = cell(1, n_catch);
            brain.CPR.muae.(chan_str_mua).catch_baseline = nan(1, n_catch);
            for iCat = 1:n_catch
                c_cur = stim.catch_cursor(iCat);
                c_end = stim.catch_trl(iCat, 2);

                i0 = bs_lt(env_t_us, c_cur - mua_offset) + 1;  % first sample >= cursor - offset
                ib = bs_lt(env_t_us, c_cur);                   % last sample  <  cursor
                i1 = bs_last(env_t_us, c_end);                 % last sample  <= trial end

                brain.CPR.muae.(chan_str_mua).catch_baseline(iCat) = mean(m.env(i0:ib), 'omitnan');
                brain.CPR.muae.(chan_str_mua).catch_env{iCat}      = single(m.env(i0:i1));
                brain.CPR.muae.(chan_str_mua).catch_t_us{iCat}     = env_t_us(i0:i1) - c_cur;
            end

            % Pooled running reference (fn_muae_reference): for every CPR cycle
            % and catch trial, the mean of the muae_ref_k nearest valid baseline
            % windows before and after it (CPR and catch windows pooled,
            % regardless of trial type and text cue), its own window left out.
            % fn_sort_MUAe_by_state normalises by it; the per-trial baselines
            % above stay as raw values.
            muae_ref = fn_muae_reference( ...
                [stim.cpr_cyle(:, 1); stim.catch_cursor], ...
                [brain.CPR.muae.(chan_str_mua).baseline(:); brain.CPR.muae.(chan_str_mua).catch_baseline(:)], ...
                muae_ref_k);
            brain.CPR.muae.(chan_str_mua).baseline_ref       = muae_ref(1:n_cyc).';
            brain.CPR.muae.(chan_str_mua).catch_baseline_ref = muae_ref(n_cyc+1:end).';
            brain.CPR.muae.(chan_str_mua).ref_k              = muae_ref_k;

            % Receptive field from the MUAe envelope (RF-mapping trials).
            [rf_pos, rf_resp, rf_base, rf_tgt, rf_sid, rf_x, rf_y, rf_P] = ...
                fn_RF_responses_muae(d, env_t_us, m.env, [], [], rf_P);
            if ~isempty(rf_pos)
                brain.RF.(chan_str_mua) = fn_characterise_RF(rf_pos, rf_resp, rf_base, ...
                    rf_tgt, rf_sid, rf_x, rf_y, ...
                    struct('response_mode','percent', 'min_response_thr', muae_rf_thr, ...
                           'signal_label','muae'));
                brain.RF.(chan_str_mua).chan_num = chan_num;
            end
        end

        clear m env_t_us muat_ts_us
    end
end


% =========================================================================
%% Format output
% =========================================================================

if ~exist('brain', 'var'); brain = struct(); end
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


% NB: all helpers below are TOP-LEVEL local functions (each after the previous
% function's closing 'end'). A function placed before 'end % RF_characterisation'
% would be NESTED in it and invisible to the main body (crashed every session
% on 2026-09-10).

function S = var_series(d, sel)
% Events of one MWorks variable for the cycle loop: S.t (double µs, sorted) and
% S.v (raw values, same orientation as d.value). MW_readData returns d sorted
% by time, so the order is d's own and a time window is one index range.
S.t = double(d.time(sel));
S.v = d.value(sel);
if ~issorted(S.t)
    [S.t, o] = sort(S.t);   % stable: simultaneous events keep their order
    S.v      = S.v(o);
end
end % var_series


function [v, t, i0, i1] = win_slice(S, a, b)
% Values and times of series S with A <= t <= B — the same elements, order and
% shape (also when empty) as d.value / d.time(d.time >= A & d.time <= B & sel).
% I0:I1 = their indices in S (empty range if none).
i0 = bs_lt(S.t, a) + 1;
i1 = bs_last(S.t, b);
v  = S.v(i0:i1);
t  = S.t(i0:i1);
end % win_slice


function ok = is_paired(A, B)
% True if series A and B are logged in pairs: equal counts and the k-th events
% of both within 2 ms (joystick: one strength/direction pair per ~10 ms sample).
ok = numel(A.t) == numel(B.t) && (isempty(A.t) || max(abs(A.t - B.t)) < 2e3);
end % is_paired


function v = last_val_before(S, t_us)
% Last logged value of series S at/before time T_US; the first logged value if
% none yet; NaN if never logged. For settings that are set once (or rarely) per
% session, e.g. the reward scale CTRL_reward_target_min/max_ml.
v = NaN;
if isempty(S.t); return; end                     % variable never logged -> NaN
k = bs_last(S.t, t_us);                          % last setting at/before t_us
if k == 0; k = 1; end                            % none yet -> first logged setting
try
    v = double(S.v{k});                          % non-numeric value -> stays NaN
catch
end
end % last_val_before


function S = ev_series(d, sel, kind)
% Time-sorted series of one MWorks variable: S.t (us, column) and S.v (double
% column for KIND='num' — non-numeric/empty values -> NaN; cellstr for 'str').
t = double(d.time(sel));   t = t(:);
v = d.value(sel);          v = v(:);
[t, o] = sort(t);          v = v(o);
if strcmp(kind, 'str')
    S.v = repmat({''}, numel(v), 1);
    for k = 1:numel(v)
        if ischar(v{k}) || isstring(v{k}); S.v{k} = char(v{k}); end
    end
else
    S.v = nan(numel(v), 1);
    for k = 1:numel(v)
        x = v{k};
        if (isnumeric(x) || islogical(x)) && ~isempty(x); S.v(k) = double(x(1)); end
    end
end
S.t = t;
end % ev_series


function tg = target_events(EV, fb_ts)
% Per DISPLAYED target (FB_TS = STIM_target_onset times, us) of one cycle: the
% logic onset, both players' outcomes, the juice commanded for a monkey hit,
% the partner's reward-equivalent, and cursor-target overlap from the task's
% arc flags. See the EV block in the behavioural section for the rules.
TGT_US = 50e3;    % target presentation (CTRL_target_duration_ms)
WIN_US = 400e3;   % outcome search window after the onset
LAG_US = 60e3;    % max display-after-logic lag (observed 8-12 ms)
n = numel(fb_ts);
tg.t_on     = nan(n, 1);
tg.outcome  = repmat({''}, n, 1);   tg.outcome2 = repmat({''}, n, 1);
tg.juice    = nan(n, 1);   tg.npulse = zeros(n, 1);   tg.ctrl = nan(n, 1);
tg.rew2     = nan(n, 1);
tg.arc      = nan(n, 1);   tg.arc_ms  = nan(n, 1);
tg.arc2     = nan(n, 1);   tg.arc2_ms = nan(n, 1);
for k = 1:n
    i = bs_last(EV.on_t, fb_ts(k));                  % last logic onset <= display onset
    if i < 1 || fb_ts(k) - EV.on_t(i) > LAG_US; continue; end
    t  = EV.on_t(i);
    nx = t + 30e6;   if i < numel(EV.on_t); nx = EV.on_t(i + 1); end
    w  = min(t + WIN_US, nx);
    tg.t_on(k) = t;
    % monkey: outcome, reward computed at the hit, juice pulses for it
    [o, t_o] = first_str(EV.out, t, w);
    tg.outcome{k} = o;
    t_ref = w;   if ~isnan(t_o); t_ref = t_o; end
    tg.ctrl(k) = last_num(EV.crew, t, t_ref + 5e3);
    if strcmp(o, 'hit') && ~isnan(tg.ctrl(k))
        sel = EV.pump.t > t_o & EV.pump.t < nx & EV.pump.v > 0 & ...
              abs(EV.pump.v - tg.ctrl(k)) < 1e-6;
        tg.juice(k) = sum(EV.pump.v(sel));   tg.npulse(k) = nnz(sel);   % 0 = command was 0
    elseif strcmp(o, 'miss')
        tg.juice(k) = 0;
    end
    % partner: outcome and reward-equivalent
    [o2, t_o2] = first_str(EV.out2, t, w);
    tg.outcome2{k} = o2;
    if strcmp(o2, 'hit')
        tg.rew2(k) = last_num(EV.crew2, t, t_o2 + 5e3);
    elseif strcmp(o2, 'miss')
        tg.rew2(k) = 0;
    end
    % cursor-target overlap (task's own arc flags)
    [tg.arc(k),  tg.arc_ms(k)]  = arc_first(EV.arc,  t, TGT_US);
    [tg.arc2(k), tg.arc2_ms(k)] = arc_first(EV.arc2, t, TGT_US);
    if ~isempty(o2) && ~(tg.arc2(k) == 1) && ~isempty(EV.aarc.t)   % computer agent
        [a, am] = arc_first(EV.aarc, t, TGT_US);
        if a == 1; tg.arc2(k) = a; tg.arc2_ms(k) = am; end
    end
end
end % target_events


function [s, ts] = first_str(S, a, b)
% First 'hit'/'miss' value of S within [a, b] (us); '' / NaN if none.
s = '';   ts = NaN;
i = bs_last(S.t, a - 1) + 1;                         % first event at/after a
while i <= numel(S.t) && S.t(i) <= b
    if any(strcmp(S.v{i}, {'hit', 'miss'})); s = S.v{i}; ts = S.t(i); return; end
    i = i + 1;
end
end % first_str


function v = last_num(S, a, b)
% Last value of S within [a, b] (us); NaN if none.
v = NaN;
i = bs_last(S.t, b);
if i >= 1 && S.t(i) >= a; v = S.v(i); end
end % last_num


function [hit, ms] = arc_first(S, t, dur_us)
% Did the arc flag S cover the target in [t, t + dur_us]? HIT = 1/0 (NaN if the
% flag was never logged), MS = first time it did (ms from t; 0 = at onset).
hit = NaN;   ms = NaN;
if isempty(S.t); return; end
i = bs_last(S.t, t);
if i >= 1 && S.v(i) == 1; hit = 1; ms = 0; return; end
hit = 0;
j = i + 1;
while j <= numel(S.t) && S.t(j) <= t + dur_us
    if S.v(j) == 1; hit = 1; ms = (S.t(j) - t) / 1e3; return; end
    j = j + 1;
end
end % arc_first


function i = bs_last(t, x)
% Index of the last element of sorted T with T <= X (0 if none), i.e. the
% number of elements <= X. Binary search.
lo = 1;   hi = numel(t);   i = 0;
while lo <= hi
    mid = floor((lo + hi) / 2);
    if t(mid) <= x; i = mid; lo = mid + 1; else; hi = mid - 1; end
end
end % bs_last


function i = bs_lt(t, x)
% Number of elements of sorted T with T < X (the last one's index; 0 if none).
lo = 1;   hi = numel(t);   i = 0;
while lo <= hi
    mid = floor((lo + hi) / 2);
    if t(mid) < x; i = mid; lo = mid + 1; else; hi = mid - 1; end
end
end % bs_lt


function t1 = first_one(S, a, b)
% Time of the first value == 1 of series S (ev_series 'num') within [a, b] (us);
% NaN if none. Catch trials: cursor onset from STIM_catcharc(2)_onset.
t1 = NaN;
i  = bs_lt(S.t, a) + 1;                              % first event at/after a
while i <= numel(S.t) && S.t(i) <= b
    if S.v(i) == 1; t1 = S.t(i); return; end
    i = i + 1;
end
end % first_one
