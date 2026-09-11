%% CPR main analysis — MUAe (v4)
%
% Slim, MUAe-only driver. No spike-sorting / MUAt code: it extracts the
% multi-unit envelope, gates channels, aligns behaviour to each state, and
% stores the per-state envelope. For the sorted-spike or MUAt path use
% PHY_main_analysis_v3 (which carries the signal_source switch).
%
% Required per recording:
%   .MWK2 / .H5           behaviour (MWorks)
%   PL2 or *_wb.mat or mua_ch*.mat   neural input (auto: cache -> wb -> PL2)
%
% Pipeline:
%   1. PHY_preprocessing_v3(..., {'muae'})  MUAe envelope + RF + per-cycle store
%   2. fn_muae_quality                      RF-fit / onset-response inclusion gate
%   3. sort_states (local)                  behaviour + state boundaries (no spikes)
%   4. fn_sort_MUAe_by_state                per-state envelope (baseline-normalised)
%
% Outputs (per recording, under dest_dir):
%   summary_<rec>_muae.mat          phy  (brain.CPR.muae, brain.RF, muae.include)
%   state_responses_<rec>_muae.mat  state (behaviour + state.muae / muae_mean)
%
% Felix Schneider, CNL
%
% Version history
%   4.0  (fxs 2026-07-03)  MUAe-only rewrite — all spike code removed.
%   4.1  (fxs 2026-08-11)  CNL04R deploy: /Users/cnl paths; rec_lst = 9 clean
%                          sessions (060,062,076-079,090-092). persist_wb=false
%                          in fn_ensure_wideband (no 142 GB/session wb writes).
%                          Skipped (naming/data issues): 063 no h5; 067 pl2
%                          block1/rec_067 vs h5 block3/rec067; 069 pl2 fxs vs h5 ann.
%   4.2  (fxs 2026-09-10)  29-session rec_lst; rec_behav_only (PL2 copies
%                          truncated) processed behaviour-only; any failed
%                          session retried behaviour-only; progress/error logs;
%                          sort_states adds partner outcome/reward, reward
%                          scale and the time-based per-target fields
%                          (reward_npulse, reward_ctrl, arc_*, target_logic_ts,
%                          state.reward_src = 'IO_rewardA').
%   4.3  (fxs 2026-09-11)  sort_states pads joystick arrays (legacy summaries
%                          with a split direction/strength pair crashed:
%                          rec077); minutes per recording in the progress log.
%                          Uses PHY_preprocessing_v3 2.3 (~25x faster MUAe).
%
% Run notes
%   import_flag = true re-reads the .mwk2 (~15-20 min/session) and rewrites the
%   .h5. Needed once for h5 files written before 2026-09-11 (they lack the
%   CTRL_ variables: reward scale, per-target reward, arc flags). The extended
%   felix_nhp_solo.cfg keeps them, so later runs can use import_flag = false.
%   Progress: <out_root>/preprocessing_progress.log; errors (full report):
%   <out_root>/preprocessing_errors.log.

clear; close all;

addpath('/Users/cnl/Desktop/CPR/code');
addpath(genpath('/Users/cnl/Documents/GitLab/matlab4mworks'));

% =========================================================================
%% Configuration
% =========================================================================

cfg_pth      = '/Users/cnl/Desktop/CPR/code/felix_nhp_solo.cfg';
source_dir   = '/Users/cnl/Documents/DATA/Nilan/';
out_root     = '/Users/cnl/Documents/DATA/Nilan/muae/';
preproc_flag = true;     % true = run preprocessing+quality; false = load saved summary
reprocess    = true;     % re-process even if outputs exist (bypass resumable skip); set false for normal runs
import_flag  = true;     % true = import .mwk2 -> .h5; false = load existing .h5 (see Run notes)

% Per-cycle baseline normalisation applied to state.muae:
%   'divide'   env/baseline   (>= 0, gain-like — best for DSI/CV tuning metrics)
%   'percent'  100*(env-bl)/bl (signed modulation)
%   'subtract' env - bl        'none' raw envelope
muae_norm    = 'percent';


% =========================================================================
%% Recording list
% =========================================================================

rec_lst = {
    '20250923_nil_CPR_block1_phy4_rec058_fxs', ...
    '20250924_nil_CPR_block1_phy4_rec059_fxs', ...
    '20250925_nil_CPR_block1_phy4_rec060_fxs', ...
    '20250926_nil_CPR_block1_phy4_rec061_fxs', ...
    '20250930_nil_CPR_block1_phy4_rec062_ann', ...
    '20251204_nil_CPR_block1_phy4_rec068_fxs', ...
    '20251205_nil_CPR_block1_phy4_rec069_fxs', ...   % behav-only
    '20251211_nil_CPR_block1_phy4_rec070_ann', ...   % behav-only
    '20260114_nil_CPR_block1_phy4_rec071_ann', ...
    '20260128_nil_CPR_block1_phy4_rec074_ann', ...
    '20260205_nil_CPR_block1_phy4_rec075_ann', ...
    '20260313_nil_CPR_block1_phy4_rec076_fxs', ...   % behav-only
    '20260320_nil_CPR_block1_phy4_rec077_fxs', ...   % behav-only
    '20260409_nil_CPR_block1_phy4_rec078_fxs', ...   % behav-only
    '20260429_nil_CPR_block1_phy4_rec079_fxs', ...   % behav-only
    '20260430_nil_CPR_block1_phy4_rec080_ann', ...
    '20260515_nil_CPR_block1_phy4_rec081_ann', ...   % behav-only
    '20260520_nil_CPR_block1_phy4_rec082_ann', ...
    '20260521_nil_CPR_block1_phy4_rec083_ann', ...
    '20260522_nil_CPR_block1_phy4_rec084_fxs', ...   % behav-only
    '20260528_nil_CPR_block1_phy4_rec086_fxs', ...   % behav-only
    '20260529_nil_CPR_block1_phy4_rec087_fxs', ...   % behav-only
    '20260603_nil_CPR_block1_phy4_rec088_fxs', ...
    '20260610_nil_CPR_block1_phy4_rec089_ann', ...
    '20260612_nil_CPR_block1_phy4_rec090_ann', ...
    '20260617_nil_CPR_block1_phy4_rec091_ann', ...
    '20260618_nil_CPR_block1_phy4_rec092_fxs', ...
    '20260625_nil_CPR_block1_phy4_rec093_fxs', ...
    '20260626_nil_CPR_block1_phy4_rec094_fxs'};

% Sessions whose PL2 cannot be synced (checked 2026-09-10) are processed
% BEHAVIOUR-ONLY (no sync, no PL2 read, no MUAe): the PL2 copies on cnl04r AND
% on the lab server (DATA/fxs/CPR_electrophysiology/Nilan/pl2) are full-size but
% zero-filled after the first 0.3-14 GB (the surviving part ends before the
% first CPR trial), so the PL2 footer/index is missing and PL2GetFileIndex
% reports 0 events on 'Strobed'/EVT01-16 -> MW_getSyncParam fails. MWorks DID
% send the strobes and the Omniplex acknowledged them during recording
% ('Strobe_ts ... return OPX_ts' in the .mwk2), and each copy's header still
% points to a footer, i.e. the ORIGINAL recordings were closed properly — the
% originals on the recording PC may be intact. Move a session out of this list
% once its PL2 is replaced to get its MUAe. (rec081 is NOT here: its PL2 is
% intact; it only failed because its .mwk2 name has a stray space, which
% PHY_preprocessing_v3 now resolves.)
rec_behav_only = {
    '20251205_nil_CPR_block1_phy4_rec069_fxs', ...
    '20251211_nil_CPR_block1_phy4_rec070_ann', ...
    '20260313_nil_CPR_block1_phy4_rec076_fxs', ...
    '20260409_nil_CPR_block1_phy4_rec078_fxs', ...
    '20260429_nil_CPR_block1_phy4_rec079_fxs', ...
    '20260522_nil_CPR_block1_phy4_rec084_fxs', ...
    '20260528_nil_CPR_block1_phy4_rec086_fxs', ...
    '20260529_nil_CPR_block1_phy4_rec087_fxs'};

% Plain-text logs next to the outputs (readable remotely): one line per
% recording start/finish, and the full MATLAB error report for failures.
log_prog = fullfile(out_root, 'preprocessing_progress.log');
log_err  = fullfile(out_root, 'preprocessing_errors.log');
write_log(log_prog, sprintf('=== run started (%d recordings, %d behaviour-only) ===', numel(rec_lst), numel(rec_behav_only)));

% =========================================================================
%% Main processing loop
% =========================================================================

failed = {};   % recordings that errored -- skipped, not fatal

for iRec = 1:numel(rec_lst)

    rec_name = rec_lst{iRec};
    rec_info = split(rec_name, '_');
    behav_only = ismember(rec_name, rec_behav_only);
    write_log(log_prog, sprintf('%s  started%s', rec_name, repmat(' (behaviour-only)', 1, double(behav_only))));
    t_rec = tic;

    dest_dir     = fullfile(out_root, sprintf('%s_%s_%s', rec_info{1}, rec_info{6}, rec_info{4}), '/');
    summary_file = fullfile(dest_dir, ['summary_'         rec_name '_muae.mat']);
    state_file   = fullfile(dest_dir, ['state_responses_' rec_name '_muae.mat']);

    fprintf('\n%s\n  Recording %d / %d:  %s\n%s\n', ...
        repmat('=', 1, 70), iRec, numel(rec_lst), rec_name, repmat('=', 1, 70));

    % Resumable: skip recordings whose outputs already exist.
    if preproc_flag && ~reprocess && isfile(summary_file) && isfile(state_file)
        fprintf('  Already processed -- skipping.\n');
        continue
    end
    try

    if ~behav_only
        % --- Stage 1: MUAe preprocessing + quality gate ------------------
        if preproc_flag
            fprintf('  [1/2]  MUAe preprocessing + quality...\n');
            phy = PHY_preprocessing_v3(rec_name, source_dir, dest_dir, cfg_pth, import_flag, {'muae'});
            phy = fn_muae_quality(phy);
            save(summary_file, 'phy', '-v7.3');
            fprintf('         Saved: %s\n', summary_file);
        else
            fprintf('  [1/2]  Loading saved summary...\n');
            load(summary_file, 'phy');
        end

        % --- Stage 2: state alignment + per-state envelope ---------------
        fprintf('  [2/2]  State alignment + MUAe...\n');
        state = sort_states(phy);
        state = fn_sort_MUAe_by_state(state, phy.brain, phy.stim, 'norm', muae_norm);
    else
        [phy, state] = behaviour_only(rec_name, source_dir, dest_dir, cfg_pth, import_flag);
        save(summary_file, 'phy', '-v7.3');
    end

    save(state_file, 'state', '-v7.3');
    fprintf('         Saved: %s\n', state_file);
    write_log(log_prog, sprintf('%s  DONE%s (%.0f min)', rec_name, ...
        repmat(' (behaviour-only)', 1, double(behav_only)), toc(t_rec) / 60));


    catch ME
        warning('Recording %s FAILED (%s).', rec_name, ME.message);
        write_log(log_prog, sprintf('%s  FAILED after %.0f min: %s', rec_name, toc(t_rec) / 60, ME.message));
        write_log(log_err,  sprintf('%s\n%s\n', rec_name, getReport(ME, 'extended', 'hyperlinks', 'off')));
        % Behaviour must not be lost with the neural part: retry behaviour-only.
        ok_b = false;
        if ~behav_only
            try
                [phy, state] = behaviour_only(rec_name, source_dir, dest_dir, cfg_pth, import_flag);
                save(summary_file, 'phy', '-v7.3');
                save(state_file,   'state', '-v7.3');
                ok_b = true;
                write_log(log_prog, sprintf('%s  DONE (behaviour-only fallback, %.0f min)', rec_name, toc(t_rec) / 60));
            catch ME2
                write_log(log_prog, sprintf('%s  behaviour-only fallback FAILED: %s', rec_name, ME2.message));
                write_log(log_err,  sprintf('%s (behaviour-only fallback)\n%s\n', rec_name, ...
                                            getReport(ME2, 'extended', 'hyperlinks', 'off')));
            end
        end
        if ~ok_b; failed{end+1} = rec_name; end   %#ok<SAGROW>
    end

end

if isempty(failed)
    fprintf('\nAll recordings processed.\n');
else
    fprintf('\nDone. %d recording(s) FAILED and were skipped:\n', numel(failed));
    fprintf('   %s\n', failed{:});
end
write_log(log_prog, sprintf('=== run finished: %d failed ===', numel(failed)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function state = sort_states(in)
% SORT_STATES  Align behaviour to each RDP direction-change (state) within a
% CPR cycle: stimulus parameters, boundaries, feedback/outcome, and joystick
% trajectories. Behaviour only — the envelope is attached afterwards by
% fn_sort_MUAe_by_state. (Slimmed from PHY_sort_spikes_by_state; no spikes.)

clear state
state_cnt = 0;
stim      = in.stim;
joy       = in.joy;

for iCyc = 1:numel(stim.rdp_dir)
    for iState = 2:numel(stim.rdp_dir{iCyc})   % discard inherited direction from previous cycle

        state_cnt = state_cnt + 1;

        % --- Task flag ---------------------------------------------------
        state.cpr_solo(state_cnt) = stim.cpr_solo{iCyc};

        % --- Stimulus parameters -----------------------------------------
        state.rdp_dir(state_cnt)    = stim.rdp_dir{iCyc}(iState);
        state.rdp_dir_ts(state_cnt) = stim.rdp_dir_ts{iCyc}(iState);

        % Coherence is inherited from the most recent coherence event
        state.rdp_coh(state_cnt) = stim.rdp_coh{iCyc}( ...
            find(state.rdp_dir_ts(state_cnt) >= stim.rdp_coh_ts{iCyc}, 1, 'last'));

        % State boundaries: [onset timestamp, offset timestamp]
        is_last_state = iState == numel(stim.rdp_dir{iCyc});
        if is_last_state
            state_end = double(stim.cpr_cyle(iCyc, 2));
        else
            state_end = stim.rdp_dir_ts{iCyc}(iState + 1);
        end
        state.boundaries(state_cnt, :) = [stim.rdp_dir_ts{iCyc}(iState), state_end];
        state.dur_s(state_cnt)         = double(state_end - stim.rdp_dir_ts{iCyc}(iState)) / 1e6;
        state.cIdx(state_cnt)          = iCyc;

        % --- Feedback and outcome ----------------------------------------
        trg_idx = stim.feedback_ts{iCyc} > stim.rdp_dir_ts{iCyc}(iState) & ...
                  stim.feedback_ts{iCyc} < state_end;

        state.feedback_state_ts_raw{state_cnt} = stim.feedback_ts{iCyc}(trg_idx);

        dt_us = stim.feedback_ts{iCyc}(trg_idx) - stim.rdp_dir_ts{iCyc}(iState);
        state.feedback_state_ms{state_cnt}    = dt_us ./ 1e3;
        state.feedback_state_smple{state_cnt} = ceil((dt_us ./ 1e3) ./ (1000/120));

        % outcome / reward_ind (and the partner's) are built per displayed target
        % by the time-based extraction, so they are parallel to feedback_ts by
        % construction. The raw logs reward_cum/score_cum can be SHORTER for a
        % cycle (a target whose counter entry was not logged); pad to the feedback
        % length before indexing so it becomes NaN instead of an out-of-bounds crash.
        nfb = numel(stim.feedback_ts{iCyc});
        state.outcome{state_cnt}    = pad_index(stim.outcome{iCyc},    trg_idx, nfb);
        state.reward_cum{state_cnt} = pad_index(stim.reward_cum{iCyc}, trg_idx, nfb);
        state.reward_ind{state_cnt} = pad_index(stim.reward_ind{iCyc}, trg_idx, nfb);   % per-target juice (mL)
        state.score_cum{state_cnt}  = pad_index(stim.score_cum{iCyc},  trg_idx, nfb);
        state.reward_cyc_start(state_cnt) = stim.reward_cyc_start{iCyc};                % cumulative at cycle onset

        % --- Partner (player 2) outcome + reward-equivalent, and reward scale --
        % Partner fields mirror the monkey's, so per-player hit rate / reward can
        % be computed identically. reward_min/max_ml give the reward scale in
        % force, for normalising reward to [0 1] across sessions.
        if isfield(stim, 'outcome2')
            state.outcome2{state_cnt}    = pad_index(stim.outcome2{iCyc},    trg_idx, nfb);
            state.reward2_ind{state_cnt} = pad_index(stim.reward2_ind{iCyc}, trg_idx, nfb);
            state.reward2_cyc_start(state_cnt) = stim.reward2_cyc_start{iCyc};
        end
        if isfield(stim, 'reward_min_ml')
            state.reward_min_ml(state_cnt) = stim.reward_min_ml{iCyc};
            state.reward_max_ml(state_cnt) = stim.reward_max_ml{iCyc};
        end

        % --- Per-target truth (time-based, PHY_preprocessing_v3 2026-09-10) ----
        % outcome / outcome2 / reward_ind / reward2_ind above already come from
        % the time-based extraction; these add the audit fields (all aligned to
        % the state's targets): juice pulses per target, the reward the task
        % computed at the hit, cursor-target overlap for both players (the task's
        % arc flag within the 50 ms target window; in dyads a 'miss' can hide a
        % later overlap, because whoever covers the target first ends it), and
        % the logic onset. Old stim structs lack them -> fields absent.
        if isfield(stim, 'reward_npulse')
            state.reward_npulse{state_cnt}   = pad_index(stim.reward_npulse{iCyc},   trg_idx, nfb);
            state.reward_ctrl{state_cnt}     = pad_index(stim.reward_ctrl{iCyc},     trg_idx, nfb);
            state.arc_hit{state_cnt}         = pad_index(stim.arc_hit{iCyc},         trg_idx, nfb);
            state.arc_ms{state_cnt}          = pad_index(stim.arc_ms{iCyc},          trg_idx, nfb);
            state.arc2_hit{state_cnt}        = pad_index(stim.arc2_hit{iCyc},        trg_idx, nfb);
            state.arc2_ms{state_cnt}         = pad_index(stim.arc2_ms{iCyc},         trg_idx, nfb);
            state.target_logic_ts{state_cnt} = pad_index(stim.target_logic_ts{iCyc}, trg_idx, nfb);
            state.reward_ind_counter{state_cnt} = pad_index(stim.reward_ind_counter{iCyc}, trg_idx, nfb);
        end

        % --- Joystick trajectories ---------------------------------------
        % tilt/dir are parallel to *_ts (PHY_preprocessing_v3 >= 2.3). Summaries
        % from older versions can have direction one sample short in a cycle
        % (pair split at the cycle edge) -> pad_index gives NaN, no crash.
        js_monk_idx = joy.js_monk_ts{iCyc} > stim.rdp_dir_ts{iCyc}(iState) & ...
                      joy.js_monk_ts{iCyc} < state_end;
        js_hum_idx  = joy.js_hum_ts{iCyc}  > stim.rdp_dir_ts{iCyc}(iState) & ...
                      joy.js_hum_ts{iCyc}  < state_end;
        n_mj = numel(joy.js_monk_ts{iCyc});
        n_hj = numel(joy.js_hum_ts{iCyc});

        state.js_monk_tlt{state_cnt} = pad_index(joy.js_monk_tlt{iCyc}, js_monk_idx, n_mj);
        state.js_monk_dir{state_cnt} = pad_index(joy.js_monk_dir{iCyc}, js_monk_idx, n_mj);
        state.js_hum_tlt{state_cnt}  = pad_index(joy.js_hum_tlt{iCyc},  js_hum_idx,  n_hj);
        state.js_hum_dir{state_cnt}  = pad_index(joy.js_hum_dir{iCyc},  js_hum_idx,  n_hj);

    end % iState
end % iCyc

% Provenance of the per-target reward: 'IO_rewardA' = juice commanded per target,
% extracted by time from the reward-board pulses (PHY_preprocessing_v3,
% 2026-09-10). Absent -> older index-paired counter reconstruction.
if isfield(stim, 'reward_npulse'); state.reward_src = 'IO_rewardA'; end

end % sort_states


function y = pad_index(x, idx, n)
% Apply logical index IDX (length N) to X, first padding X to length N if it is
% shorter (numeric -> NaN, cell -> ''). Guards sort_states against per-cycle
% arrays that are shorter than the array IDX was built from (reward/score vs
% feedback, joystick direction vs strength samples). Equal lengths: X(IDX).
if numel(x) < n
    if iscell(x); x(end+1:n) = {''}; else; x(end+1:n) = NaN; end
end
y = x(idx);
end % pad_index


function write_log(fn, msg)
% Append a time-stamped line (or block) to a plain-text log; never fatal.
try
    fid = fopen(fn, 'a');
    if fid < 0; return; end
    fprintf(fid, '[%s] %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), msg);
    fclose(fid);
catch
end
end % write_log


function [phy, state] = behaviour_only(rec_name, source_dir, dest_dir, cfg_pth, import_flag)
% Behaviour-only processing: MWorks import + behavioural block (no sync, no PL2,
% no MUAe) and the state alignment. The state carries every behavioural field
% (outcomes, per-target juice, joystick, ...) but no muae / muae_mean, which is
% how the analysis pipeline recognises a behaviour-only session.
fprintf('  [1/2]  Behaviour-only preprocessing (no PL2)...\n');
phy   = PHY_preprocessing_v3(rec_name, source_dir, dest_dir, cfg_pth, import_flag, {'behaviour'});
fprintf('  [2/2]  State alignment (behaviour only)...\n');
state = sort_states(phy);
end % behaviour_only
