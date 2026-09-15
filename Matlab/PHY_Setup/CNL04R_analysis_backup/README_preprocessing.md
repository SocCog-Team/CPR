# CPR MUAe preprocessing on cnl04r

*Updated 2026-09-15. Code is in `/Users/cnl/Desktop/CPR/code`. The versions before the pooled MUAe baseline are in `_backup_2026-09-15_pre_pooled_baseline/`; those before the 2026-09-10 speed-up are in `_backup_2026-09-10_pre_speedup/`.*

## How to run

1. **Pre-flight check.** Run `preflight_pooled_baseline` in MATLAB on cnl04r once before a full run; it takes about 5 min. It processes rec094 in memory with `import_flag = false` and writes only `preflight_pooled_baseline_rec094.txt` and `.mat` to `/Users/cnl/Documents/DATA/Nilan/muae/`. The report ends with one verdict per check, and all four should read `PASS`:
   1. Loading the `.h5` gives the same summary as the 2026-09-11 run from the `.mwk2`, apart from the new fields.
   2. The old per-cycle normalisation is reproduced exactly.
   3. The pooled normalisation is the exact re-expression of the per-cycle one.
   4. The stored references match a brute-force computation, and each catch trial's arc count matches its `TRIAL_type`.
2. **Full run.** Run `PHY_main_analysis_v4` with these settings:
   - `import_flag = false` if check 1 passed, otherwise `true`. All 29 `.h5` files were rewritten on 2026-09-11 and contain every variable the preprocessing reads (checked 2026-09-15).
   - `reprocess = true`
   - `preproc_flag = true`
3. **Progress.** Follow `/Users/cnl/Documents/DATA/Nilan/muae/preprocessing_progress.log`, which has one line per recording start and end, with minutes. `preprocessing_errors.log` in the same folder holds the full MATLAB error report of any failure.
4. **Expected time.** With `import_flag = true`, reading the `.mwk2` dominates at about 15–20 min per session. With `import_flag = false`, a session loads from the `.h5` in seconds. Slicing the cached envelopes takes about 15 s per 64-channel session; loading the channel caches and saving the summary and state files add a few minutes.
5. **Afterwards.** The laptop watcher copies the state and summary files to `~/Desktop/muae/` and validates the per-target fields against `~/Desktop/muae/target_truth`. Then rebuild the v11 shards (`run_build = true`). The new state files carry `state.muae_baseline = 'pooled'`. `helper_muae_renormalise` leaves such states unchanged, but the v11 driver warns until `cfg.muae_baseline` is switched to the preprocessed normalisation.

## What runs

`PHY_main_analysis_v4` processes each recording in `rec_lst` in one of two modes.

**Full mode (MUAe)**
1. `PHY_preprocessing_v3(rec, …, {'muae'})`
2. `fn_muae_quality`, then save `summary_<rec>_muae.mat`
3. `sort_states`
4. `fn_sort_MUAe_by_state` (normalised to the pooled reference), then save `state_responses_<rec>_muae.mat`

**Behaviour-only mode** is used for the sessions in `rec_behav_only`. It is also the automatic fallback when any session fails.
1. `PHY_preprocessing_v3(rec, …, {'behaviour'})`. This does no MWorks–Plexon sync and never reads the PL2.
2. `sort_states`. The resulting state has no `muae*` fields, which is how the v11 analysis recognises a behaviour-only session.

`PHY_preprocessing_v3` works in stages:

| Stage | What it does |
|---|---|
| Import | `.mwk2` → event struct `d` (variables in `var_import`), and the `.h5` cache is rewritten. With `import_flag = false`, it reads the `.h5` instead. A failed `.h5` write only warns. |
| Sync | `syncParam_<rec>.mat`, or `MW_getSyncParam` on the PL2. Skipped for behaviour-only. |
| Behaviour | Cycle loop over CPR trials. Produces `stim.*` and `joy.*` per cycle and `stim.cpr_cyle` (cycle onset and end). Per-target outcome, juice and overlap are extracted by time (`target_events`). Catch trials with a cursor go into `stim.catch_*`. |
| MUAe | Per channel, the cached envelope (`mua_ch###.mat`, otherwise computed from the PL2 wideband) is sliced per CPR cycle and per catch trial, each with the 300 ms baseline window before onset. The pooled reference is computed from these windows (`fn_muae_reference`). The RF is estimated from the RF-mapping trials of the envelope (`fn_RF_responses_muae` → `fn_characterise_RF`). |

## MUAe baseline normalisation (2026-09-15)

**Why.** Normalising each cycle to its own pre-cycle baseline carried trial-type differences of that baseline into every state response. The baseline differs between solo, dyad and catch trials, because fixation timing, joystick hold and the text cue differ. In MT this inflated dyad − solo from +0.12 to +0.56 percentage points (all directions).

**Definition.** Per channel, every baseline window enters one time-ordered pool:
- CPR cycles: mean envelope over [`TRIAL_start` − 300 ms, `TRIAL_start`);
- catch trials: mean envelope over [cursor onset − 300 ms, cursor onset).

A window's reference is the mean of the 10 nearest valid windows before it and the 10 nearest valid windows after it, regardless of trial type and text cue. The window itself is left out, and at the session edges fewer windows are used. A window is valid if its mean is finite and > 0. `fn_sort_MUAe_by_state` then computes `100 · (envelope − reference) / reference`.

**Choice of K.** The comparison used 19 sessions with 450 MT and 153 IPS channels (`reference_window_choice.py`, in the SfN handover folder). The criterion was the SD of normalised cycle responses within solo and within dyad cycles: lower means the reference removes more drift without adding noise. Values are relative to K = 10:

| Reference | MT | IPS |
|---|---|---|
| Own per-cycle baseline | +44 % | +29 % |
| Whole-session mean | +54 % | +57 % |
| K = 3 | +4.7 % | −1.0 % |
| K = 5 | +2.8 % | −0.7 % |
| **K = 10** | **0** | **0** |
| K = 15 | +3.5 % | +3.2 % |
| K = 20 | +6.1 % | +6.3 % |
| K = 30 | +17 % | +15 % |

- Using fewer windows at the session edges was slightly better than reaching further into the session (MT −0.3 %, p = 0.03).
- The mean was better than the median (MT +1.5 % for the median).
- Solo, dyad and catch trials are randomly interleaved (median run length 1). A reference therefore holds each trial type at its overall share, for every K.
- With every pooled reference, MT dyad − solo was +0.02 to +0.15 pp and not significant, against +0.57 pp (17 of 19 sessions) with the own baseline.

**Fields.**

| Where | Field | Meaning |
|---|---|---|
| `stim` | `catch_trl` | [start end] of each catch trial with a cursor (µs) |
| | `catch_cursor` | Cursor onset: the first `STIM_catcharc_onset` or `STIM_catcharc2_onset` = 1 in the trial (µs) |
| | `catch_ncur`, `catch_type` | Arcs shown (1 or 2); `TRIAL_type` |
| `brain.CPR.muae.(ch)` | `baseline`, `catch_baseline` | Raw window means (mV) |
| | `catch_env`, `catch_t_us` | Catch-trial envelope from 300 ms before cursor onset to trial end; times relative to cursor onset |
| | `baseline_ref`, `catch_baseline_ref` | Pooled reference per cycle and per catch trial (mV) |
| | `ref_k` | Windows on each side (10) |
| `state` | `muae_baseline`, `muae_ref_k` | `'pooled'` and 10 |

`fn_sort_MUAe_by_state(…, 'baseline', 'cycle')` reproduces the old normalisation, and `'ref_k'` recomputes the reference from the stored windows. Summaries written before `PHY_preprocessing_v3` 2.4 have no catch windows, so their pooled reference uses CPR cycles only, with a warning. The inclusion test in `fn_muae_quality` still compares each cycle's onset response with that cycle's own baseline, so channel inclusion does not change.

## Per-target fields in `state`

All fields are aligned to the state's displayed targets (`feedback_state_ts_raw`).

| Field | Meaning |
|---|---|
| `outcome`, `outcome2` | Monkey and partner `hit` / `miss`: the first `TRIAL_outcome(2)` after the target onset |
| `reward_ind` | Juice commanded for the target, in mL. These are the `IO_rewardA` pulses after the hit that match `CTRL_reward_target_ml`. A miss is 0. |
| `reward_npulse` | Pulses per hit. 2 means a re-reward after a fixation break (a task bug). |
| `reward_ctrl` | `CTRL_reward_target_ml` at the hit |
| `reward2_ind` | Partner reward-equivalent (`CTRL_reward2_target_ml`). The partner receives no juice. |
| `arc_hit`, `arc_ms`, `arc2_hit`, `arc2_ms` | Whether the cursor arc covered the target within its 50 ms, and when (task arc flags) |
| `target_logic_ts` | Logic onset of the target (`INFO_TargetCounter` increment). The display follows 8–12 ms later. |
| `reward_ind_counter` | Old reconstruction from the `INFO_Juice_ml` counter. Kept for diagnostics only; it is wrong for 13–30 % of targets in rec088–094. |
| `reward_min_ml`, `reward_max_ml` | Reward scale per state |
| `state.reward_src = 'IO_rewardA'` | Marks the time-based extraction |

The validation against the raw `.mwk2` is described in `~/Desktop/muae/target_truth/README.md` on the laptop.

## Data caveats

- **`.h5` cache.** `MW_writeH5` keeps only variables defined in `MW_readFile.cfg` or `felix_nhp_solo.cfg`. It stores floats as single precision and renames `IO_rewardA` to `IO_rewardA_ml`; the preprocessing accepts both names. Any variable newly added to `var_import` must also be added to `felix_nhp_solo.cfg`. Lines in that file must be of the form `"NAME"  0  type  ""`; comment lines are not allowed.
- **Truncated PL2 copies.** The copies of rec069, 070, 076, 078, 079, 084, 086 and 087 are truncated, and so is rec085's (not in `rec_lst`); see the IT report. The originals were deleted from the recording PC for space reasons and cannot be recovered (confirmed 2026-09-11). These sessions therefore stay behaviour-only for good.
- **Last files from the recording PC (checked 2026-09-11).** rec091, rec093, rec094 and rec095 were copied from the recording PC to an external drive. The copies on cnl04r and on the lab server are intact and byte-identical to them at every GiB, so nothing was replaced. The server holds rec095 twice, as `…_CPR_block1_…_rec095_ann.pl2` and `…_CPRagent_block1_…_rec095_ann.pl2`; the two files are identical.
- **rec095 has no behaviour file (checked 2026-09-15).** Its PL2 is in `DATA/Nilan/pl2/`, and its `muae/20260708_rec095_block1/` folder is empty. There is no `.mwk2` or `.h5` on cnl04r or on the lab server, so the recording fails if it is added to `rec_lst`.
- **33-channel sessions.** rec077, rec089, rec090, rec091 and rec096 block 2 were recorded with one 32-channel V-probe in MT (channels 1–32). Channel 33 is a headstage input with no electrode, recorded as a noise control; the v11 analysis drops it. Their 33 caches are complete.
- **Channel list.** The channel list is taken from the local `mua_ch###.mat` files when any exist; otherwise it comes from the PL2. After an interrupted extraction, delete that session's `mua_ch*.mat` files to force a full recompute.
- **Time base.** The envelope time base assumes one continuous PL2 recording that starts at t = 0.
- **Juice volume.** The juice is the *commanded* volume. The pump calibration is not logged.

## Changes 2026-09-15

| File | Change |
|---|---|
| `PHY_preprocessing_v3` 2.4 | Catch-trial table `stim.catch_*`. Per channel: catch envelope and baseline windows, and the pooled reference (`baseline_ref`, `catch_baseline_ref`, `ref_k`). Existing fields are unchanged. |
| `fn_muae_reference` 1.0 | New: the pooled running reference |
| `fn_sort_MUAe_by_state` 1.2 | Normalises by the pooled reference by default. New options `'baseline'` and `'ref_k'`; provenance in `state.muae_baseline` and `state.muae_ref_k`. |
| `preflight_pooled_baseline` | New: the one-session check before a full run |
| `PHY_main_analysis_v4` | Unchanged. Its call `fn_sort_MUAe_by_state(…, 'norm', muae_norm)` now uses the pooled reference; the comment above `muae_norm` still says "per-cycle". |

### Tests

The code could not be run in MATLAB on the laptop on 2026-09-15, because the licence server was unreachable. The Python comparison that chose K used the same definition. `preflight_pooled_baseline` is the MATLAB test.

## Changes 2026-09-10/11

| File | Change |
|---|---|
| `PHY_preprocessing_v3` 2.2 | Partner outcome and reward, reward scale, behaviour-only mode, `.mwk2` name fallback (rec081), time-based per-target outcome, juice and overlap |
| `PHY_preprocessing_v3` 2.3 | **Speed.** Each variable is extracted from the event table once; cycle windows are found by binary search instead of about 20 whole-table masks per cycle. MUAe cycle slices also use binary search. RF-mapping trials are parsed once per session instead of once per channel. `cpr_cyle` is built in one place. **Fixes.** Joystick direction is taken at the strength sample indices. A split direction/strength pair at a cycle edge had made direction one sample short and crashed rec077 in `sort_states`. The pump variable is found under both names, and a failed `.h5` write no longer stops the run. |
| `fn_RF_responses_muae` 1.1 | Returns and reuses the parsed presentation table across channels |
| `fn_sort_MUAe_by_state` 1.1 | State windows found by binary search |
| `PHY_main_analysis_v4` 4.2 / 4.3 | 29 sessions, `rec_behav_only` list, behaviour-only fallback, logs with minutes per recording, new state fields. `sort_states` tolerates the joystick length mismatch in older summaries. |
| `felix_nhp_solo.cfg` | Adds `CTRL_reward_target_ml`, `CTRL_reward_target_min_ml`, `CTRL_reward_target_max_ml`, `CTRL_reward2_target_ml`, `CTRL_reward_power`, `CTRL_arc_flag`, `CTRL_arc2_flag` and `AGNT_arc_flag`, so that the `.h5` keeps them |

### Tests

The tests ran in local MATLAB R2021a on rec058, rec077 and rec094. Old (2.2) and new (2.3) output were compared with `isequaln`.

- **Behaviour-only.** Identical on rec058 and rec094.
- **Sorted + MUAe with synthetic neural inputs.** These were sync parameters, a spike channel and 2–8 cached envelope channels. Output was identical, including the RF fits, `sort_states` and `fn_sort_MUAe_by_state`.
- **rec077.** Output was identical except for the one split joystick pair in cycle 75. The old `sort_states` reproduces the crash of 2026-09-10; the new one runs on both the old and the new output (1,422 states).
- **`.h5` round trip.** With the extended config, all 8 `CTRL_` / `AGNT_arc` variables survive `MW_writeH5` → `MW_readData`. The old config drops them.
- **Speed.** Per MUAe channel, 5.6 s → 0.22 s (64 channels: 6 min → 15 s per session). Behaviour stage: 14 s → 2 s for rec094.
