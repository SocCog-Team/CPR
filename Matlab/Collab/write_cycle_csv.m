function write_cycle_csv(phy, iCycle, output_dir)
% WRITE_CYCLE_CSV  Export one CPR experimental cycle to a frame-aligned CSV.
%
%   WRITE_CYCLE_CSV(PHY, ICYCLE, OUTPUT_DIR) reads raw stimulus and joystick
%   data for cycle ICYCLE from the PHY struct, resamples everything onto the
%   screen frame timestamps, and writes a CSV to OUTPUT_DIR.
%
% -------------------------------------------------------------------------
%   OUTPUT FILE
%       <rec_num>_<block>_cycle<NN>.csv  (filename matches header JSON)
%
%   COLUMNS
%       frame_idx           0-based integer frame counter
%       time_s              Time from cycle onset (seconds)
%       stim_direction_deg  RDP direction [0, 360) — step-held, not interpolated
%       stim_coherence      RDP coherence [0, 1]   — step-held, not interpolated
%       resp_dir_monk_deg   Monkey joystick direction [0, 360)
%       resp_conf_monk      Monkey joystick tilt [0, 1]
%       resp_dir_hum_deg    Human  joystick direction [0, 360)  [NaN in solo]
%       resp_conf_hum       Human  joystick tilt [0, 1]         [NaN in solo]
%
% -------------------------------------------------------------------------
%   RESAMPLING METHODS
%       Stimulus            Piecewise-constant (step at each change event).
%                           Values are never altered between change points.
%
%       Joystick direction  Circular linear interpolation:
%                             1. Convert degrees -> radians
%                             2. Decompose into sin / cos components
%                             3. Linearly interpolate each component
%                             4. Recover angle with atan2, convert to [0,360)
%
%       Joystick tilt       Scalar linear interpolation (interp1).
%
%       Edge behaviour      Both joystick signals are clamped to the value of
%                           the first / last available sample when the frame
%                           grid extends beyond the joystick timestamps.
%                           Change this by editing the INTERP OPTIONS section.
%
% -------------------------------------------------------------------------
%   SOLO vs DYAD
%       In solo cycles (cpr_solo == true) the absent partner's columns are
%       filled with NaN. Condition is read from phy.stim.cpr_solo{iCycle} and
%       phy.stim.cpr_dyad{iCycle}.
%
% -------------------------------------------------------------------------
%   REQUIREMENTS
%       MATLAB R2019b or later (writematrix / writecell).
%
% -------------------------------------------------------------------------
%   EXAMPLE
%       for iCycle = 1:203
%           write_cycle_csv(phy, iCycle, fullfile('data', phy.exp.rec_num));
%       end
%
%   See also: write_experiment_headers, interp1, atan2

% =========================================================================
% 0.  SETUP & VALIDATION
% =========================================================================

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Convenience aliases
e   = phy.exp;
stm = phy.stim;
joy = phy.joy;

is_solo = stm.cpr_solo{iCycle};
is_dyad = stm.cpr_dyad{iCycle};

% Warn if condition is ambiguous
if is_solo && is_dyad
    warning('write_cycle_csv:ambiguousCondition', ...
        'Cycle %d: cpr_solo and cpr_dyad are both true. Treating as dyad.', iCycle);
    is_solo = false;
end

% =========================================================================
% 1.  FRAME GRID
% =========================================================================

frme_ts   = double(stm.frme_ts{iCycle});    % [1 x nFrames] absolute timestamps (s)
cycle_ts  = double(stm.cpr_cyle(iCycle,:)); % [t_start, t_end]

t_start   = cycle_ts(1);
nFrames   = numel(frme_ts);

frame_idx = (0 : nFrames-1)';               % 0-based index column
time_s    = (frme_ts - t_start)';           % time relative to cycle onset (s)

% =========================================================================
% 2.  STIMULUS COLUMNS  (piecewise-constant step expansion)
%     Stimulus values are held constant between change events.
%     We do NOT interpolate stimulus — only joystick is resampled.
% =========================================================================

stim_dir = expand_stepwise( ...
    double(stm.rdp_dir{iCycle}), ...
    double(stm.rdp_dir_ts{iCycle}), ...
    frme_ts);                               % [nFrames x 1]

stim_coh = expand_stepwise( ...
    double(stm.rdp_coh{iCycle}), ...
    double(stm.rdp_coh_ts{iCycle}), ...
    frme_ts);                               % [nFrames x 1]

% =========================================================================
% 3.  JOYSTICK RESAMPLING
% =========================================================================

% --- 3a. Monkey ----------------------------------------------------------
if ~isempty(joy.js_monk_ts{iCycle})
    monk_dir = interp_circular( ...
        double(joy.js_monk_dir{iCycle}), ...
        double(joy.js_monk_ts{iCycle}), ...
        frme_ts);                           % circular interp, degrees [0,360)

    monk_tlt = interp_scalar( ...
        double(joy.js_monk_tlt{iCycle}), ...
        double(joy.js_monk_ts{iCycle}), ...
        frme_ts);                           % linear interp, clamped at edges
else
    monk_dir = nan(nFrames, 1);
    monk_tlt = nan(nFrames, 1);
end

% --- 3b. Human -----------------------------------------------------------
if is_dyad && ~isempty(joy.js_hum_ts{iCycle})
    hum_dir = interp_circular( ...
        double(joy.js_hum_dir{iCycle}), ...
        double(joy.js_hum_ts{iCycle}), ...
        frme_ts);

    hum_tlt = interp_scalar( ...
        double(joy.js_hum_tlt{iCycle}), ...
        double(joy.js_hum_ts{iCycle}), ...
        frme_ts);
else
    % Solo cycle or missing data: absent partner -> NaN
    hum_dir = nan(nFrames, 1);
    hum_tlt = nan(nFrames, 1);
end

% =========================================================================
% 4.  FEEDBACK COLUMNS  (sparse — NaN on non-event frames)
%
%   feedback_outcome : 1 = hit, 0 = miss, NaN = no feedback on this frame
%   reward_cum       : cumulative reward at the feedback frame, NaN elsewhere
%
%   Each feedback timestamp is snapped to the nearest frame index.
% =========================================================================

feedback_ts  = double(stm.feedback_ts{iCycle});   % [1 x nEvents] absolute timestamps
outcome_cell = stm.outcome{iCycle};               % {1 x nEvents} cell of 'hit'/'miss'
reward_cv    = double(stm.reward_cum{iCycle});     % [1 x nEvents] cumulative reward

[fb_outcome, fb_reward] = expand_feedback( ...
    feedback_ts, outcome_cell, reward_cv, frme_ts);

% =========================================================================
% 5.  ASSEMBLE TABLE AND WRITE CSV
% =========================================================================

% Column headers (must match README / header JSON schema)
col_names = { ...
    'frame_idx', ...
    'time_s', ...
    'stim_direction_deg', ...
    'stim_coherence', ...
    'resp_dir_monk_deg', ...
    'resp_conf_monk', ...
    'resp_dir_hum_deg', ...
    'resp_conf_hum', ...
    'feedback_outcome', ...
    'reward_cum'};

% Numeric matrix: one row per frame
data_matrix = [ ...
    frame_idx, ...
    time_s, ...
    stim_dir, ...
    stim_coh, ...
    monk_dir, ...
    monk_tlt, ...
    hum_dir, ...
    hum_tlt, ...
    fb_outcome, ...
    fb_reward ];

% Build output filename to match the naming in the header JSON
out_name = sprintf('%s_%s_cycle%02d.csv', e.rec_num, e.block, iCycle);
out_path = fullfile(output_dir, out_name);

% Write header row then data
fid = fopen(out_path, 'w');
if fid == -1
    error('write_cycle_csv:fileOpenFailed', ...
        'Could not open file for writing: %s', out_path);
end
fprintf(fid, '%s\n', strjoin(col_names, ','));
fclose(fid);

% Append numeric data (writematrix appends when file exists)
writematrix(data_matrix, out_path, 'WriteMode', 'append');

fprintf('[write_cycle_csv] Cycle %03d written -> %s  (%d frames)\n', ...
    iCycle, out_name, nFrames);

end % main function


% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function vals_at_frames = expand_stepwise(values, change_ts, frame_ts)
% EXPAND_STEPWISE  Map event-based (value, timestamp) pairs onto a frame grid.
%
%   For each frame, the active value is the last one whose change timestamp
%   is <= the frame timestamp (piecewise-constant / zero-order hold).
%
%   INPUTS
%       values      [1 x nEvents] array of stimulus values
%       change_ts   [1 x nEvents] absolute timestamps of each value onset
%       frame_ts    [1 x nFrames] absolute frame timestamps
%
%   OUTPUT
%       vals_at_frames  [nFrames x 1] column vector

nFrames        = numel(frame_ts);
vals_at_frames = nan(nFrames, 1);

for iFrame = 1:nFrames
    % Find all change events that have occurred by this frame
    active_idx = find(change_ts <= frame_ts(iFrame), 1, 'last');
    if ~isempty(active_idx)
        vals_at_frames(iFrame) = values(active_idx);
    end
    % If no event yet, value stays NaN (frame precedes first stimulus onset)
end

end % expand_stepwise


% -------------------------------------------------------------------------

function [outcome_at_frames, reward_at_frames] = expand_feedback( ...
    feedback_ts, outcome_cell, reward_cum, frame_ts)
% EXPAND_FEEDBACK  Map sparse feedback events onto the frame grid.
%
%   Each feedback event is assigned to the nearest frame in time.
%   All other frames receive NaN.
%
%   Outcome encoding:
%       1   = hit
%       0   = miss
%       NaN = no feedback event on this frame
%
%   INPUTS
%       feedback_ts     [1 x nEvents] absolute timestamps of feedback events
%       outcome_cell    {1 x nEvents} cell array of 'hit' / 'miss' strings
%       reward_cum      [1 x nEvents] cumulative reward at each event
%       frame_ts        [1 x nFrames] absolute frame timestamps
%
%   OUTPUT
%       outcome_at_frames   [nFrames x 1] 1 / 0 / NaN
%       reward_at_frames    [nFrames x 1] cumulative reward / NaN

nFrames           = numel(frame_ts);
outcome_at_frames = nan(nFrames, 1);
reward_at_frames  = nan(nFrames, 1);

for iEvent = 1:numel(feedback_ts)

    % Snap to nearest frame
    [~, nearest_frame] = min(abs(frame_ts - feedback_ts(iEvent)));

    % Encode outcome as numeric (1 = hit, 0 = miss)
    if strcmpi(outcome_cell{iEvent}, 'hit')
        outcome_at_frames(nearest_frame) = 1;
    else
        outcome_at_frames(nearest_frame) = 0;
    end

    reward_at_frames(nearest_frame) = reward_cum(iEvent);
end

end % expand_feedback


% -------------------------------------------------------------------------

function dir_at_frames = interp_circular(dir_deg, joy_ts, frame_ts)
% INTERP_CIRCULAR  Resample angular joystick direction onto frame timestamps.
%
%   Uses unit-circle decomposition to avoid wrap-around artefacts:
%       1. Convert degrees to radians
%       2. Compute sin and cos components
%       3. Linearly interpolate each component independently
%       4. Recover angle with atan2 and convert back to [0, 360) degrees
%
%   Edge behaviour: clamped (first/last available sample is held).
%
%   INPUTS
%       dir_deg     [1 x N] joystick direction samples in degrees
%       joy_ts      [1 x N] joystick sample timestamps (seconds)
%       frame_ts    [1 x F] frame timestamps (seconds)
%
%   OUTPUT
%       dir_at_frames   [F x 1] interpolated direction in degrees [0, 360)

% Convert to radians and decompose
dir_rad = deg2rad(dir_deg);
s = sin(dir_rad);
c = cos(dir_rad);

% Linear interpolation with edge clamping
% 'extrap' is intentionally omitted so MATLAB clamps to boundary values
s_interp = interp1(joy_ts(:), s(:), frame_ts(:), 'linear', 'extrap');
c_interp = interp1(joy_ts(:), c(:), frame_ts(:), 'linear', 'extrap');

% --- INTERP OPTIONS ---
% To use nearest-neighbour instead:
%   s_interp = interp1(joy_ts(:), s(:), frame_ts(:), 'nearest', 'extrap');
%   c_interp = interp1(joy_ts(:), c(:), frame_ts(:), 'nearest', 'extrap');
% To fill out-of-range frames with NaN instead of clamping:
%   s_interp = interp1(joy_ts(:), s(:), frame_ts(:), 'linear');
%   c_interp = interp1(joy_ts(:), c(:), frame_ts(:), 'linear');
% ----------------------

% Recover angle and convert to [0, 360)
dir_at_frames = mod(rad2deg(atan2(s_interp, c_interp)), 360);

end % interp_circular


% -------------------------------------------------------------------------

function tlt_at_frames = interp_scalar(tlt, joy_ts, frame_ts)
% INTERP_SCALAR  Resample scalar joystick tilt onto frame timestamps.
%
%   Standard linear interpolation. Tilt is a non-circular magnitude so no
%   unit-circle decomposition is needed.
%
%   Edge behaviour: clamped (first/last available sample is held).
%
%   INPUTS
%       tlt         [1 x N] joystick tilt samples
%       joy_ts      [1 x N] joystick sample timestamps (seconds)
%       frame_ts    [1 x F] frame timestamps (seconds)
%
%   OUTPUT
%       tlt_at_frames   [F x 1] interpolated tilt values

tlt_at_frames = interp1(joy_ts(:), tlt(:), frame_ts(:), 'linear', 'extrap');

% --- INTERP OPTIONS ---
% To fill out-of-range frames with NaN instead of clamping:
%   tlt_at_frames = interp1(joy_ts(:), tlt(:), frame_ts(:), 'linear');
% ----------------------

end % interp_scalar
