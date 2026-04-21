function cpr_write_cycle_frames(phy, iCycle, output_dir)
% CPR_WRITE_CYCLE_FRAMES  Export one CPR cycle to a frame-aligned CSV file.
%
%   CPR_WRITE_CYCLE_FRAMES(PHY, ICYCLE, OUTPUT_DIR) resamples all stimulus
%   and behavioural signals onto the display frame grid for cycle ICYCLE
%   and writes the result as a CSV to:
%
%       <output_dir>/cycle_frames/<rec>_<block>_cycle<NNN>_frames.csv
%
%   Each row corresponds to one display frame (~120 Hz). Columns cover
%   stimulus state, joystick responses, and sparse feedback events.
%
% -------------------------------------------------------------------------
%   OUTPUT COLUMNS
%
%       frame_idx              0-based integer frame counter
%       time_s                 Time from cycle onset (seconds)
%
%       stim_direction_deg     RDP direction [0, 360) — piecewise-constant,
%                                never interpolated between change events
%       stim_coherence         RDP coherence [0, 1] — piecewise-constant
%
%       resp_dir_monkey_deg    Monkey joystick direction [0, 360)
%       resp_conf_monkey       Monkey joystick tilt magnitude [0, 1]
%       resp_dir_human_deg     Human joystick direction [dyad only; NaN in solo]
%       resp_conf_human        Human joystick tilt      [dyad only; NaN in solo]
%
%       feedback_outcome_code  1=hit  0=miss  -1=FixationBreak  NaN=no event
%       feedback_outcome_str   All outcomes pipe-joined, e.g. "hit|FixationBreak"
%                                Empty string on non-feedback frames
%       reward_ind             Reward earned at this specific event (cycle-relative).
%                                NaN on non-feedback frames and on FixationBreak
%                                events without an associated reward timestamp.
%       reward_cycle_cum       Running reward total within the cycle (resets at
%                                each cycle onset). NaN on non-feedback frames.
%
% -------------------------------------------------------------------------
%   RESAMPLING METHODS
%
%       Stimulus              Piecewise-constant zero-order hold: the value
%                             at each frame is the last value whose change
%                             timestamp precedes that frame. Not altered.
%
%       Joystick direction    Circular linear interpolation:
%                               1. Degrees -> radians
%                               2. Decompose into sin / cos components
%                               3. Linearly interpolate each component
%                               4. Recover angle via atan2 -> [0, 360)
%                             Avoids wrap-around artefacts at 0/360 boundary.
%
%       Joystick tilt         Scalar linear interpolation (interp1).
%
%       Edge behaviour        Both joystick channels clamp to the nearest
%                             available sample when the frame grid extends
%                             beyond the joystick timestamp range.
%                             Edit the EDGE BEHAVIOUR comments in
%                             interp_joystick_direction / interp_joystick_tilt
%                             to switch to NaN-padding instead.
%
% -------------------------------------------------------------------------
%   FEEDBACK / OUTCOME ENCODING
%
%       Each feedback event is snapped to its nearest display frame.
%       If two events land on the same frame, priority determines the code:
%           hit (1)  >  miss (0)  >  FixationBreak (-1)
%       All outcomes at a frame are listed in feedback_outcome_str, regardless
%       of priority, separated by '|'.
%
%       Cycles terminated by a fixation break contain a 'FixationBreak' event
%       as the final feedback entry. It may co-occur with a 'hit' or 'miss'
%       from the same reward target event — both are recorded.
%
% -------------------------------------------------------------------------
%   INPUTS
%       phy         Master experiment struct (see cpr_write_session_header).
%       iCycle      Cycle index (1-based integer).
%       output_dir  Root output folder path.
%
% -------------------------------------------------------------------------
%   REQUIREMENTS
%       MATLAB R2019b or later (writetable).
%
% -------------------------------------------------------------------------
%   USAGE
%       for iCycle = 1:203
%           cpr_write_cycle_frames(phy, iCycle, fullfile('export', phy.exp.rec_num))
%       end
%
%   See also: cpr_write_session_header, cpr_write_unit_spikes
%
% -------------------------------------------------------------------------
%   Created    : 2026-04-16
%   AI Source  : Generated with Claude (Anthropic) — https://claude.ai
%   Reference  : Schneider et al. (2025) eLife https://doi.org/10.7554/eLife.101021.2

% =========================================================================
% 0.  SETUP
% =========================================================================

frames_dir = fullfile(output_dir, 'cycle_frames');
if ~exist(output_dir,  'dir'), mkdir(output_dir);  end
if ~exist(frames_dir,  'dir'), mkdir(frames_dir);  end

e   = phy.exp;
stm = phy.stim;
joy = phy.joy;

is_solo = stm.cpr_solo{iCycle};
is_dyad = stm.cpr_dyad{iCycle};

if is_solo && is_dyad
    warning('cpr_write_cycle_frames:ambiguousCondition', ...
        'Cycle %d: cpr_solo and cpr_dyad are both true. Treating as dyad.', iCycle);
    is_solo = false;
end

% =========================================================================
% 1.  FRAME GRID
% =========================================================================

% frme_ts: absolute timestamps (s) for every display frame in this cycle
frme_ts  = double(stm.frme_ts{iCycle});     % [1 x nFrames]
t_start  = double(stm.cpr_cyle(iCycle, 1)); % absolute cycle onset (s)
nFrames  = numel(frme_ts);

frame_idx = (0 : nFrames-1)';               % 0-based frame counter [nFrames x 1]
time_s    = (frme_ts(:) - t_start);         % time relative to cycle onset [nFrames x 1]

% =========================================================================
% 2.  STIMULUS COLUMNS  (piecewise-constant — zero-order hold)
%
%   RDP direction and coherence change at discrete events. Between events
%   the value is constant. We map each change onto the frame grid by
%   finding the last event timestamp <= each frame timestamp.
%   Stimulus values are NEVER interpolated.
% =========================================================================

stim_dir = expand_stimulus_to_frames( ...
    double(stm.rdp_dir{iCycle}), ...
    double(stm.rdp_dir_ts{iCycle}), ...
    frme_ts);                               % [nFrames x 1]

stim_coh = expand_stimulus_to_frames( ...
    double(stm.rdp_coh{iCycle}), ...
    double(stm.rdp_coh_ts{iCycle}), ...
    frme_ts);                               % [nFrames x 1]

% =========================================================================
% 3.  JOYSTICK RESAMPLING
%
%   Joystick samples are asynchronous to display frames. Both direction and
%   tilt are resampled onto the frame grid. Missing (empty) channels are
%   filled with NaN — this covers solo cycles and any recording gaps.
% =========================================================================

% --- 3a. Monkey joystick -------------------------------------------------
if ~isempty(joy.js_monk_ts{iCycle})
    resp_dir_monkey = interp_joystick_direction( ...
        double(joy.js_monk_dir{iCycle}), ...
        double(joy.js_monk_ts{iCycle}), ...
        frme_ts);

    resp_conf_monkey = interp_joystick_tilt( ...
        double(joy.js_monk_tlt{iCycle}), ...
        double(joy.js_monk_ts{iCycle}), ...
        frme_ts);
else
    resp_dir_monkey  = nan(nFrames, 1);
    resp_conf_monkey = nan(nFrames, 1);
end

% --- 3b. Human joystick --------------------------------------------------
%   Present only in dyadic cycles. NaN-filled for solo cycles or when the
%   human joystick channel has no data.
if is_dyad && ~isempty(joy.js_hum_ts{iCycle})
    resp_dir_human = interp_joystick_direction( ...
        double(joy.js_hum_dir{iCycle}), ...
        double(joy.js_hum_ts{iCycle}), ...
        frme_ts);

    resp_conf_human = interp_joystick_tilt( ...
        double(joy.js_hum_tlt{iCycle}), ...
        double(joy.js_hum_ts{iCycle}), ...
        frme_ts);
else
    resp_dir_human  = nan(nFrames, 1);
    resp_conf_human = nan(nFrames, 1);
end

% =========================================================================
% 4.  FEEDBACK COLUMNS  (sparse — NaN / empty on non-event frames)
%
%   There are typically a handful of feedback events per cycle (~6).
%   Each is snapped to its nearest frame. When multiple events land on the
%   same frame, the priority code reflects the most informative outcome
%   (hit > miss > FixationBreak), while the string column preserves all.
%
%   Reward values are expressed relative to cycle onset, not the session.
%   reward_ind  = reward earned at this specific feedback event (difference
%                 between consecutive cumulative values within the cycle).
%   reward_cum  = running total within the cycle at each feedback event.
%
%   Outcome codes:
%       1  = hit           (reward target acquired)
%       0  = miss          (reward target missed)
%      -1  = FixationBreak (cycle ended early; fixation requirement violated)
%     NaN  = no feedback event on this frame
%
%   IMPORTANT — FixationBreak handling:
%   A fixation break terminates the cycle and is stored in outcome{iCycle}
%   but may NOT have a corresponding entry in feedback_ts{iCycle} (because
%   it is a cycle-end signal, not a reward-target event). If outcome_cell
%   contains more entries than feedback_ts, the extra outcomes (typically
%   a trailing FixationBreak) are snapped to the last frame.
% =========================================================================

feedback_ts  = double(stm.feedback_ts{iCycle});   % [1 x nTimestamps]
outcome_cell = stm.outcome{iCycle};               % {1 x nOutcomes} — may be longer

% Compute cycle-relative cumulative reward.
% reward_cum stores session-wide cumulative values. We subtract the first
% value of this cycle's own reward_cum as the baseline, because that value
% is recorded at (or just before) the first feedback event of the cycle.
% For a miss as the first event this equals the pre-cycle session value,
% giving reward_cv_cycle(1) = 0 and reward_ind = 0 for that miss.
% For an empty cycle (no feedback) both arrays stay empty and all columns
% remain NaN.
reward_cv_session = double(stm.reward_cum{iCycle});   % session-wide cumulative
if isempty(reward_cv_session)
    reward_cv_cycle = [];                             % no feedback events this cycle
else
    reward_cv_cycle = reward_cv_session - reward_cv_session(1); % cycle-relative
end

[fb_code, fb_str, fb_reward_cum, fb_reward_ind, fb_reward_session] = ...
    expand_feedback_to_frames(feedback_ts, outcome_cell, ...
        reward_cv_cycle, reward_cv_session, frme_ts);

% =========================================================================
% 5.  ASSEMBLE TABLE AND WRITE CSV
%
%   writetable is used instead of writematrix because the outcome string
%   column (feedback_outcome_str) requires mixed numeric/text output.
%   Column names in the table become the CSV header row automatically.
% =========================================================================

T = table( ...
    frame_idx,          ...
    time_s,             ...
    stim_dir,           ...
    stim_coh,           ...
    resp_dir_monkey,    ...
    resp_conf_monkey,   ...
    resp_dir_human,     ...
    resp_conf_human,    ...
    fb_code,            ...
    fb_str,             ...
    fb_reward_ind,      ...
    fb_reward_cum,      ...
    fb_reward_session,  ...
    'VariableNames', {  ...
        'frame_idx',             ...
        'time_s',                ...
        'stim_direction_deg',    ...
        'stim_coherence',        ...
        'resp_dir_monkey_deg',   ...
        'resp_conf_monkey',      ...
        'resp_dir_human_deg',    ...
        'resp_conf_human',       ...
        'feedback_outcome_code', ...
        'feedback_outcome_str',  ...
        'reward_ind',            ...
        'reward_cycle_cum',      ...
        'reward_session_cum'     ...
    });

out_name = sprintf('%s_%s_cycle%03d_frames.csv', e.rec_num, e.block, iCycle);
out_path = fullfile(frames_dir, out_name);

writetable(T, out_path);

fprintf('[cpr_write_cycle_frames] Cycle %03d -> %s  (%d frames)\n', ...
    iCycle, out_name, nFrames);

end % main function


% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function vals = expand_stimulus_to_frames(values, change_ts, frame_ts)
% EXPAND_STIMULUS_TO_FRAMES  Map event-driven stimulus onto the frame grid.
%
%   Implements zero-order hold (piecewise-constant): the active value at
%   each frame is the most recent value whose onset timestamp <= frame time.
%   Frames before the first event receive NaN.
%
%   INPUTS
%       values      [1 x nEvents] stimulus values (direction or coherence)
%       change_ts   [1 x nEvents] onset timestamps of each value (seconds)
%       frame_ts    [1 x nFrames] display frame timestamps (seconds)
%
%   OUTPUT
%       vals        [nFrames x 1] stimulus value at each frame

nFrames = numel(frame_ts);
vals    = nan(nFrames, 1);

for iFrame = 1:nFrames
    % Last event at or before this frame
    idx = find(change_ts <= frame_ts(iFrame), 1, 'last');
    if ~isempty(idx)
        vals(iFrame) = values(idx);
    end
end

end % expand_stimulus_to_frames


% -------------------------------------------------------------------------

function dir_at_frames = interp_joystick_direction(dir_deg, joy_ts, frame_ts)
% INTERP_JOYSTICK_DIRECTION  Resample joystick direction onto the frame grid.
%
%   Direction is a circular (angular) variable. Naive linear interpolation
%   across the 360/0 wraparound produces artefacts (e.g. sweeping through
%   180 degrees when the true movement is a small step across 360 -> 2 deg).
%   This is resolved by decomposing into sin/cos components, interpolating
%   each independently, then recovering the angle via atan2.
%
%   Steps:
%       1. Convert degrees to radians
%       2. Compute sin and cos components
%       3. Linearly interpolate each component to the frame grid
%       4. Recover angle: atan2(sin, cos), convert to [0, 360) degrees
%
%   INPUTS
%       dir_deg     [1 x N] joystick direction samples (degrees)
%       joy_ts      [1 x N] joystick sample timestamps (seconds)
%       frame_ts    [1 x F] display frame timestamps (seconds)
%
%   OUTPUT
%       dir_at_frames   [F x 1] direction at each frame in degrees [0, 360)

rad    = deg2rad(dir_deg(:));
s_comp = sin(rad);
c_comp = cos(rad);

% EDGE BEHAVIOUR: 'extrap' clamps to the boundary value outside the joystick
% timestamp range. Replace 'extrap' with NaN to pad edges with NaN instead:
%   s_interp = interp1(joy_ts(:), s_comp, frame_ts(:), 'linear');
s_interp = interp1(joy_ts(:), s_comp, frame_ts(:), 'linear', 'extrap');
c_interp = interp1(joy_ts(:), c_comp, frame_ts(:), 'linear', 'extrap');

% Recover angle in [0, 360)
dir_at_frames = mod(rad2deg(atan2(s_interp, c_interp)), 360);

end % interp_joystick_direction


% -------------------------------------------------------------------------

function tlt_at_frames = interp_joystick_tilt(tlt, joy_ts, frame_ts)
% INTERP_JOYSTICK_TILT  Resample joystick tilt onto the frame grid.
%
%   Tilt is a scalar confidence magnitude — standard linear interpolation
%   is appropriate. No circular decomposition needed.
%
%   INPUTS
%       tlt         [1 x N] joystick tilt samples [0, 1]
%       joy_ts      [1 x N] joystick sample timestamps (seconds)
%       frame_ts    [1 x F] display frame timestamps (seconds)
%
%   OUTPUT
%       tlt_at_frames   [F x 1] tilt at each frame

% EDGE BEHAVIOUR: 'extrap' clamps to boundary value. Replace with NaN to
% pad edges instead:
%   tlt_at_frames = interp1(joy_ts(:), tlt(:), frame_ts(:), 'linear');
tlt_at_frames = interp1(joy_ts(:), tlt(:), frame_ts(:), 'linear', 'extrap');

end % interp_joystick_tilt


% -------------------------------------------------------------------------

function [fb_code, fb_str, fb_reward_cum, fb_reward_ind, fb_reward_session] = ...
    expand_feedback_to_frames(feedback_ts, outcome_cell, ...
        reward_cv_cycle, reward_cv_session, frame_ts)
% EXPAND_FEEDBACK_TO_FRAMES  Map sparse feedback events onto the frame grid.
%
%   Each event is snapped to the nearest display frame. Non-event frames
%   receive NaN (numeric) or empty string (text).
%
%   When multiple events snap to the same frame:
%     - fb_code       reflects the highest-priority outcome
%     - fb_str        lists ALL outcomes, pipe-separated
%
%   Priority: hit (1)  >  miss (0)  >  FixationBreak (-1)
%
%   FIXATION BREAK HANDLING
%   A FixationBreak has no feedback timestamp (it terminates the cycle, it
%   is not a reward event). Any outcome without a matching timestamp is
%   silently placed on the last frame with NaN reward columns.
%
%   REWARD COLUMNS
%   reward_cv_cycle is already cycle-relative (session offset subtracted by
%   the caller). Individual reward per event is computed here as the
%   difference between consecutive cycle-relative cumulative values.
%
%   INPUTS
%       feedback_ts       [1 x nTS]      feedback event timestamps (seconds)
%       outcome_cell      {1 x nOut}     outcome strings — may be longer than nTS
%       reward_cv_cycle   [1 x nTS]      cycle-relative cumulative reward
%       reward_cv_session [1 x nTS]      session-wide cumulative reward (raw)
%       frame_ts          [1 x nFrames]  display frame timestamps (seconds)
%
%   OUTPUTS
%       fb_code           [nFrames x 1] double — priority-encoded outcome or NaN
%       fb_str            [nFrames x 1] string — pipe-joined outcomes or ""
%       fb_reward_cum     [nFrames x 1] double — cycle-cumulative reward or NaN
%       fb_reward_ind     [nFrames x 1] double — individual reward per event or NaN
%       fb_reward_session [nFrames x 1] double — session-cumulative reward or NaN

nFrames      = numel(frame_ts);
fb_code      = nan(nFrames, 1);
fb_reward_cum = nan(nFrames, 1);
fb_reward_ind     = nan(nFrames, 1);
fb_reward_session = nan(nFrames, 1);
fb_str            = strings(nFrames, 1);   % "" by default

nTimestamps = numel(feedback_ts);
nOutcomes   = numel(outcome_cell);


% -------------------------------------------------------------------------
% Individual reward per event = diff of cycle-relative cumulative values.
%   hit  -> positive increment (reward earned)
%   miss -> 0 increment (reward_cum is flat between consecutive misses)
%
% Guard: reward_cv_cycle may be empty (cycle has feedback timestamps but
% no reward data — data inconsistency) or shorter than nTimestamps. In
% either case, pad with zeros so every event gets a defined reward_ind.
% -------------------------------------------------------------------------
nRewards = numel(reward_cv_cycle);

if nRewards < nTimestamps
    if nRewards > 0
        warning('cpr_write_cycle_frames:rewardLengthMismatch', ...
            'reward_cum has %d entries but feedback_ts has %d. Padding missing reward values with 0.', nRewards, nTimestamps);
    end
    % Pad to match timestamp count; missing entries treated as flat (no reward).
    % Guard: reward_cv_cycle(end) would crash on an empty array, so branch
    % explicitly — pad with zeros when empty, last value otherwise.
    if nRewards == 0
        pad_value = 0;
    else
        pad_value = reward_cv_cycle(end);
    end
    reward_cv_cycle_padded = [reward_cv_cycle(:); ...
        repmat(pad_value, nTimestamps - nRewards, 1)];
else
    reward_cv_cycle_padded = reward_cv_cycle(:);
end

% diff against 0 baseline: first event gets its full cycle-relative value;
% subsequent events get the increment since the previous event.
reward_ind_events = diff([0; reward_cv_cycle_padded(1:nTimestamps)]);  % [nTimestamps x 1]

% -------------------------------------------------------------------------
% Pass 1: events that have a matching timestamp (indices 1..nTimestamps)
% -------------------------------------------------------------------------
for iEvent = 1:nTimestamps

    [~, nearest_frame] = min(abs(frame_ts - feedback_ts(iEvent)));

    outcome_str = char(outcome_cell{iEvent});
    priority    = get_outcome_priority(outcome_str);

    % Priority-encoded outcome code
    if isnan(fb_code(nearest_frame)) || priority > fb_code(nearest_frame)
        fb_code(nearest_frame) = priority;
    end

    % Pipe-joined outcome string (all outcomes preserved)
    if fb_str(nearest_frame) == ""
        fb_str(nearest_frame) = outcome_str;
    else
        fb_str(nearest_frame) = fb_str(nearest_frame) + "|" + outcome_str;
    end

    % Reward: use most-recent event value if multiple events share a frame
    fb_reward_cum(nearest_frame)     = reward_cv_cycle_padded(iEvent);
    fb_reward_ind(nearest_frame)     = reward_ind_events(iEvent);
    % Guard: reward_cv_session may be shorter than nTimestamps in the
    % mismatch case; leave as NaN rather than crashing.
    if iEvent <= numel(reward_cv_session)
        fb_reward_session(nearest_frame) = reward_cv_session(iEvent);
    end

end

% -------------------------------------------------------------------------
% Pass 2: FixationBreak (and any other outcome without a timestamp)
%   A FixationBreak terminates the cycle and has no feedback timestamp.
%   It is silently added to the last frame. Reward columns stay NaN since
%   no reward event is associated with a fixation break.
% -------------------------------------------------------------------------
if nOutcomes > nTimestamps
    last_frame = nFrames;   % snap to final frame of the cycle

    for iEvent = (nTimestamps + 1) : nOutcomes
        outcome_str = char(outcome_cell{iEvent});
        priority    = get_outcome_priority(outcome_str);

        if isnan(fb_code(last_frame)) || priority > fb_code(last_frame)
            fb_code(last_frame) = priority;
        end

        if fb_str(last_frame) == ""
            fb_str(last_frame) = outcome_str;
        else
            fb_str(last_frame) = fb_str(last_frame) + "|" + outcome_str;
        end
        % reward_cum and reward_ind remain NaN — no reward event occurred
    end
end

end % expand_feedback_to_frames




% -------------------------------------------------------------------------

function p = get_outcome_priority(outcome_str)
% GET_OUTCOME_PRIORITY  Map an outcome string to a numeric priority code.
%
%   Priority order: hit (1) > miss (0) > FixationBreak (-1)
%   Unknown outcomes are assigned -2 with a warning.

switch lower(strtrim(outcome_str))
    case 'hit',            p =  1;
    case 'miss',           p =  0;
    case 'fixationbreak',  p = -1;
    otherwise
        p = -2;
        warning('cpr_write_cycle_frames:unknownOutcome', ...
            'Unrecognised outcome string: "%s". Assigned code -2.', outcome_str);
end

end % get_outcome_priority
