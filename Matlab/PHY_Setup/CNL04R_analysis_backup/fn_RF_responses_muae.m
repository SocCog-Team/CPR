function [pos, resp, base, is_target, stim_id, x_dva, y_dva] = ...
        fn_RF_responses_muae(d, env_t_us, env, win_offset, target_pos)
% FN_RF_RESPONSES_MUAE  Per-presentation MUAe activity for the RF-mapping task.
%
% Mirrors the sorted-spike RF_mapping trial parsing, but the activity in each
% window is the MEAN ENVELOPE instead of a spike count.  Output feeds
% fn_characterise_RF (response_mode='percent').
%
% Uses the CONTINUOUS envelope (which spans the whole recording, including the
% RF-mapping trials that precede the CPR block), sliced by uniform-grid index
% for speed.
%
% INPUT
%   d           struct   MWorks event/time/value struct
%   env_t_us    [1xM]    envelope sample times, MWorks µs (uniform spacing)
%   env         [1xM]    MUAe envelope amplitude (single/double)
%   win_offset  [1x2]    response window [start end] µs after stimulus onset
%                        (default [50e3 150e3], matching RF_mapping)
%   target_pos  char     fixation target label (default 'X0_Y0')
%
% OUTPUT  (one entry per RF-mapping stimulus presentation)
%   pos        [1xP]   grid stimulus id (index into stim_id)
%   resp       [1xP]   mean envelope in the response window
%   base       [1xP]   mean envelope in the pre-stimulus baseline window
%   is_target  [1xP]   logical, true for the [0,0] presentations
%   stim_id    [RxC]   grid position -> stimulus id map
%   x_dva,y_dva        grid coordinates (dva)
%
% Felix Schneider, CNL  (drafted with Claude)
%
% Version history
%   1.0  (2026-07-03)  Initial version.

if nargin < 4 || isempty(win_offset); win_offset = [50e3 150e3]; end
if nargin < 5 || isempty(target_pos); target_pos = 'X0_Y0';      end

% --- grid (must match RF_mapping) ---------------------------------------
x_dva   = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
y_dva   = fliplr([-12 -9 -6 -3 0 3 6 9 12]);
stim_id = reshape(1:(numel(x_dva)*numel(y_dva)), [numel(y_dva), numel(x_dva)]);

% --- event indices ------------------------------------------------------
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

% --- envelope uniform-grid parameters (for fast windowed means) ---------
env   = double(env);
t0_us = double(env_t_us(1));
dt_us = double(env_t_us(2) - env_t_us(1));

% --- iterate RF-mapping trials ------------------------------------------
maxP      = 20000;
pos       = nan(1, maxP);
resp      = nan(1, maxP);
base      = nan(1, maxP);
is_target = false(1, maxP);
cnt       = 0;

for iTrl = 1:numel(trl.start)

    if ~strcmp(trl.task{iTrl}, 'RF_mapping')
        continue
    end

    trlIdx   = d.time >= trl.start(iTrl) & d.time <= trl.end(iTrl);
    outcome  = getTrialData(d.value, trlIdx, idx.outcome);
    ttype    = getTrialData(d.value, trlIdx, idx.type);
    rdp_time = getTrialData(d.time,  trlIdx, idx.type);

    if ~iscell(ttype)
        continue
    end

    if strcmp(outcome, 'fixation break')
        nStim = numel(ttype) - 1;
    else
        nStim = numel(ttype);
    end

    % Baseline: trial start -> first RDP onset (stimulus-free).
    bl_val = win_mean_env(env, t0_us, dt_us, trl.start(iTrl), rdp_time(1));

    for iStim = 1:nStim
        p_split = split(ttype{iStim}, '_');
        rdp_x   = str2double(p_split{1}(2:end));
        rdp_y   = str2double(p_split{2}(2:end));

        cnt            = cnt + 1;
        pos(cnt)       = stim_id(rdp_y == y_dva, rdp_x == x_dva);
        base(cnt)      = bl_val;
        resp(cnt)      = win_mean_env(env, t0_us, dt_us, ...
                                      rdp_time(iStim) + win_offset(1), ...
                                      rdp_time(iStim) + win_offset(2));
        is_target(cnt) = strcmp(ttype{iStim}, target_pos);
    end
end

pos       = pos(1:cnt);
resp      = resp(1:cnt);
base      = base(1:cnt);
is_target = is_target(1:cnt);

end % fn_RF_responses_muae


function mu = win_mean_env(env, t0_us, dt_us, ta, tb)
% Mean envelope between MWorks times ta and tb (µs), via uniform-grid index.
N  = numel(env);
i1 = min(max(round((double(ta) - t0_us) / dt_us) + 1, 1), N);
i2 = min(max(round((double(tb) - t0_us) / dt_us) + 1, 1), N);
if i2 < i1
    mu = NaN;
else
    mu = mean(env(i1:i2), 'omitnan');
end
end % win_mean_env
