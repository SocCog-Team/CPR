function [pos, resp, base, is_target, stim_id, x_dva, y_dva, P] = ...
        fn_RF_responses_muae(d, env_t_us, env, win_offset, target_pos, P)
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
% The trial parsing does not depend on the channel: the first call returns the
% presentation table P; pass it back on every further channel of the session
% so the event table is parsed only once.
%
% INPUT
%   d           struct   MWorks event/time/value struct (only read to build P)
%   env_t_us    [1xM]    envelope sample times, MWorks µs (uniform spacing)
%   env         [1xM]    MUAe envelope amplitude (single/double)
%   win_offset  [1x2]    response window [start end] µs after stimulus onset
%                        (default [50e3 150e3], matching RF_mapping)
%   target_pos  char     fixation target label (default 'X0_Y0')
%   P           struct   (optional) presentation table from a previous call on
%                        the same session ([] = parse d)
%
% OUTPUT  (one entry per RF-mapping stimulus presentation)
%   pos        [1xP]   grid stimulus id (index into stim_id)
%   resp       [1xP]   mean envelope in the response window
%   base       [1xP]   mean envelope in the pre-stimulus baseline window
%   is_target  [1xP]   logical, true for the [0,0] presentations
%   stim_id    [RxC]   grid position -> stimulus id map
%   x_dva,y_dva        grid coordinates (dva)
%   P          struct  presentation table (windows in MWorks µs)
%
% Felix Schneider, CNL  (drafted with Claude)
%
% Version history
%   1.0  (2026-07-03)  Initial version.
%   1.1  (2026-09-10)  Parse the RF trials once per session (P, returned and
%                      reused across channels) with binary search on
%                      per-variable event series, instead of masking the whole
%                      event table per trial and channel. Output unchanged.

if nargin < 4 || isempty(win_offset); win_offset = [50e3 150e3]; end
if nargin < 5 || isempty(target_pos); target_pos = 'X0_Y0';      end
if nargin < 6 || isempty(P);          P = parse_presentations(d, win_offset, target_pos); end

stim_id   = P.stim_id;
x_dva     = P.x_dva;
y_dva     = P.y_dva;
pos       = P.pos;
is_target = P.is_target;

% --- envelope uniform-grid parameters (for fast windowed means) ---------
env   = double(env);
t0_us = double(env_t_us(1));
dt_us = double(env_t_us(2) - env_t_us(1));

% Baseline once per trial (trial start -> first RDP onset, stimulus-free),
% response once per presentation.
bl = nan(1, numel(P.bl_a));
for k = 1:numel(P.bl_a)
    bl(k) = win_mean_env(env, t0_us, dt_us, P.bl_a(k), P.bl_b(k));
end
base = bl(P.pres_trl);
resp = nan(1, numel(P.pos));
for k = 1:numel(P.pos)
    resp(k) = win_mean_env(env, t0_us, dt_us, P.rs_a(k), P.rs_b(k));
end

end % fn_RF_responses_muae


function P = parse_presentations(d, win_offset, target_pos)
% RF-mapping presentations of one session: grid id, target flag, and the
% baseline (per trial) and response (per presentation) windows in MWorks µs.
% A trial's task is the first INFO_task logged after its start; a trial that
% ended in a fixation break loses its last presentation.

% --- grid (must match RF_mapping) ---------------------------------------
P.x_dva   = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
P.y_dva   = fliplr([-12 -9 -6 -3 0 3 6 9 12]);
P.stim_id = reshape(1:(numel(P.x_dva)*numel(P.y_dva)), [numel(P.y_dva), numel(P.x_dva)]);

% --- per-variable event series (d is time-sorted) -------------------------
t_start = double(d.time(d.event == 'TRIAL_start'));
t_end   = double(d.time(d.event == 'TRIAL_end'));
task    = event_series(d, d.event == 'INFO_task');
outc    = event_series(d, d.event == 'TRIAL_outcome');
ttyp    = event_series(d, d.event == 'TRIAL_type');

maxP        = 20000;
P.pos       = nan(1, maxP);
P.is_target = false(1, maxP);
P.pres_trl  = zeros(1, maxP);             % presentation -> row of bl_a / bl_b
P.rs_a      = nan(1, maxP);
P.rs_b      = nan(1, maxP);
P.bl_a      = nan(1, numel(t_start));
P.bl_b      = nan(1, numel(t_start));
cnt  = 0;
ntrl = 0;

for iTrl = 1:min(numel(t_start), numel(t_end))
    ts = t_start(iTrl);
    te = t_end(iTrl);

    k = count_le(task.t, ts) + 1;                        % first INFO_task after the start
    if k > numel(task.t) || ~strcmp(task.v{k}, 'RF_mapping')
        continue
    end

    [i0, i1] = win_range(ttyp.t, ts, te);
    if i1 < i0 || ~ischar(ttyp.v{i0})                    % no (text) stimulus in this trial
        continue
    end
    ttype    = ttyp.v(i0:i1);
    rdp_time = ttyp.t(i0:i1);
    [o0, o1] = win_range(outc.t, ts, te);
    outcome  = outc.v(o0:o1);

    if ~isempty(outcome) && all(strcmp(outcome, 'fixation break'))
        nStim = numel(ttype) - 1;
    else
        nStim = numel(ttype);
    end

    ntrl         = ntrl + 1;
    P.bl_a(ntrl) = ts;
    P.bl_b(ntrl) = rdp_time(1);

    for iStim = 1:nStim
        p_split = split(ttype{iStim}, '_');
        rdp_x   = str2double(p_split{1}(2:end));
        rdp_y   = str2double(p_split{2}(2:end));

        cnt              = cnt + 1;
        P.pos(cnt)       = P.stim_id(rdp_y == P.y_dva, rdp_x == P.x_dva);
        P.is_target(cnt) = strcmp(ttype{iStim}, target_pos);
        P.pres_trl(cnt)  = ntrl;
        P.rs_a(cnt)      = rdp_time(iStim) + win_offset(1);
        P.rs_b(cnt)      = rdp_time(iStim) + win_offset(2);
    end
end

P.pos       = P.pos(1:cnt);
P.is_target = P.is_target(1:cnt);
P.pres_trl  = P.pres_trl(1:cnt);
P.rs_a      = P.rs_a(1:cnt);
P.rs_b      = P.rs_b(1:cnt);
P.bl_a      = P.bl_a(1:ntrl);
P.bl_b      = P.bl_b(1:ntrl);
end % parse_presentations


function S = event_series(d, sel)
% Times (double µs, sorted) and raw values of one MWorks variable.
S.t = double(d.time(sel));
S.v = d.value(sel);
if ~issorted(S.t)
    [S.t, o] = sort(S.t);
    S.v      = S.v(o);
end
end % event_series


function [i0, i1] = win_range(t, a, b)
% Index range i0:i1 of sorted T with A <= T <= B (i1 < i0 if none).
i0 = count_lt(t, a) + 1;
i1 = count_le(t, b);
end % win_range


function c = count_le(t, x)
% Number of elements of sorted T that are <= X (binary search).
lo = 1;   hi = numel(t);   c = 0;
while lo <= hi
    mid = floor((lo + hi) / 2);
    if t(mid) <= x; c = mid; lo = mid + 1; else; hi = mid - 1; end
end
end % count_le


function c = count_lt(t, x)
% Number of elements of sorted T that are < X (binary search).
lo = 1;   hi = numel(t);   c = 0;
while lo <= hi
    mid = floor((lo + hi) / 2);
    if t(mid) < x; c = mid; lo = mid + 1; else; hi = mid - 1; end
end
end % count_lt


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
