function mua = fn_compute_MUA(sig, fs, varargin)
% FN_COMPUTE_MUA  Multi-unit activity from a single wideband/filtered channel.
%
% Computes two complementary drift-robust multi-unit measures from one
% channel of continuous data, following Stark & Abeles (2007):
%
%   MUAe  — analog multi-unit ENVELOPE.  The signal is (optionally)
%           band-passed into the MUA band, full-wave rectified, low-pass
%           filtered, and downsampled.  A graded, continuous estimate of
%           pooled spiking power around the electrode.  Robust to drift
%           because it needs no unit isolation.
%
%   MUAt  — multi-unit THRESHOLD crossings ("hash").  Threshold crossings
%           on the band-passed signal, relative to a robust noise estimate.
%           Yields spike-like event times that drop straight into a
%           spike-timestamp pipeline with no downstream changes.
%
% This function is deliberately I/O-free and works on a single channel so
% it can be unit-tested in isolation.  Time base of all returned times is
% the LOCAL data clock (t = 0 at sig(1), seconds).  Synchronisation to
% MWorks/behaviour is applied by the caller (see fn_extract_MUA_session).
%
% INPUT
%   sig   [1xN] or [Nx1]   one channel of continuous data (int16/single/double)
%                          - if already band-passed (e.g. datafilt, 333-5000 Hz)
%                            leave 'rebandpass' = false (default)
%                          - if raw broadband, set 'rebandpass' = true to
%                            apply muae_band first
%   fs    scalar           sampling rate of sig, Hz (e.g. 40000)
%
% OPTIONAL NAME-VALUE PAIRS (defaults follow Stark & Abeles 2007 / lab pipeline)
%   'rebandpass'  logical  false  apply muae_band band-pass before rectifying
%                                 (set true when sig is raw broadband)
%   'muae_band'   [lo hi]  [300 5000]  MUA band-pass (Hz), used if rebandpass
%   'rect'        char     'square'    rectification: 'square' (RMS envelope,
%                                 matches lab MUA code) or 'abs' (mean rectified)
%   'env_lp_Hz'   scalar   200    envelope low-pass cutoff (Hz)
%   'fs_out'      scalar   1000   MUAe output sampling rate (Hz)
%   'filt_order'  scalar   4      Butterworth order for band-pass / env LP
%   'do_muat'     logical  true   also detect threshold crossings (MUAt)
%   'thr_factor'  scalar   4      MUAt threshold = thr_factor * robust sigma
%   'polarity'    char     'neg'  'neg' | 'pos' | 'both' crossing direction
%                                 ('neg' matches the negthr spike-sorting path)
%   'refractory_ms' scalar 0.5    dead time enforced between MUAt events (ms)
%   'sigma_win_s' scalar   []     if set, estimate MUAt noise in sliding
%                                 windows of this length (s) so the threshold
%                                 tracks slow drift; [] = single global sigma
%
% OUTPUT  mua  struct
%   .env          [1xM]   MUAe envelope at fs_out (same units as sig; sqrt-ed
%                         back to amplitude when rect='square')
%   .env_fs       scalar  fs_out
%   .env_t_s      [1xM]   envelope sample times, local clock (s), t0 = sig(1)
%   .muat_idx     [1xK]   MUAt crossing sample indices into sig (at fs)
%   .muat_t_s     [1xK]   MUAt crossing times, local clock (s)
%   .sigma        scalar or [1xK]  robust noise estimate used for thresholding
%   .thr          scalar or [1xK]  applied threshold(s)
%   .params       struct  all parameters + provenance
%
% NOTES
%   * MUAe amplitude is arbitrary and drifts with the electrode; normalise
%     downstream (this pipeline uses a per-cycle 300 ms pre-stimulus baseline).
%   * 'square'+sqrt gives an RMS-like envelope (matches the lab's existing
%     Stark-Abeles code); 'abs' gives the classic mean-rectified envelope.
%   * Robust sigma = median(|x|)/0.6745 (Quiroga et al. 2004) is insensitive
%     to the spikes themselves, unlike std().
%
% Felix Schneider, CNL  (drafted with Claude)
%
% Version history
%   1.0  (2026-07-03)  Initial version. MUAe (Stark & Abeles 2007) + MUAt.

% =========================================================================
%% Parse inputs
% =========================================================================

p = inputParser;
p.addRequired('sig',  @(x) isnumeric(x) && isvector(x));
p.addRequired('fs',   @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('rebandpass',    false,      @(x) islogical(x) && isscalar(x));
p.addParameter('muae_band',     [300 5000], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('rect',          'square',   @(x) any(strcmpi(x, {'square','abs'})));
p.addParameter('env_lp_Hz',     200,        @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('fs_out',        1000,       @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('filt_order',    4,          @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('do_muat',       true,       @(x) islogical(x) && isscalar(x));
p.addParameter('thr_factor',    4,          @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('polarity',      'neg',      @(x) any(strcmpi(x, {'neg','pos','both'})));
p.addParameter('refractory_ms', 0.5,        @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('sigma_win_s',   [],         @(x) isempty(x) || (isscalar(x) && x > 0));
p.parse(sig, fs, varargin{:});
o = p.Results;

% Work in double, row vector, for numerical stability of filtfilt.
x = double(sig(:)');
N = numel(x);

nyq = fs / 2;
if o.env_lp_Hz >= nyq
    error('fn_compute_MUA:badCutoff', 'env_lp_Hz (%g) must be < Nyquist (%g).', o.env_lp_Hz, nyq);
end
if o.fs_out > fs
    error('fn_compute_MUA:badFsOut', 'fs_out (%g) cannot exceed fs (%g).', o.fs_out, fs);
end

% =========================================================================
%% 1. (Optional) band-pass into the MUA band
%     datafilt is already 333-5000 Hz, so this is skipped by default.
% =========================================================================

if o.rebandpass
    if o.muae_band(2) >= nyq
        o.muae_band(2) = nyq * 0.99;   % guard against invalid Wn
    end
    [b_bp, a_bp] = butter(o.filt_order, o.muae_band / nyq, 'bandpass');
    x = filtfilt(b_bp, a_bp, x);
end

% =========================================================================
%% 2. MUAt — threshold crossings on the band-passed signal
%     Done BEFORE rectification, on the (bandpassed) spike-band signal.
% =========================================================================

mua.muat_idx = [];
mua.muat_t_s = [];
mua.sigma    = NaN;
mua.thr      = NaN;

if o.do_muat
    refr_samples = max(1, round(o.refractory_ms / 1000 * fs));

    if isempty(o.sigma_win_s)
        % Single global robust noise estimate.
        sigma = median(abs(x)) / 0.6745;
        thr   = o.thr_factor * sigma;
        [idx, thr_used] = local_detect_crossings(x, thr, o.polarity, refr_samples);
        mua.sigma = sigma;
        mua.thr   = thr_used;
    else
        % Drift-tracking: per-window sigma, threshold follows slow changes.
        win = max(1, round(o.sigma_win_s * fs));
        edges = 1 : win : N;
        sigma_vec = nan(1, numel(edges));
        idx = [];
        for k = 1 : numel(edges)
            a = edges(k);
            b = min(a + win - 1, N);
            seg = x(a:b);
            s   = median(abs(seg)) / 0.6745;
            sigma_vec(k) = s;
            [seg_idx, ~] = local_detect_crossings(seg, o.thr_factor * s, o.polarity, refr_samples);
            idx = [idx, seg_idx + (a - 1)]; %#ok<AGROW>
        end
        % De-duplicate crossings straddling window edges (enforce refractory).
        idx = sort(idx);
        idx = local_enforce_refractory(idx, refr_samples);
        mua.sigma = sigma_vec;
        mua.thr   = o.thr_factor * sigma_vec;
    end

    mua.muat_idx = idx;
    mua.muat_t_s = (idx - 1) / fs;   % local clock, t0 at sig(1)
end

% =========================================================================
%% 3. MUAe — rectify, low-pass, downsample  (Stark & Abeles 2007)
% =========================================================================

switch lower(o.rect)
    case 'square'
        r = x .^ 2;
    case 'abs'
        r = abs(x);
end

% Envelope low-pass (zero-phase).
[b_lp, a_lp] = butter(o.filt_order, o.env_lp_Hz / nyq, 'low');
env = filtfilt(b_lp, a_lp, r);

% Downsample to fs_out. env is already band-limited to env_lp_Hz << fs_out/2,
% so resampling introduces no aliasing.
[pp, qq] = rat(o.fs_out / fs);
env = resample(env, pp, qq);
fs_out_actual = fs * pp / qq;

% Undo the squaring to return to amplitude units (RMS-like envelope).
if strcmpi(o.rect, 'square')
    env(env < 0) = 0;      % clamp tiny negative ripples from resampling
    env = sqrt(env);
end

mua.env      = env;
mua.env_fs   = fs_out_actual;
mua.env_t_s  = (0 : numel(env) - 1) / fs_out_actual;   % local clock

% =========================================================================
%% 4. Provenance
% =========================================================================

mua.params            = o;
mua.params.fs_in      = fs;
mua.params.n_samples  = N;
mua.params.fs_out     = fs_out_actual;
mua.params.created    = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
mua.params.func       = mfilename;

end % fn_compute_MUA


% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================

function [idx, thr_used] = local_detect_crossings(x, thr, polarity, refr_samples)
% Detect falling/rising threshold crossings, then enforce a dead time.
switch lower(polarity)
    case 'neg'
        below = x < -thr;
        cross = find(~below(1:end-1) & below(2:end)) + 1;   % downward through -thr
        thr_used = -thr;
    case 'pos'
        above = x > thr;
        cross = find(~above(1:end-1) & above(2:end)) + 1;   % upward through +thr
        thr_used = thr;
    case 'both'
        out   = abs(x) > thr;
        cross = find(~out(1:end-1) & out(2:end)) + 1;
        thr_used = thr;
end
idx = local_enforce_refractory(cross, refr_samples);
end % local_detect_crossings


function idx = local_enforce_refractory(cross, refr_samples)
% Greedy dead-time: keep a crossing only if it is >= refr_samples after the
% last KEPT crossing (correct, unlike a plain diff() mask).
if isempty(cross)
    idx = [];
    return
end
cross = sort(cross(:)');
keep  = false(1, numel(cross));
last  = -Inf;
for i = 1 : numel(cross)
    if cross(i) - last >= refr_samples
        keep(i) = true;
        last    = cross(i);
    end
end
idx = cross(keep);
end % local_enforce_refractory
