function out = fn_characterise_RF(pos, resp, base, is_target, stim_id, x_dva, y_dva, params)
% FN_CHARACTERISE_RF  Signal-agnostic receptive-field characterisation.
%
% Same RF computation as the sorted-spike path, but driven by generic
% per-presentation ACTIVITY values, so it works for any input signal:
%   sorted spikes / MUAt -> activity = spike COUNT in the window
%   MUAe                 -> activity = MEAN ENVELOPE in the window
% Because the maths is identical, a sorted-vs-MUAe RF comparison is fair.
%
% Combines (as in the original RF_characterisation):
%   PART A  Fixation-position response  — Wilcoxon signed-rank of baseline vs
%           response activity across the [0,0] presentations.
%   PART B  RF mapping — mean baseline-corrected activity map over the grid,
%           2-D Gaussian fit (centre/size), RF area, distance to fixation, and
%           2-sigma-ellipse ∩ cursor-aperture overlap on a fine grid.
%
% INPUT  (one entry per RF-mapping stimulus presentation)
%   pos        [1xP]   grid stimulus id (index into stim_id)
%   resp       [1xP]   activity in the response window (count or mean envelope)
%   base       [1xP]   activity in the baseline window (same units as resp)
%   is_target  [1xP]   logical, true for the [0,0] fixation-position stimulus
%   stim_id    [RxC]   grid position -> stimulus id map
%   x_dva      [1xC]   x grid coordinates (dva)
%   y_dva      [1xR]   y grid coordinates (dva)
%   params     struct  fields (all optional, sensible defaults below):
%     .response_mode  'diff' (resp-base) | 'percent' (100*(resp-base)/base)
%     .min_response_thr  responsiveness gate on the peak map value
%     .cursor_radius, .rf_sd_boundary, .overlap_grid_step, .alpha
%     .signal_label   char, stored in out.signal for provenance
%
% OUTPUT  out  struct with the SAME field names as the original
%              RF_characterisation (out.fixpos.*, out.rf.*, out.params.*) so
%              downstream code reads sorted and MUAe RFs identically.
%
% Felix Schneider, CNL  (drafted with Claude)
%
% Version history
%   1.0  (2026-07-03)  Initial version — generalised from RF_characterisation.

% =========================================================================
%% Parameters / defaults
% =========================================================================

def = struct('response_mode','diff', 'min_response_thr',5, ...
             'cursor_radius',1.75, 'rf_sd_boundary',2, ...
             'overlap_grid_step',0.05, 'alpha',0.05, 'signal_label','');
if nargin < 8 || isempty(params); params = struct(); end
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(params, fn{i}) || isempty(params.(fn{i}))
        params.(fn{i}) = def.(fn{i});
    end
end

cursor_radius   = params.cursor_radius;
cursor_area     = pi * cursor_radius^2;
rf_sd_boundary  = params.rf_sd_boundary;
rf_boundary_thr = rf_sd_boundary^2;
overlap_step    = params.overlap_grid_step;

pos       = pos(:).';
resp      = double(resp(:).');
base      = double(base(:).');
is_target = logical(is_target(:).');

% Per-presentation baseline-corrected response.
switch lower(params.response_mode)
    case 'percent'
        safe_base = base;  safe_base(safe_base == 0) = NaN;
        r_corr    = 100 * (resp - base) ./ safe_base;
    otherwise   % 'diff'
        r_corr    = resp - base;
end

n_stim = stim_id(end, end);

% =========================================================================
%% Initialise output with safe defaults
% =========================================================================

out.signal = params.signal_label;

out.fixpos.n_presentations  = 0;
out.fixpos.baseline_vals    = [];
out.fixpos.response_vals    = [];
out.fixpos.mean_baseline    = NaN;
out.fixpos.mean_response    = NaN;
out.fixpos.delta            = NaN;
out.fixpos.wilcoxon_W       = NaN;
out.fixpos.p_value          = NaN;
out.fixpos.significant      = false;

out.rf.responsive            = false;
out.rf.fit_ok                = false;
out.rf.spike_map             = [];    % kept name for downstream compatibility
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

out.params = params;

% =========================================================================
%% PART A — Fixation-position response
% =========================================================================

if any(is_target)
    bl = base(is_target);
    rp = resp(is_target);
    np = numel(rp);

    if any(bl ~= rp)
        [pval, hsig, wst] = signrank(bl', rp', 'alpha', params.alpha);
        wstat = wst.signedrank;
    else
        pval = 1; hsig = false; wstat = 0;
    end

    out.fixpos.n_presentations = np;
    out.fixpos.baseline_vals   = bl;
    out.fixpos.response_vals   = rp;
    out.fixpos.mean_baseline   = mean(bl, 'omitnan');
    out.fixpos.mean_response   = mean(rp, 'omitnan');
    out.fixpos.delta           = mean(rp, 'omitnan') - mean(bl, 'omitnan');
    out.fixpos.wilcoxon_W      = wstat;
    out.fixpos.p_value         = pval;
    out.fixpos.significant     = logical(hsig);
end

% =========================================================================
%% PART B — RF map and Gaussian fit
% =========================================================================

[X, Y]        = meshgrid(x_dva, y_dva);
dist_from_fix = sqrt(X.^2 + Y.^2);
cursor_mask   = dist_from_fix <= cursor_radius;
xy_flat       = [X(:), Y(:)];

% Mean baseline-corrected activity at each grid position.
act_map = nan(size(stim_id));
for iPos = 1:n_stim
    vals = r_corr(pos == iPos);
    if ~isempty(vals)
        act_map(stim_id == iPos) = mean(vals, 'omitnan');
    end
end
out.rf.spike_map          = act_map;
out.rf.mean_act_in_cursor = mean(act_map(cursor_mask), 'omitnan');

% Quality gate — skip flat maps.
if max(act_map(:), [], 'omitnan') < params.min_response_thr
    return
end
out.rf.responsive = true;

% 2-D Gaussian fit: p = [amplitude, x0, y0, sigma_x, sigma_y, offset]
gauss2d = @(p, xy) p(1) .* exp( -( (xy(:,1)-p(2)).^2 ./ (2.*p(4).^2) + ...
                                   (xy(:,2)-p(3)).^2 ./ (2.*p(5).^2) ) ) + p(6);
p0  = [max(act_map(:),[],'omitnan'), 0, 0, 6, 6, 0];
lb  = [0, min(x_dva), min(y_dva), 0.5, 0.5, -Inf];
ub  = [Inf, max(x_dva), max(y_dva), 20, 20, Inf];
opt = optimoptions('lsqcurvefit', 'Display', 'off');

z = act_map(:);  valid = ~isnan(z);
if sum(valid) >= 6
    try
        pf = lsqcurvefit(gauss2d, p0, xy_flat(valid,:), z(valid), lb, ub, opt);
        out.rf.RF_x_dva = pf(2); out.rf.RF_y_dva = pf(3);
        out.rf.RF_sigma_x_dva = pf(4); out.rf.RF_sigma_y_dva = pf(5);
        out.rf.fit_ok = true;
    catch ME
        warning('fn_characterise_RF: fit failed (%s).', ME.message);
    end
end

if out.rf.fit_ok
    x0 = out.rf.RF_x_dva; y0 = out.rf.RF_y_dva;
    sx = out.rf.RF_sigma_x_dva; sy = out.rf.RF_sigma_y_dva;

    out.rf.RF_area_dva2    = pi * (rf_sd_boundary*sx) * (rf_sd_boundary*sy);
    out.rf.dist_to_fix_dva = sqrt(x0^2 + y0^2);

    ellipse = ((X-x0)./sx).^2 + ((Y-y0)./sy).^2;
    inside  = ellipse <= rf_boundary_thr;
    out.rf.mean_spk_inside_RF  = mean(act_map(inside),  'omitnan');
    out.rf.mean_spk_outside_RF = mean(act_map(~inside), 'omitnan');

    x_lo = min(x0 - rf_sd_boundary*sx, -cursor_radius) - overlap_step;
    x_hi = max(x0 + rf_sd_boundary*sx,  cursor_radius) + overlap_step;
    y_lo = min(y0 - rf_sd_boundary*sy, -cursor_radius) - overlap_step;
    y_hi = max(y0 + rf_sd_boundary*sy,  cursor_radius) + overlap_step;
    [Xg, Yg] = meshgrid(x_lo:overlap_step:x_hi, y_lo:overlap_step:y_hi);

    in_rf     = ((Xg-x0)./sx).^2 + ((Yg-y0)./sy).^2 <= rf_boundary_thr;
    in_cursor = Xg.^2 + Yg.^2 <= cursor_radius^2;
    in_both   = in_rf & in_cursor;

    out.rf.overlap_flag          = any(in_both(:));
    out.rf.overlap_area_dva2     = sum(in_both(:)) * overlap_step^2;
    out.rf.overlap_pct_of_RF     = 100 * out.rf.overlap_area_dva2 / out.rf.RF_area_dva2;
    out.rf.overlap_pct_of_cursor = 100 * out.rf.overlap_area_dva2 / cursor_area;
end

end % fn_characterise_RF
