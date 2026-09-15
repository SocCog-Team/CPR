%% Pre-flight check for the pooled MUAe baseline — run once before the full run
%
% Tests PHY_preprocessing_v3 2.4, fn_sort_MUAe_by_state 1.2 and
% fn_muae_reference 1.0 on ONE catch-trial session (rec094), in memory, with
% import_flag = false (about 5 min). Nothing is written into the session folders.
%   1. The preprocessing runs, and loading the .h5 gives the same exp / stim /
%      joy / brain as the 2026-09-11 run from the .mwk2. Differences are listed;
%      numeric differences up to 1e-6 (relative) are only counted, because the
%      .h5 stores floats as single precision. Added fields show as 'new field'.
%   2. fn_sort_MUAe_by_state with 'baseline' 'cycle' reproduces the 2026-09-11
%      state responses.
%   3. The pooled responses equal the exact re-expression of the per-cycle
%      ones, f * pct + 100 * (f - 1) with f = baseline / baseline_ref.
%   4. The stored references equal a brute-force computation, and the arc count
%      of each catch trial agrees with TRIAL_type (Catch0 = one arc, Catch1 = two).
% Report:  <out_root>/preflight_pooled_baseline_rec094.txt (appended per run)
% Values for the cross-check on the laptop:  preflight_pooled_baseline_rec094.mat
%
% Felix Schneider, CNL  (drafted with Claude), 2026-09-15

clear; close all;
addpath('/Users/cnl/Desktop/CPR/code');
addpath(genpath('/Users/cnl/Documents/GitLab/matlab4mworks'));

rec        = '20260626_nil_CPR_block1_phy4_rec094_fxs';
cfg_pth    = '/Users/cnl/Desktop/CPR/code/felix_nhp_solo.cfg';
source_dir = '/Users/cnl/Documents/DATA/Nilan/';
out_root   = '/Users/cnl/Documents/DATA/Nilan/muae/';
rec_info   = split(rec, '_');
dest_dir   = fullfile(out_root, sprintf('%s_%s_%s', rec_info{1}, rec_info{6}, rec_info{4}), '/');
tol        = 1e-6;   % relative tolerance for numeric differences
verdict    = {'CHECK', 'PASS'};

diary(fullfile(out_root, 'preflight_pooled_baseline_rec094.txt'));
fprintf('\n===== Pre-flight check, pooled MUAe baseline: %s (%s)\n', rec, char(datetime('now')));
try
    % ---- 1. preprocessing from the .h5 vs the 2026-09-11 summary (.mwk2) --------
    t0  = tic;
    phy = PHY_preprocessing_v3(rec, source_dir, dest_dir, cfg_pth, false, {'muae'});
    phy = fn_muae_quality(phy);
    fprintf('\n[1] Preprocessing from the .h5 took %.1f min. Differences to the 2026-09-11 summary:\n', toc(t0) / 60);
    S = load(fullfile(dest_dir, ['summary_' rec '_muae.mat']), 'phy');
    [n1, r1] = compare(S.phy, phy, 'phy', tol, 0, 0);
    clear S
    fprintf('    -> %d differences; %d numeric fields differ at rounding level only\n', n1, r1);

    % ---- 2. 'cycle' normalisation reproduces the 2026-09-11 state ----------------
    S  = load(fullfile(dest_dir, ['state_responses_' rec '_muae.mat']), 'state');
    st = rmfield(S.state, intersect(fieldnames(S.state), ...
                 {'muae', 'muae_t', 'muae_mean', 'muae_fs', 'muae_baseline', 'muae_ref_k'}));
    st_cyc = fn_sort_MUAe_by_state(st, phy.brain, phy.stim, 'norm', 'percent', 'baseline', 'cycle');
    fprintf('\n[2] ''cycle'' normalisation vs the 2026-09-11 state:\n');
    [n2, r2] = compare(S.state.muae_mean, st_cyc.muae_mean, 'state.muae_mean', tol, 0, 0);
    [n2, r2] = compare(S.state.muae,      st_cyc.muae,      'state.muae',      tol, n2, r2);
    [n2, r2] = compare(S.state.muae_t,    st_cyc.muae_t,    'state.muae_t',    tol, n2, r2);
    clear S
    st_cyc = rmfield(st_cyc, 'muae_t');
    fprintf('    -> %d differences; %d numeric fields differ at rounding level only\n', n2, r2);

    % ---- 3. pooled = exact re-expression of the per-cycle normalisation ---------
    st_pool = fn_sort_MUAe_by_state(st, phy.brain, phy.stim, 'norm', 'percent');
    chans   = setdiff(fieldnames(phy.brain.CPR.muae), {'include'});
    d_mean  = 0;   d_seg = 0;   n_pool_only = 0;
    for i = 1:numel(chans)
        ch = chans{i};   M = phy.brain.CPR.muae.(ch);
        f  = double(M.baseline(st.cIdx)) ./ double(M.baseline_ref(st.cIdx));
        e  = f .* st_cyc.muae_mean.(ch) + 100 * (f - 1);
        g  = st_pool.muae_mean.(ch);
        v  = isfinite(e);
        d_mean      = max([d_mean, abs(e(v) - g(v))]);
        n_pool_only = n_pool_only + nnz(~v & isfinite(g));
        for s = find(v)
            x     = f(s) * double(st_cyc.muae.(ch){s}(:)) + 100 * (f(s) - 1);
            d_seg = max([d_seg; abs(x - double(st_pool.muae.(ch){s}(:)))]);
        end
    end
    fprintf('\n[3] Pooled vs re-expressed per-cycle responses: max |diff| muae_mean %.2g, muae segments %.2g (single precision)\n', d_mean, d_seg);
    fprintf('    %d states normalised with the pooled reference only (own baseline invalid); state.muae_baseline = ''%s'', ref_k = %d\n', ...
            n_pool_only, st_pool.muae_baseline, st_pool.muae_ref_k);

    % ---- 4. references (brute force) and catch trials -----------------------------
    t_all = [phy.stim.cpr_cyle(:, 1); phy.stim.catch_cursor];
    d_ref = 0;   n_nan = 0;
    for i = 1:numel(chans)
        M  = phy.brain.CPR.muae.(chans{i});
        v  = [M.baseline(:); M.catch_baseline(:)];
        r  = [M.baseline_ref(:); M.catch_baseline_ref(:)];
        ok = isfinite(v) & v > 0;
        for j = 1:numel(v)
            b  = find(ok & t_all < t_all(j));   [~, o] = sort(t_all(b));   b = b(o);
            a  = find(ok & t_all > t_all(j));   [~, o] = sort(t_all(a));   a = a(o);
            bf = mean(v([b(max(1, end - M.ref_k + 1):end); a(1:min(M.ref_k, end))]));
            if isnan(bf) ~= isnan(r(j))
                n_nan = n_nan + 1;
            elseif ~isnan(bf)
                d_ref = max(d_ref, abs(bf - r(j)) / bf);
            end
        end
    end
    ty  = phy.stim.catch_type;   nc = phy.stim.catch_ncur;
    bad = nnz((startsWith(ty, 'Catch1') & nc ~= 2) | (startsWith(ty, 'Catch0') & nc ~= 1));
    fprintf('\n[4] References vs brute force: max relative diff %.2g, NaN mismatches %d\n', d_ref, n_nan);
    fprintf('    Catch trials with cursor: %d (one arc %d, two %d); arc count vs TRIAL_type mismatches %d; no TRIAL_type %d\n', ...
            numel(nc), nnz(nc == 1), nnz(nc == 2), bad, nnz(cellfun(@isempty, ty)));
    fprintf('    Cursor onset after TRIAL_start: median %.0f ms\n', median((phy.stim.catch_cursor - phy.stim.catch_trl(:, 1)) / 1e3));

    % ---- values for the laptop cross-check ----------------------------------------
    R.chan               = cellfun(@(c) phy.brain.CPR.muae.(c).chan_num, chans);
    R.cyc_on             = phy.stim.cpr_cyle(:, 1);
    R.catch_trl          = phy.stim.catch_trl;
    R.catch_cursor       = phy.stim.catch_cursor;
    R.catch_ncur         = phy.stim.catch_ncur;
    R.ref_k              = phy.brain.CPR.muae.(chans{1}).ref_k;
    R.baseline           = stack_rows(phy.brain.CPR.muae, chans, 'baseline');
    R.baseline_ref       = stack_rows(phy.brain.CPR.muae, chans, 'baseline_ref');
    R.catch_baseline     = stack_rows(phy.brain.CPR.muae, chans, 'catch_baseline');
    R.catch_baseline_ref = stack_rows(phy.brain.CPR.muae, chans, 'catch_baseline_ref');
    R.cIdx               = st.cIdx;
    R.muae_mean_cycle    = stack_rows(st_cyc.muae_mean,  chans, '');
    R.muae_mean_pooled   = stack_rows(st_pool.muae_mean, chans, '');
    save(fullfile(out_root, 'preflight_pooled_baseline_rec094.mat'), '-struct', 'R', '-v7.3');

    fprintf('\nVerdict: [1] %s  [2] %s  [3] %s  [4] %s\n', verdict{1 + (n1 == 0)}, verdict{1 + (n2 == 0)}, ...
            verdict{1 + (d_mean < 1e-9 && d_seg < 1e-3)}, verdict{1 + (d_ref < 1e-12 && n_nan == 0 && bad == 0)});
catch ME
    fprintf('\nPRE-FLIGHT ERROR\n%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
end
diary off


function [nd, nr] = compare(a, b, path, tol, nd, nr)
% List where A (2026-09-11) and B (this run) differ; the first 40 differences
% are printed. ND counts differences, NR numeric fields that differ by at most
% TOL (relative) only.
if isstruct(a) && isstruct(b)
    fa = fieldnames(a);   fb = fieldnames(b);
    % Loop over indices: for a one-field struct setdiff returns a 1x0 cell,
    % and a for loop over its transpose (0x1) runs once with an empty value.
    added = setdiff(fb, fa);   lost = setdiff(fa, fb);   both = intersect(fa, fb);
    for k = 1:numel(added); fprintf('    new field   %s.%s\n', path, added{k}); end
    for k = 1:numel(lost)
        if nd < 40; fprintf('  ! missing     %s.%s\n', path, lost{k}); end
        nd = nd + 1;
    end
    if numel(a) ~= numel(b)
        if nd < 40; fprintf('  ! size        %s (%d vs %d elements)\n', path, numel(a), numel(b)); end
        nd = nd + 1;   return
    end
    for i = 1:numel(a)
        for k = 1:numel(both)
            p = sprintf('%s.%s', path, both{k});
            if numel(a) > 1; p = sprintf('%s(%d).%s', path, i, both{k}); end
            [nd, nr] = compare(a(i).(both{k}), b(i).(both{k}), p, tol, nd, nr);
        end
    end
elseif iscell(a) && iscell(b)
    if ~isequal(size(a), size(b))
        if nd < 40; fprintf('  ! size        %s (%s vs %s)\n', path, mat2str(size(a)), mat2str(size(b))); end
        nd = nd + 1;   return
    end
    for i = 1:numel(a)
        [nd, nr] = compare(a{i}, b{i}, sprintf('%s{%d}', path, i), tol, nd, nr);
    end
elseif (isnumeric(a) || islogical(a)) && (isnumeric(b) || islogical(b))
    if isequaln(a, b); return; end
    if ~isequal(size(a), size(b))
        if nd < 40; fprintf('  ! size        %s (%s vs %s)\n', path, mat2str(size(a)), mat2str(size(b))); end
        nd = nd + 1;   return
    end
    x = double(a(:));   y = double(b(:));
    if ~isequal(isnan(x), isnan(y))
        if nd < 40; fprintf('  ! NaN         %s (%d vs %d NaN)\n', path, nnz(isnan(x)), nnz(isnan(y))); end
        nd = nd + 1;   return
    end
    d = max(abs(x - y));   s = max(abs(x));
    if d <= tol * max(s, 1)
        nr = nr + 1;
    else
        if nd < 40; fprintf('  ! value       %s (max |diff| %.3g, max |value| %.3g)\n', path, d, s); end
        nd = nd + 1;
    end
elseif ~isequaln(a, b)
    if nd < 40; fprintf('  ! differs     %s (%s vs %s)\n', path, class(a), class(b)); end
    nd = nd + 1;
end
end % compare


function X = stack_rows(S, chans, fld)
% Per-channel vectors S.(ch) (FLD empty) or S.(ch).(FLD) as rows of [nCh x n].
X = cell(numel(chans), 1);
for i = 1:numel(chans)
    if isempty(fld); x = S.(chans{i}); else; x = S.(chans{i}).(fld); end
    X{i} = double(x(:)).';
end
X = vertcat(X{:});
end % stack_rows
