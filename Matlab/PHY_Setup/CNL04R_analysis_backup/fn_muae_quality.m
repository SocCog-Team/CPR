function phy = fn_muae_quality(phy, varargin)
% FN_MUAE_QUALITY  Inclusion / quality gate for MUAe channels.
%
% The MUAe analogue of PHY_quality_assessment (which only handles spikes/MUAt).
% A channel is INCLUDED if it has a valid RF fit OR a significant envelope
% onset response — with a minimum number of usable cycles. Manual drift-range
% curation is intentionally omitted: MUAe is drift-robust, so the whole session
% is used.
%
% Onset test: per CPR cycle, mean envelope in the response window (default
% 0-200 ms after cycle onset) vs the pre-cycle baseline (brain.CPR.muae.
% (ch).baseline). Wilcoxon signed-rank across cycles; "significant" = a
% significant CHANGE in either direction (increase OR suppression), so
% suppressed units are kept and not cherry-picked away. onset_incr records
% the sign (true = increase).
%
% INPUT
%   phy   struct with phy.brain.CPR.muae.(chXXX_mua) and (optionally)
%         phy.brain.RF.(chXXX_mua) from PHY_preprocessing_v3.
%
% OPTIONAL NAME-VALUE
%   'onset_win_ms' [0 200]   response window (ms rel. cycle onset)
%   'alpha'        0.05      significance level
%   'min_cycles'   30        minimum usable cycles to assess a channel
%
% OUTPUT
%   phy   with phy.brain.CPR.muae.include — a table (parallel to
%         brain.CPR.spks.include) with columns:
%           unit_ID  n_cycles  onset_p  onset_sig  onset_incr  mod_pct  rf_fit_ok  inclusion_flag
%
% Felix Schneider, CNL  (drafted with Claude)
%
% Version history
%   1.0  (2026-07-03)  Initial version.

p = inputParser;
p.addParameter('onset_win_ms', [0 200], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('alpha',        0.05,    @(x) isscalar(x) && x > 0);
p.addParameter('min_cycles',   30,      @(x) isscalar(x) && x >= 0);
p.parse(varargin{:});
o = p.Results;

if ~isfield(phy.brain, 'CPR') || ~isfield(phy.brain.CPR, 'muae')
    warning('fn_muae_quality:noMUAe', 'brain.CPR.muae missing — nothing to assess.');
    return
end

chans = setdiff(fieldnames(phy.brain.CPR.muae), {'include'});

vn = {'unit_ID','n_cycles','onset_p','onset_sig','onset_incr','mod_pct','rf_fit_ok','inclusion_flag'};
vt = {'string','double','double','logical','logical','double','logical','logical'};
T  = table('Size', [numel(chans) numel(vn)], 'VariableTypes', vt, 'VariableNames', vn);

for i = 1:numel(chans)
    ch = chans{i};
    M  = phy.brain.CPR.muae.(ch);
    nc = numel(M.env);

    % Per-cycle response = mean envelope in the onset window.
    resp = nan(1, nc);
    base = double(M.baseline);
    for iCyc = 1:nc
        if isempty(M.env{iCyc}); continue; end
        t_ms = double(M.t_us{iCyc}) / 1e3;                 % µs -> ms rel. cycle onset
        sel  = t_ms >= o.onset_win_ms(1) & t_ms <= o.onset_win_ms(2);
        if any(sel); resp(iCyc) = mean(double(M.env{iCyc}(sel))); end
    end

    valid = ~isnan(resp) & ~isnan(base);
    nv    = sum(valid);

    if nv >= o.min_cycles && any(resp(valid) ~= base(valid))
        pval       = signrank(resp(valid), base(valid), 'alpha', o.alpha);
        mb         = mean(base(valid));
        onset_sig  = pval < o.alpha;                          % STIMULUS-DRIVEN: significant change either way
        onset_incr = onset_sig && mean(resp(valid)) > mb;     % ...and it is an INCREASE (else a suppression)
        mod_pct    = 100 * (mean(resp(valid)) - mb) / mb;     % signed: + increase, - suppression
    else
        pval = NaN; onset_sig = false; onset_incr = false; mod_pct = NaN;
    end

    rf_ok = isfield(phy.brain, 'RF') && isfield(phy.brain.RF, ch) && ...
            isfield(phy.brain.RF.(ch), 'rf') && phy.brain.RF.(ch).rf.fit_ok;

    T.unit_ID(i)        = string(ch);
    T.n_cycles(i)       = nv;
    T.onset_p(i)        = pval;
    T.onset_sig(i)      = onset_sig;
    T.onset_incr(i)     = onset_incr;
    T.mod_pct(i)        = mod_pct;
    T.rf_fit_ok(i)      = rf_ok;
    T.inclusion_flag(i) = (rf_ok || onset_sig) && nv >= o.min_cycles;
end

phy.brain.CPR.muae.include = T;

fprintf('  MUAe quality: %d/%d channels included (RF-fit or onset response, incl. suppression) | %d incr / %d suppr.\n', ...
        sum(T.inclusion_flag), height(T), ...
        sum(T.inclusion_flag & T.onset_incr), sum(T.inclusion_flag & T.onset_sig & ~T.onset_incr));

end % fn_muae_quality
