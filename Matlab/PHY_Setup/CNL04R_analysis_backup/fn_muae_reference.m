function ref = fn_muae_reference(t_us, v, k)
% FN_MUAE_REFERENCE  Pooled running reference level for the MUAe normalisation.
%
% Every baseline window (mean envelope over the 300 ms before a CPR cycle onset
% or a catch-trial cursor onset) gets as its reference the mean of the K nearest
% valid windows before it and the K nearest valid windows after it in time.
% Windows are pooled over trial types (solo, dyad, catch) and text-cue status,
% and the window itself is left out. The reference therefore follows slow drift
% of the envelope, but not the level of one trial or of one trial type
% (fixation timing, joystick hold and text cue differ between solo, dyad and
% catch trials). At the session edges fewer windows are used.
%
% Choice of K (2026-09-15; 19 sessions, 450 MT / 153 IPS channels). SD of the
% normalised cycle responses within solo and within dyad, relative to K = 10:
%   MT   K = 3 +4.7 %, 5 +2.8 %, 15 +3.5 %, 20 +6.1 %, 30 +17 %;
%        whole-session mean +54 %; own per-cycle baseline +44 %
%   IPS  K = 3-10 within 1 %, 15 +3 %; whole-session mean +57 %
% Truncating at the edges was better than extending, and the mean better than
% the median. Solo, dyad and catch trials are randomly interleaved (median run
% length 1), so each condition enters a reference at its overall share, for any K.
%
% INPUT
%   t_us  [n x 1]  window times (MWorks µs): cycle onsets / cursor onsets
%   v     [n x 1]  window means (mV). NaN or <= 0 marks an invalid window: it is
%                  not used for other references but still gets its own.
%   k     scalar   valid windows used on each side (default 10)
%
% OUTPUT
%   ref   [n x 1]  reference level per window (mV); NaN if no other valid window
%
% Felix Schneider, CNL  (drafted with Claude)
%
% Version history
%   1.0  (2026-09-15)  Initial version.

if nargin < 3 || isempty(k); k = 10; end
t_us = double(t_us(:));
v    = double(v(:));
ref  = nan(numel(t_us), 1);

[ts, o]  = sort(t_us);                       % time order; ties keep input order
vs       = v(o);
ok       = isfinite(ts) & isfinite(vs) & vs > 0;
val      = vs(ok);                           % valid windows in time order
n_val    = numel(val);
n_before = cumsum(ok) - ok;                  % valid windows before each window

for j = 1:numel(ts)
    if ~isfinite(ts(j)); continue; end
    nb  = n_before(j);                       % val(1:nb) lie before window j
    na  = n_val - nb - ok(j);                % val(nb+ok(j)+1:end) lie after it
    sel = [nb - min(k, nb) + 1 : nb, nb + ok(j) + 1 : nb + ok(j) + min(k, na)];
    if ~isempty(sel)
        ref(o(j)) = mean(val(sel));
    end
end

end % fn_muae_reference
