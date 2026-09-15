function state = fn_sort_MUAe_by_state(state, brain, stim, varargin)
% FN_SORT_MUAE_BY_STATE  Add per-state MUAe segments to an existing state struct.
%
% Companion to the state segmentation (sort_states in PHY_main_analysis_v4, or
% PHY_sort_spikes_by_state): it reuses the state boundaries already computed
% (state.boundaries, state.cIdx, state.dur_s) and slices the per-cycle MUAe
% envelope (brain.CPR.muae, from PHY_preprocessing_v3) into the same states,
% normalised to a baseline reference level.
%
% Baseline reference ('baseline'):
%   'pooled' (default)  pooled running reference of the cycle (fn_muae_reference):
%                       mean of the ref_k nearest valid baseline windows before
%                       and after the cycle, CPR cycles and catch trials pooled,
%                       the cycle's own window left out. Taken from
%                       brain.CPR.muae.(ch).baseline_ref (PHY_preprocessing_v3
%                       >= 2.4) when built with the same ref_k; otherwise
%                       computed here from the stored windows (older summaries
%                       have no catch windows: CPR cycles only, with a warning).
%   'cycle'             the cycle's own pre-cycle baseline (versions <= 1.1). It
%                       differs systematically between solo, dyad and catch
%                       trials (fixation timing, joystick hold, text cue), and
%                       dividing by it carries that difference into every state
%                       response (MT dyad - solo: +0.56 pp vs +0.12 pp pooled).
%
% INPUT
%   state   struct   state segmentation; needs .boundaries [nState x 2]
%                    (onset/offset, MWorks µs), .cIdx, .dur_s
%   brain   struct   with brain.CPR.muae.(chXXX_mua) from PHY_preprocessing_v3
%                    (per-cycle .env / .t_us with t_us increasing, .baseline;
%                    >= 2.4 also .baseline_ref, .ref_k, .catch_baseline)
%   stim    struct   with stim.cpr_cyle [nCyc x 2] (cycle onset/end, MWorks µs);
%                    >= 2.4 also stim.catch_cursor (catch-trial cursor onsets)
%
% OPTIONAL NAME-VALUE PAIRS
%   'norm'      char     'percent' (default) = 100*(env-ref)/ref
%                        'subtract' = env - ref
%                        'divide'   = env / ref
%                        'none'     = raw envelope
%   'baseline'  char     'pooled' (default) | 'cycle', see above
%   'ref_k'     scalar   windows on each side for 'pooled' (default 10)
%
% OUTPUT  state, augmented with, for every channel key CH = 'chXXX_mua':
%   state.muae.(CH){s}       normalised MUAe segment for state s
%   state.muae_t.(CH){s}     sample times (µs) RELATIVE TO state onset
%   state.muae_mean.(CH)(s)  mean normalised MUAe over the whole state
%   state.muae_fs.(CH)       envelope sample rate (Hz)
% and the provenance
%   state.muae_baseline      'pooled' | 'cycle'
%   state.muae_ref_k         ref_k ('pooled' only)
%
% Downstream (poster) reads these via compute_muae(state, CH, win_ms, ref),
% the MUAe analogue of compute_fr — the uid is the channel key 'chXXX_mua'.
%
% Felix Schneider, CNL  (drafted with Claude)
%
% Version history
%   1.0  (2026-07-03)  Initial version.
%   1.1  (2026-09-11)  State window found by binary search on the (increasing)
%                      cycle time vector; only the state's samples are cast to
%                      double. Output unchanged.
%   1.2  (2026-09-15)  Normalised to the pooled running reference by default
%                      ('baseline' 'pooled', 'ref_k' 10); 'cycle' reproduces 1.1.
%                      Provenance in state.muae_baseline / state.muae_ref_k.

p = inputParser;
p.addParameter('norm',     'percent', @(x) any(strcmpi(x, {'percent','subtract','divide','none'})));
p.addParameter('baseline', 'pooled',  @(x) any(strcmpi(x, {'pooled','cycle'})));
p.addParameter('ref_k',    10,        @(x) isscalar(x) && x >= 1 && x == round(x));
p.parse(varargin{:});
norm_mode = lower(p.Results.norm);
base_mode = lower(p.Results.baseline);
ref_k     = p.Results.ref_k;

if ~isfield(brain, 'CPR') || ~isfield(brain.CPR, 'muae')
    warning('fn_sort_MUAe_by_state:noMUAe', ...
        'brain.CPR.muae missing — run PHY_preprocessing_v3 with signal_source ''muae'' first.');
    return
end

cyc      = double(stim.cpr_cyle);
n_cyc    = size(cyc, 1);
chans    = setdiff(fieldnames(brain.CPR.muae), {'include'});
n_state  = size(state.boundaries, 1);
warned   = false;

for iCh = 1:numel(chans)
    ch = chans{iCh};
    M  = brain.CPR.muae.(ch);

    % Reference level of every cycle (mV).
    if strcmp(base_mode, 'cycle')
        ref = double(M.baseline);
    elseif isfield(M, 'baseline_ref') && isfield(M, 'ref_k') && M.ref_k == ref_k
        ref = double(M.baseline_ref);
    else
        t_pool = cyc(:, 1);
        v_pool = double(M.baseline(:));
        if isfield(M, 'catch_baseline') && isfield(stim, 'catch_cursor')
            t_pool = [t_pool; double(stim.catch_cursor(:))];
            v_pool = [v_pool; double(M.catch_baseline(:))];
        elseif ~warned
            warning('fn_sort_MUAe_by_state:cprOnlyPool', ...
                'Summary has no catch-trial baseline windows (PHY_preprocessing_v3 < 2.4) — pooled reference from CPR cycles only.');
            warned = true;
        end
        ref = fn_muae_reference(t_pool, v_pool, ref_k);
        ref = ref(1:n_cyc);
    end

    state.muae_fs.(ch)   = M.env_fs;
    state.muae.(ch)      = cell(1, n_state);
    state.muae_t.(ch)    = cell(1, n_state);
    state.muae_mean.(ch) = nan(1, n_state);

    for s = 1:n_state
        iCyc = state.cIdx(s);

        % Skip states whose cycle has no stored MUAe (defensive).
        if iCyc > numel(M.env) || isempty(M.env{iCyc})
            continue
        end

        c_on     = cyc(iCyc, 1);
        t_cyc    = M.t_us{iCyc};            % µs rel to cycle onset (increasing)
        baseline = ref(iCyc);

        % State window expressed relative to cycle onset (same frame as t_cyc):
        % samples with st_onset_rel <= t <= st_end_rel form one index range.
        st_onset_rel = state.boundaries(s, 1) - c_on;
        st_end_rel   = state.boundaries(s, 2) - c_on;
        i0  = count_lt(t_cyc, st_onset_rel) + 1;
        i1  = count_le(t_cyc, st_end_rel);
        seg = double(M.env{iCyc}(i0:i1));   % amplitude

        % Baseline normalisation.
        switch norm_mode
            case 'percent'
                seg = 100 * (seg - baseline) / baseline;
            case 'subtract'
                seg = seg - baseline;
            case 'divide'
                seg = seg / baseline;
            case 'none'
                % raw
        end

        state.muae.(ch){s}      = single(seg);
        state.muae_t.(ch){s}    = t_cyc(i0:i1) - st_onset_rel;   % µs rel to state onset
        state.muae_mean.(ch)(s) = mean(seg, 'omitnan');
    end
end

state.muae_baseline = base_mode;
if strcmp(base_mode, 'pooled')
    state.muae_ref_k = ref_k;
elseif isfield(state, 'muae_ref_k')
    state = rmfield(state, 'muae_ref_k');
end

end % fn_sort_MUAe_by_state


function c = count_le(t, x)
% Number of elements of increasing T that are <= X (binary search).
lo = 1;   hi = numel(t);   c = 0;
while lo <= hi
    mid = floor((lo + hi) / 2);
    if t(mid) <= x; c = mid; lo = mid + 1; else; hi = mid - 1; end
end
end % count_le


function c = count_lt(t, x)
% Number of elements of increasing T that are < X (binary search).
lo = 1;   hi = numel(t);   c = 0;
while lo <= hi
    mid = floor((lo + hi) / 2);
    if t(mid) < x; c = mid; lo = mid + 1; else; hi = mid - 1; end
end
end % count_lt
