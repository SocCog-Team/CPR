function state = fn_sort_MUAe_by_state(state, brain, stim, varargin)
% FN_SORT_MUAE_BY_STATE  Add per-state MUAe segments to an existing state struct.
%
% Companion to PHY_sort_spikes_by_state: it reuses the state boundaries that
% function already computed (state.boundaries, state.cIdx, state.dur_s) and
% slices the per-cycle MUAe envelope (brain.CPR.muae, from PHY_add_MUA) into the
% same states, applying the per-cycle baseline normalisation.
%
% Call it AFTER PHY_sort_spikes_by_state so `state` already carries the
% direction-change segmentation.
%
% INPUT
%   state   struct   output of PHY_sort_spikes_by_state; needs .boundaries
%                    [nState x 2] (onset/offset, MWorks µs), .cIdx, .dur_s
%   brain   struct   with brain.CPR.muae.(chXXX_mua) from PHY_add_MUA
%   stim    struct   with stim.cpr_cyle [nCyc x 2] (cycle onset/end, MWorks µs)
%
% OPTIONAL NAME-VALUE PAIRS
%   'norm'   char   'percent' (default) = 100*(env-baseline)/baseline
%                   'subtract' = env - baseline
%                   'divide'   = env / baseline
%                   'none'     = raw envelope
%
% OUTPUT  state, augmented with, for every channel key CH = 'chXXX_mua':
%   state.muae.(CH){s}       normalised MUAe segment for state s
%   state.muae_t.(CH){s}     sample times (µs) RELATIVE TO state onset
%   state.muae_mean.(CH)(s)  mean normalised MUAe over the whole state
%   state.muae_fs.(CH)       envelope sample rate (Hz)
%
% Downstream (poster) reads these via compute_muae(state, CH, win_ms, ref),
% the MUAe analogue of compute_fr — the uid is the channel key 'chXXX_mua'.
%
% Felix Schneider, CNL  (drafted with Claude)
%
% Version history
%   1.0  (2026-07-03)  Initial version.

p = inputParser;
p.addParameter('norm', 'percent', @(x) any(strcmpi(x, {'percent','subtract','divide','none'})));
p.parse(varargin{:});
norm_mode = lower(p.Results.norm);

if ~isfield(brain, 'CPR') || ~isfield(brain.CPR, 'muae')
    warning('fn_sort_MUAe_by_state:noMUAe', ...
        'brain.CPR.muae missing — run PHY_preprocessing_v3 with signal_source ''muae'' first.');
    return
end

cyc      = double(stim.cpr_cyle);
chans    = fieldnames(brain.CPR.muae);
n_state  = size(state.boundaries, 1);

for iCh = 1:numel(chans)
    ch = chans{iCh};

    state.muae_fs.(ch)   = brain.CPR.muae.(ch).env_fs;
    state.muae.(ch)      = cell(1, n_state);
    state.muae_t.(ch)    = cell(1, n_state);
    state.muae_mean.(ch) = nan(1, n_state);

    for s = 1:n_state
        iCyc = state.cIdx(s);

        % Skip states whose cycle has no stored MUAe (defensive).
        if iCyc > numel(brain.CPR.muae.(ch).env) || isempty(brain.CPR.muae.(ch).env{iCyc})
            continue
        end

        c_on     = cyc(iCyc, 1);
        env_cyc  = double(brain.CPR.muae.(ch).env{iCyc});     % amplitude, per cycle
        t_cyc    = brain.CPR.muae.(ch).t_us{iCyc};            % µs rel to cycle onset
        baseline = brain.CPR.muae.(ch).baseline(iCyc);

        % State window expressed relative to cycle onset (same frame as t_cyc).
        st_onset_rel = state.boundaries(s, 1) - c_on;
        st_end_rel   = state.boundaries(s, 2) - c_on;

        sel = t_cyc >= st_onset_rel & t_cyc <= st_end_rel;
        seg = env_cyc(sel);

        % Per-cycle baseline normalisation.
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
        state.muae_t.(ch){s}    = t_cyc(sel) - st_onset_rel;   % µs rel to state onset
        state.muae_mean.(ch)(s) = mean(seg, 'omitnan');
    end
end

end % fn_sort_MUAe_by_state
