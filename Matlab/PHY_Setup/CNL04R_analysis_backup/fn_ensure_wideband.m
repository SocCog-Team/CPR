function raw = fn_ensure_wideband(dest_dir, exp_info, chan_idx, pl2_fqn, pl2)
% FN_ENSURE_WIDEBAND  Return one channel's wideband, extracting it on demand.
%
% Lazy wrapper around the PLX_writeWideband logic:
%   * if the per-channel wb file already exists in DEST_DIR -> load and return it;
%   * else extract that channel from the PL2 (Plexon SDK), SAVE the wb file
%     (matching PLX_writeWideband's naming/fields), and return it. On the first
%     extraction of a session it also writes timestamps.mat (fragment times).
%
% INPUT
%   dest_dir   char     session folder (where wb / timestamps live)
%   exp_info   cellstr  split(fname,'_') -> {date, monkey, task, block, .., rec, ..}
%                       (wb name uses {1} date, {2} monkey, {6} rec, {4} block)
%   chan_idx   scalar   analog channel index (1..64), as in PLX_writeWideband
%   pl2_fqn    char     full path to the .pl2 file (only needed for extraction)
%   pl2        struct   PL2GetFileIndex(pl2_fqn) (only needed for extraction; [] if
%                       you already know the wb file exists)
%
% OUTPUT
%   raw   struct  fields fname, Fs, ad_raw_mv, gain_rec, conv_coeff_mv
%                 (== PLX_writeWideband); [] if the channel was not recorded.
%
% Felix Schneider, CNL  (drafted with Claude)
%
% Version history
%   1.0  (2026-07-03)  Initial version — PLX_writeWideband as a lazy function.

wb_name = sprintf('%s_%s_%s_%s_ch%02d_wb', ...
                  exp_info{1}, exp_info{2}, exp_info{6}, exp_info{4}, chan_idx);
wb_file = fullfile(dest_dir, [wb_name '.mat']);

% --- fast path: wb already on disk --------------------------------------
if isfile(wb_file)
    S   = load(wb_file, 'raw');
    raw = S.raw;
    return
end

% --- otherwise extract from the PL2 -------------------------------------
if isempty(pl2) || isempty(pl2_fqn) || ~isfile(pl2_fqn)
    error('fn_ensure_wideband:noSource', ...
          'wb file %s.mat missing and no PL2 available to extract it.', wb_name);
end

addpath('/Users/cnl/Desktop/CPR/PlexonMatlabOfflineFilesSDK/');

ad = PL2Ad(pl2_fqn, pl2.AnalogChannels{chan_idx}.Name);
if isempty(ad.Values)
    raw = [];                      % channel present but not recorded
    return
end

raw               = struct();
raw.fname         = wb_name;
raw.Fs            = ad.ADFreq;
raw.ad_raw_mv     = ad.Values;
raw.gain_rec      = pl2.AnalogChannels{chan_idx}.TotalGain;
raw.conv_coeff_mv = pl2.AnalogChannels{chan_idx}.CoeffToConvertToUnits;

% Reconstruct per-sample timestamps once per session (from PL2 fragments).
ts_file = fullfile(dest_dir, 'timestamps.mat');
if ~isfile(ts_file)
    timestamps = zeros(1, sum(ad.FragCounts));
    p = 0;
    for f = 1:numel(ad.FragTs)
        n                    = ad.FragCounts(f);
        timestamps(p+1:p+n)  = ad.FragTs(f) + (0:n-1) / ad.ADFreq;
        p                    = p + n;
    end
    save(ts_file, 'timestamps', '-v7.3');
end

save(wb_file, 'raw', '-v7.3');

end % fn_ensure_wideband
