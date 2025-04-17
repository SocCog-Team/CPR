addpath('/Users/fschneider/Documents/MATLAB/PlexonMatlabOfflineFilesSDK/');

% Get channel information
pl2file                     = '/Volumes/T7_Shield/plexon/Data_Nilan/fxs-CPR20250312-nil-010-01+01-sorted.pl2';
[pl2]                       = PL2GetFileIndex(pl2file);
wideband_data               = struct();
dest_dir                    = '/Users/fschneider/Desktop/plx_wb_20250312_Nilan/';
filter_cutoff               = [300 5000];

for iChan = 1:33
    % Only include channels that are enabled
    if pl2.EventChannels{iChan}.Enabled
        disp(['Process Channel: ' pl2.AnalogChannels{iChan}.Name])
   
        % Import analog data
        [ad]                = PL2Ad(pl2file, pl2.AnalogChannels{iChan}.Name);
        
        if iChan == 1
            % Reconstruct time vector: Each fragment starts at FragTs(i) and contains FragCounts(i) samples at ADFreq
            timestamps      = [];
            
            for f = 1:length(ad.FragTs)
                t0          = ad.FragTs(f);
                n_samples   = ad.FragCounts(f);
                ts          = t0 + (0:n_samples-1) / ad.ADFreq;
                timestamps  = [timestamps, ts];
            end
            save([dest_dir 'timestamps.mat'], 'timestamps', '-v7.3');
        end
        
        % Save signal in structure for each channel
        out                 = [];
        out.name            = pl2.AnalogChannels{iChan}.Name;
        out.Fs              = ad.ADFreq;
        out.ad_raw          = ad.Values;
        
        % Filter raw data
        % filtered_data     = bessel_filter_chunked(ad.Values, ad.ADFreq, filter_cutoff, 'bandpass');
        % out.values_filt   = filtered_data;
        
        save([dest_dir num2str(pl2.AnalogChannels{iChan}.Name) '.mat'], 'out', '-v7.3');
    end
end

fprintf('Done!');

%%
% idx = 1e3:1e6;
% figure
% plot(wb(3).timestamps(idx),wb(3).values_filt(idx))

% [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(pl2file)
% [n, ts] = plx_ts(filename, channel, unit) % PLX spike times!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filtered = bessel_filter_chunked(data, Fs, cutoff, filter_type, chunk_size)
% Filters large data in chunks with optional overlap padding to avoid edge artifacts
%
% INPUTS:
%   data        - raw signal
%   Fs          - sampling frequency in Hz
%   cutoff      - cutoff frequency in Hz
%   filter_type - 'low', 'high', 'bandpass'
%   chunk_size  - number of samples per chunk (e.g., 1e6)
%
% OUTPUT:
%   filtered    - filtered signal

if nargin < 5
    chunk_size = 1e6; % Default: 1 million samples per chunk
end

% Manual digital Bessel bandpass design
[b, a] = besself(6, 1);  % 6-pole analog lowpass prototype

if strcmpi(filter_type, 'bandpass')
    Wn = 2 * pi * cutoff;           % Hz → rad/s
    wo = sqrt(Wn(1) * Wn(2));       % center frequency
    bw = Wn(2) - Wn(1);             % bandwidth
    [b_bp, a_bp] = lp2bp(b, a, wo, bw);  % Bandpass transform
elseif strcmpi(filter_type, 'low')
    Wn = 2 * pi * cutoff;
    [b_bp, a_bp] = lp2lp(b, a, Wn);
elseif strcmpi(filter_type, 'high')
    Wn = 2 * pi * cutoff;
    [b_bp, a_bp] = lp2hp(b, a, Wn);
else
    error('filter_type must be ''low'', ''high'', or ''bandpass''.');
end

% Now apply bilinear transform with correct syntax
[b_digital, a_digital] = bilinear(b_bp, a_bp, Fs);

% Preallocate output
filtered = zeros(size(data));
N = length(data);

% Overlap length (padding)
overlap = 1000;  % ~25 ms at 40 kHz

% Process in chunks
idx = 1;
while idx <= N
    start_idx = max(1, idx - overlap);
    end_idx = min(N, idx + chunk_size + overlap - 1);
    chunk = data(start_idx:end_idx);
    
    % Filter with zero-phase
    chunk_filtered = filtfilt(b_digital, a_digital, chunk);
    
    % Assign valid central portion to output
    valid_start = idx;
    valid_end = min(N, idx + chunk_size - 1);
    out_range = (valid_start:valid_end) - start_idx + 1;
    filtered(valid_start:valid_end) = chunk_filtered(out_range);
    
    idx = idx + chunk_size;
end
end
