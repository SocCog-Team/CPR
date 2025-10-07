% This script writes wideband electrophysiology data from raw Plexon files.
% It requires 
%   .PL2 Plexon data file
%   .PlexonMatlabOfflineFilesSDK library
%
% Felix Schneider, CNL
%
% Version history
%   1.0     (fxs 2025-09-05) Initial version.

addpath('/Users/cnl/Desktop/CPR/PlexonMatlabOfflineFilesSDK/');
cd /Users/cnl/Documents/DATA/Nilan/pl2/

% Compile list of all data files
pl2_files = dir(fullfile('/Users/cnl/Documents/DATA/Nilan/pl2/', '*.pl2'));
pl2_files = pl2_files(end) % tmp fix - add option to process specific file here

cd /Users/cnl/Documents/DATA/Nilan/spike_sorting/

%%% File-Loop %%%
for iFile = 1:length(pl2_files)
    disp(['Processing: ' pl2_files(iFile).name])
    [pl2]                       = PL2GetFileIndex([pl2_files(iFile).folder '/' pl2_files(iFile).name]);
    wideband_data               = struct();
    exp_info                    = split(pl2_files(iFile).name,'_');
    dest_dir                    = [exp_info{1} '_' exp_info{6}];
    
    if ~isfolder(dest_dir)
        mkdir(dest_dir);
    else 
        continue
    end
    
    %%% Channel-Loop %%%
    for iChan = 1:64

        % Only include enabled channels
        if pl2.AnalogChannels{iChan}.Enabled
            disp(['Processing Channel: ' pl2.AnalogChannels{iChan}.Name])
            
            % Import analog data
            [ad]                = PL2Ad([pl2_files(iFile).folder '/' pl2_files(iFile).name], pl2.AnalogChannels{iChan}.Name);
            
            % Skip empty channels
            if isempty(ad.Values)
                continue
            end
            
            % Reconstruct timestamps
            if iChan == 1
                % Each fragment starts at FragTs(i) and contains FragCounts(i) samples at ADFreq
                timestamps      = [];
                
                for f = 1:length(ad.FragTs)
                    t0          = ad.FragTs(f);
                    n_samples   = ad.FragCounts(f);
                    ts          = t0 + (0:n_samples-1) / ad.ADFreq;
                    timestamps  = [timestamps, ts];
                end
                save(['/Users/cnl/Documents/DATA/Nilan/spike_sorting/' dest_dir '/timestamps.mat'], 'timestamps', '-v7.3');
            end
            
            % Save raw signal in structure for each channel
            raw                 = [];
            raw.fname           = [exp_info{1} '_' exp_info{2} '_' exp_info{6} '_' exp_info{4} '_ch' pad(num2str(iChan), 2, 'left', '0') '_wb'];
            raw.Fs              = ad.ADFreq;
            raw.ad_raw_mv       = ad.Values;
            raw.gain_rec        = pl2.AnalogChannels{iChan}.TotalGain;
            raw.conv_coeff_mv   = pl2.AnalogChannels{iChan}.CoeffToConvertToUnits;
            
            save(['/Users/cnl/Documents/DATA/Nilan/spike_sorting/' dest_dir '/' raw.fname '.mat'], 'raw', '-v7.3');
        end
    end
end

fprintf('Done!');

