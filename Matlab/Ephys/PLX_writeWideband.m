addpath('/Users/fschneider/Documents/MATLAB/PlexonMatlabOfflineFilesSDK/');
cd /Volumes/T7_Shield/plexon/Data_Nilan

% Compile list of all data files
pl2_files = dir(fullfile('/Volumes/T7_Shield/plexon/Data_Nilan', '*.pl2'));

%%% File-Loop %%%
for iFile = 1:length(pl2_files)
    [pl2]                       = PL2GetFileIndex([pl2_files(iFile).folder '/' pl2_files(iFile).name]);
    wideband_data               = struct();
    exp_info                    = split(pl2_files(iFile).name,'_');
    dest_dir                    = [exp_info{1} '_wb_mat'];
    
    if ~isfolder(dest_dir)
        mkdir(dest_dir);
    end
    
    %%% Channel-Loop %%%
    for iChan = 1:64
        
        % Only include enabled channels
        if pl2.AnalogChannels{iChan}.Enabled
            disp(['Process Channel: ' pl2.AnalogChannels{iChan}.Name])
            
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
                save([dest_dir '/timestamps.mat'], 'timestamps', '-v7.3');
            end
            
            % Save raw signal in structure for each channel
            raw                 = [];
            raw.name            = [exp_info{1} '_' exp_info{2} '_' exp_info{6} '_' exp_info{4} '_ch' num2str(pad(iChan, 2, 'left', '0')) '_wb'];
            raw.Fs              = ad.ADFreq;
            raw.ad_raw          = ad.Values;
            
            save([dest_dir '/' raw.name '.mat'], 'raw', '-v7.3');
        end
    end
end

fprintf('Done!');

