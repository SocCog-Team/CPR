function [gain, raw_minmax, clipping_cnt] = check_for_clipping(pl2_name)

% addpath('/Users/fschneider/Documents/MATLAB/PlexonMatlabOfflineFilesSDK/');
addpath('/Users/cnl/Desktop/CPR/PlexonMatlabOfflineFilesSDK/')
pl2                         = PL2GetFileIndex(pl2_name);

for iChan = 1:length(pl2.SpikeChannels)
    disp(['Processing channel: ' num2str(iChan)])

    clear ad
    [ad]                    = PL2Ad(pl2_name, iChan);
    gain(iChan,:)           = pl2.AnalogChannels{iChan}.TotalGain;
    raw_minmax(iChan,:)     = [min(ad.Values) max(ad.Values)];
    clipping_cnt(iChan,:)   = [sum(ad.Values == min(ad.Values)) sum(ad.Values == max(ad.Values))];
    if sum(clipping_cnt(iChan,:)) > 10
        warning("you are SCREWED")
    end
end
end