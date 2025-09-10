function [gain, raw_minmax, clipping_cnt] = PLX_check_for_clipping(source_dir, pl2_name, dest_dir)
% This function checks if the raw, wideband Plexon data clipped.
% It requires 
%   .PL2 Plexon data file
%   .PlexonMatlabOfflineFilesSDK library
%
% If clipping is detected, remove affected channels in further processing 
% or skip recording
%
% Felix Schneider, CNL
%
% Version history
%   1.0     (fxs 2025-09-05) Initial version.

addpath('/Users/cnl/Desktop/CPR/PlexonMatlabOfflineFilesSDK/')
pl2                         = PL2GetFileIndex([source_dir pl2_name]);

for iChan = 1:length(pl2.SpikeChannels)
    disp(['Processing channel: ' num2str(iChan)])

    clear ad
    [ad]                    = PL2Ad([source_dir pl2_name], iChan);
    gain(iChan,:)           = pl2.AnalogChannels{iChan}.TotalGain;
    raw_minmax(iChan,:)     = [min(ad.Values) max(ad.Values)];
    clipping_cnt(iChan,:)   = [sum(ad.Values == min(ad.Values)) sum(ad.Values == max(ad.Values))];
    if sum(clipping_cnt(iChan,:)) > 10
        warning("you are SCREWED")
    end
end

save([dest_dir '/INFO_clipping_' pl2_name(1:end-4) '.mat'], 'gain','raw_minmax','clipping_cnt', '-v7.3')

end


