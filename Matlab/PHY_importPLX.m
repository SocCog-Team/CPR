function out = PHY_importPLX(OpenedFileName, channel, plotFlag)

% OpenedFileName 	.String, File name
% channel           .Vector, Channel numbers to be imported
%
% Reference script: readall.m
% Felix Schneider, 03.04.2023

if nargin < 3
    plotFlag = true;
end

if nargin < 2 | isnan(channel)
    channel = [];
end

% Open a plx file
[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(OpenedFileName);

% Display recording information
disp(['Opened File Name: ' OpenedFileName]);
disp(['Version: ' num2str(Version)]);
disp(['Frequency : ' num2str(Freq)]);
disp(['Comment : ' Comment]);
disp(['Date/Time : ' DateTime]);
disp(['Duration : ' num2str(Duration)]);
disp(['Num Pts Per Wave : ' num2str(NPW)]);
disp(['Num Pts Pre-Threshold : ' num2str(PreThresh)]);

% Get a/d data into a cell array
% Note that analog ch numbering starts at 0, not 1 in the data, but the
% 'allad' cell array is indexed by ich+1
count                          = 0;

    % Channel Info
    [~,adnames]                 = plx_adchan_names(OpenedFileName);
    ch_id                      	= cellstr(adnames);
    
    for iCh = 0:nchannels-1
        chIdx = channel(iCh+1);
        
        %%% Wideband
        [adfreq, nad, tsad, fnad, allad] = plx_ad(OpenedFileName, chIdx);
        
        out.(['WB' num2str(chIdx)]).data = allad;
        out.(['WB' num2str(chIdx)]).fs = adfreq;
        out.(['WB' num2str(chIdx)]).nSample = nad;
        
        %%% Spikes
        
        %%% Field potentials
        
        count = count + 1;
    end
end

