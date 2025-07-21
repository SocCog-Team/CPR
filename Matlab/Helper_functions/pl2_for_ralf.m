addpath('/Users/fschneider/Documents/MATLAB/PlexonMatlabOfflineFilesSDK/');
pl2_name = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_electrophysiology/Nilan/pl2/20250312_nil_CPR_block1_phy4_rec010_fxs.pl2';

%% Import raw wideband data
[ad] = PL2Ad(pl2_name, 'WB07');

%% Get events
pl2 = PL2GetFileIndex(pl2_name);
Fs = pl2.TimestampFrequency;

%% Print all event channels
PL2Print(pl2.EventChannels);

%% This section is based on your code %%%
tmp_ts      = [];
for evChan = 1:length(pl2.EventChannels)
    if pl2.EventChannels{evChan}.NumEvents > 0 && contains(pl2.EventChannels{evChan}.Name, 'EVT')
        events  = PL2EventTs(pl2_name, evChan);
        tmp_ts  = [tmp_ts; events.Ts];
    end
end

[plx.tsTemp, sortIdx] = sort(tmp_ts);
    
plx.tsOrg(1) = plx.tsTemp(1);
for cx = 2:numel(plx.tsTemp)
    if ~((plx.tsOrg(length(plx.tsOrg)) > plx.tsTemp(cx) - 0.0005) && (plx.tsOrg(length(plx.tsOrg)) < plx.tsTemp(cx) + 0.0005))
        plx.tsOrg(length(plx.tsOrg)+1) = plx.tsTemp(cx);
    else
    end
end

%% Import sorted spikes
load([pth 'dataspikes_ch007_negthr.mat']) % ID + timestamp
spk_ts = cluster_class(:,2);

%% Plot
close all; 
figure; hold on

plot(ad.Values); % raw signal
scatter(spk_ts(1).*(Fs/1e3),0,'filled')
scatter(spk_ts(end).*(Fs/1e3),0,'filled')
scatter(plx.tsOrg(1)*Fs,0.1,'filled') % 40kHz sampling rate
scatter(plx.tsOrg(end)*Fs,0.1,'filled')
title('data.wideband.raw')
legend('data', '1st spike', 'Last spike','1st strobe','Last strobe')
xlabel('time [samples]')
ylabel('[mv]')

