%% Load file
close all

load('/Users/cnl/Documents/DATA/Nilan/spike_sorting/20250903_rec050/summary_20250903_nil_CPR_block1_phy4_rec050_ann.mat')

%%
tstep               = .001;
time                = [0:tstep:1];
figure;hold on
[all, sdf]          = FR_estimation(phy.brain.RF.raw.ch002_neg.spks_ts, time, true);

figure
plot(mean(sdf))

%% RF map
plot_RF_map(phy.brain.RF.stim_id, phy.brain.RF.stim_pos, phy.brain.RF.ch012_neg.nSpikes)

%% RF raster
plot_RF_raster(phy.brain.RF.raw.ch012_neg.spks_ts,100);

%% CPR tuning
close all

str_chan            = 'ch009_neg';
str_unit            = 'unit2';
incl                = state.include.([str_chan '_' str_unit]);
FR                  = state.spk_n.([str_chan '_' str_unit])(incl) ./ state.dur_s(incl);
stim_dir            = deg2rad(state.rdp_dir(incl));
stim_dir(FR==0)     = [];
FR(FR==0)           = [];

figure;
h = polaraxes; hold on
plot_CPRtuning(h,FR,stim_dir, true,[])

% Coherence-based tuning
FR                  = state.spk_n.([str_chan '_' str_unit]) ./ state.dur_s;
snr                 = unique(state.rdp_coh);
snr(snr==0)         = [];
stim_coh            = state.rdp_coh;
stim_dir            = deg2rad(state.rdp_dir);
col                 = cool(3);

f                   = figure;
h                   = polaraxes;
for iCoh = 1:length(snr)
    coh_idx         = stim_coh == snr(iCoh);
    plot_CPRtuning(h,FR(coh_idx & FR~=0),stim_dir(coh_idx & FR~=0),false,col(iCoh,:));
end
lg = legend(cellfun(@num2str,{snr(1) snr(2) snr(3)},'UniformOutput',false));

%% CPR cycle onset raster.

str_chan            = 'ch009_neg';
str_unit            = 'unit2';
cpr_spk_times       = [];
for iCyc = 1:length(phy.stim.cpr_cyle)
    cpr_spk_times{iCyc} = double(phy.brain.CPR.spks.(str_chan).(str_unit){iCyc});% ./1e6;
    cpr_spk_times{iCyc}(cpr_spk_times{iCyc} > 1e6) = [];
end

cycidx = phy.brain.CPR.spks.include.cyc_id{phy.brain.CPR.spks.include.unit_ID == [str_chan '_' str_unit]};

f                   = figure;
ax                  = subplot(3,3,1:6); hold on
tstep               = .001;
time                = [0:tstep:1];
[all, sdf]          = FR_estimation(cpr_spk_times(cycidx), time, true);
ax.XLim             = [0 1];
% ax.XLabel.String    = 'time [s]';
ax.YLabel.String    = 'CPR cycle [#]';
ax.FontSize         = 16;

ax                  = subplot(3,3,7:9);
pl                  = plot(mean(sdf),'k', 'LineWidth', 1.5);
ax.XLim             = [0 1000];
ax.XLabel.String    = 'time [ms]';
ax.YLabel.String    = 'FR [Hz]';
ax.Box              = 'off';
ax.FontSize         = 16;

% State-wise CPR raster
% f                   = figure;hold on
% [all, sdf]          = FR_estimation(state.spk_ts.([str_chan '_' str_unit]), time, true);

%% CPR spike - joystick correlation
addpath /Users/cnl/Desktop/CPR/CircStat2012a

str_chan            = 'ch012_neg';
str_unit            = 'unit2';
snr                 = unique(state.rdp_coh);
stim_coh            = state.rdp_coh;
snr(snr==0)         = [];
col                 = cool(3);
FR                  = state.spk_n.([str_chan '_' str_unit]) ./ state.dur_s;
tlt                 = cellfun(@mean,state.js_monk_tlt);

for iState = 1:length(state.dur_s)
    js_dev              = rad2deg(circ_dist(deg2rad(state.js_monk_dir{iState}),deg2rad(state.rdp_dir(iState)))); % Get circular distance to RDP direction
    state.js_acc{iState}= abs(1 - abs(js_dev / 180));                               % Calculate accuracy
end
acc                 = cellfun(@mean,state.js_acc);

figure; hold on
sc = scatter(tlt(FR~=0),FR(FR~=0),'k.');
% sc = scatter(acc(FR~=0),FR(FR~=0),'k.');
sc.SizeData = 50;
[r,p] = corrcoef(tlt(FR~=0),FR(FR~=0))
% [r,p] = corrcoef(acc(FR~=0),FR(FR~=0))
lsl = lsline;
lsl.Color = 'r';
lsl.LineWidth = 2;

% for iCoh = 1:length(snr)
%     coh_idx         = stim_coh == snr(iCoh);
%     sc = scatter(tlt(coh_idx & FR~=0),FR(coh_idx & FR~=0));
%     sc.Marker = '.';
%     sc.SizeData = 50;
%     sc.MarkerEdgeColor = col(iCoh,:);
%     [r,p] = corrcoef(tlt(coh_idx & FR~=0),FR(coh_idx & FR~=0))
%     lsl = lsline;
% end
%
% for iCoh = 1:length(snr)
%     lsl(iCoh).Color = col(iCoh,:);
%     lsl(iCoh).LineWidth = 2;
% end

xlabel('JS tilt [norm]')
ylabel('FR [Hz]')
title(['r: ' num2str(r(2)) ' | p: ' num2str(p(2))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gauss = fit_gaussian(spks, time)

sigma               = .005;                     % Width of gaussian/window [s]

% For every spike
for iSpk = 1:length(spks)

    % Center gaussian at spike time
    mu              = spks(iSpk);

    % Calculate gaussian
    p1              = -.5 * ((time - mu)/sigma) .^ 2;
    p2              = (sigma * sqrt(2*pi));
    gauss(iSpk,:)   = exp(p1) ./ p2;
end
end

function [all, sdf] = FR_estimation(spike_times, time, plot_flag)

sdf                 = [];
all                 = [];

for iTrial = 1:length(spike_times)

    spks            = double(spike_times{iTrial}') ./1e6;  	% Get all spikes of respective trial

    all             = [all spks];                   % Concatenate spikes of all trials
    xspikes         = repmat(spks,3,1);             % Replicate array
    yspikes      	= nan(size(xspikes));           % NaN array

    if ~isempty(yspikes)
        yspikes(1,:) = iTrial-1;                   	% Y-offset for raster plot
        yspikes(2,:) = iTrial;
    end

    % Plot trial raster
    if plot_flag
        pl           = plot(xspikes, yspikes, 'Color', 'k', 'LineWidth',1.25);
    end
    % Spike density function
    if isempty(spks)
    else

        % Fit gaussian to spikes
        gauss           = fit_gaussian(spks, time);

        % Sum over all distributions to get spike density function
        sdf(iTrial,:)	= sum(gauss,1);
    end

end
end

function plot_RF_map(stim_id, stim_pos, nSpikes)

for iUnit = 1:size(nSpikes,1)
    f                   	= figure;
    mat_sum                 = nan(size(stim_id));

    for iPos = 1:max(max(stim_id))
        mat_sum(stim_id == iPos) = sum(nSpikes(iUnit, stim_pos == iPos));
    end

    % hm                      = imagesc(mat_sum);
    hm                      = imagesc(imgaussfilt(mat_sum,.75));

    x_position              = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
    y_position              = fliplr([-12 -9 -6 -3 0 3 6 9 12]);

    ax                      = gca;
    ax.XTick                = [1 5 9];
    ax.XTickLabel           = {num2str(x_position(1)) num2str(x_position(5)) num2str(x_position(9))};
    ax.YTick                = [1 5 9];
    ax.YTickLabel           = {num2str(y_position(1)) num2str(y_position(5)) num2str(y_position(9))};
    ax.XLabel.String        = 'dva';
    ax.YLabel.String        = 'dva';
    ax.FontSize             = 16;
    colormap(jet(256)); colorbar
end
end

function sdf = plot_RF_raster(in, trl_end)
tstep       = .001;
time        =  [0:tstep:1];

f                   = figure; hold on
[all, sdf]          = FR_estimation(in(1:trl_end), time, true);
ax                  = gca;
ax.XTick            = [0:.2:2.4];
ax.YLim             = [1 trl_end];
numericLabels       = str2double(ax.XTickLabel);
newLabels           = numericLabels - .2;
ax.XTickLabel       = arrayfun(@num2str, newLabels, 'UniformOutput', false);
ax.XLabel.String    = 'time [s]';
ax.YLabel.String    = 'Trial [#]';
ax.FontSize         = 16;

for iLine = 1:(2200/200)
    ln              = line([.2*iLine .2*iLine],[0 trl_end],'LineStyle',':','Color','k');
end

% RF spike density
figure;hold on
plot(mean(sdf))
ax                  = gca;
ax.XLim             = [0 1000];
ax.YLim             = [0 ceil(max(mean(sdf)))];
ax.XTick            = 0:100:1000;
numericLabels       = str2double(ax.XTickLabel);
newLabels           = numericLabels - 200;
ax.XTickLabel       = arrayfun(@num2str, newLabels, 'UniformOutput', false);
ax.XLabel.String    = 'time [ms]';
ax.YLabel.String    = 'FR [Hz]';
ax.FontSize         = 16;
ln                  = line([200 200],[0 ceil(max(mean(sdf)))],'LineStyle',':','Color','k');

end

function plot_CPRtuning(ax,r, theta, data_flag, col)
hold on
if data_flag
    ps = polarscatter(ax,theta,r);
    ps.Marker = '.';
    ps.SizeData = 50;
    ps.MarkerEdgeColor = 'k';
end
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';

nbins = 360/30;
edges = linspace(0, 2*pi, nbins+1);
r_mean_bin = zeros(1, nbins);
theta_bin = zeros(1, nbins);

for i = 1:nbins
    idx = theta >= edges(i) & theta < edges(i+1);
    r_mean_bin(i) = mean(r(idx));
    theta_bin(i) = (edges(i)+edges(i+1))/2;
end

if isempty(col)
    polarplot(ax, theta_bin, r_mean_bin, 'r-', 'LineWidth', 2)
else
    polarplot(ax, theta_bin, r_mean_bin, 'LineStyle','-', 'LineWidth', 2, 'Color', col)
end

set(gca, 'FontSize', 16)
end
