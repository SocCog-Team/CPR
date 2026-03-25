%% Load file
clear all
close all

% source_dir        = '/Users/cnl/Documents/DATA/Nilan/spike_sorting/20250924_rec059_block1/';
% load([source_dir '/state_responses_20250924_nil_CPR_block1_phy4_rec059_fxs.mat'])
% load([source_dir '/summary_20250924_nil_CPR_block1_phy4_rec059_fxs.mat'])

% source_dir        = '/Users/cnl/Documents/DATA/Nilan/spike_sorting/20250926_rec061_block1/';
% load([source_dir '/state_responses_20250926_nil_CPR_block1_phy4_rec061_fxs.mat'])
% load([source_dir '/summary_20250926_nil_CPR_block1_phy4_rec061_fxs.mat'])
% 
dest_dir        = '/Users/cnl/Documents/DATA/Nilan/spike_sorting/20250925_rec060_block1/';
load([dest_dir '/state_responses_20250925_nil_CPR_block1_phy4_rec060_fxs.mat'])
load([dest_dir '/summary_20250925_nil_CPR_block1_phy4_rec060_fxs.mat'])

dest_dir        = '/Users/cnl/Documents/DATA/Nilan/plots/20250925_rec060_block1/';
if exist(dest_dir,'dir') ~= 7
    mkdir(dest_dir)
end

%% Solo vs dyadic behavior

addpath /Users/cnl/Desktop/CPR/CircStat2012a/
addpath /Users/cnl/Desktop/CPR/Violinplot-Matlab/

% Avg response: tilt and error
tlt_state = cellfun(@(x) mean(x),state.js_monk_tlt);
for iState = 1:length(state.js_monk_dir)
    clear js_dev
    js_dev = rad2deg(circ_dist(deg2rad(state.js_monk_dir{iState}), deg2rad(state.rdp_dir(iState))));
    err_state(iState) = mean(abs(js_dev));
end

snr = unique(state.rdp_coh);
solo_idx = cellfun(@(x) contains(x,'CPR_solo_stepfunction_neutral'), state.task);
dyad_idx = cellfun(@(x) contains(x,'CPR_dyadic_stepfunction_neutral'), state.task);

for iCoh = 1:length(snr)
    coh_idx = state.rdp_coh == snr(iCoh);
    
    [p_tlt,h_tlt,stats_tlt] = ranksum(tlt_state(coh_idx & solo_idx), tlt_state(coh_idx & dyad_idx))
    [p_err,h_err,stats_err] = ranksum(err_state(coh_idx & solo_idx), err_state(coh_idx & dyad_idx))

    tlt(1,iCoh) = median(tlt_state(coh_idx & solo_idx));
    tlt(2,iCoh) = median(tlt_state(coh_idx & dyad_idx));

    err(1,iCoh) = median(err_state(coh_idx & solo_idx));
    err(2,iCoh) = median(err_state(coh_idx & dyad_idx));
end

x = round(snr.*100);
figure
subplot(1,2,1); hold on
plot(x,tlt(1,:),'r', 'LineWidth',2)
plot(x,tlt(2,:),'k', 'LineWidth',2)
title('tilt')
ylabel('joystick tilt [norm]')
xlabel('coherence')
set(gca, 'XTick', x)
set(gca, 'fontsize', 20)

subplot(1,2,2); hold on
plot(x,err(1,:),'r', 'LineWidth',2)
plot(x,err(2,:),'k', 'LineWidth',2)
title('error')
ylabel('joystick error [deg]')
xlabel('coherence')
set(gca, 'XTick', x)
legend('solo', 'dyad')
set(gca, 'fontsize', 20)

print([dest_dir '/avg_behavior'],'-dpng')
print([dest_dir '/avg_behavior'],'-dsvg')


% clear clus_idx clus_ts
% clus = unique(phy.brain.RF.raw.ch008_neg.spks_id{1});
% clus_idx = cellfun(@(x) x == clus(1), phy.brain.RF.raw.ch015_neg.spks_id, 'UniformOutput', false);
% 
% for iCyc = 1:length(clus_idx)
%     if sum(clus_idx{iCyc}) > 0
%         clus_ts{iCyc} = phy.brain.RF.raw.ch015_neg.spks_ts{iCyc}(clus_idx{iCyc});
%     else
%         clus_ts{iCyc} = [];
%     end
% end

%% List of good units
phy.brain.CPR.spks.include(phy.brain.CPR.spks.include.inclusion_flag,:)
sum(phy.brain.CPR.spks.include.inclusion_flag) / size(phy.brain.CPR.spks.include,1)

%% RF map
plot_RF_map(phy.brain.RF.stim_id, phy.brain.RF.stim_pos, phy.brain.RF.ch046_neg.nSpikes)

%% RF raster
plot_RF_raster(phy.brain.RF.raw.ch046_neg.spks_ts,100);

%% Unit plots
addpath /Users/cnl/Desktop/CPR/CircStat2012a
lst                 = phy.brain.CPR.spks.include.unit_ID(phy.brain.CPR.spks.include.inclusion_flag);
lst(contains(lst, 'unit0')) = [];
snr                 = unique(state.rdp_coh);
stim_coh            = state.rdp_coh;

for iUnit = 1:length(lst)
    close all
    plot_onset_raster(phy,char(lst(iUnit)), dest_dir)
    plot_cpr_tuning(state, char(lst(iUnit)), dest_dir)
    plot_spkjoy_correlation(state,char(lst(iUnit)), dest_dir)
end

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
    ps.MarkerFaceAlpha = .5;
    ps.MarkerEdgeColor = col;
end
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';

nbins = 360/20;
edges = linspace(0, 2*pi, nbins+1);
r_mean_bin = zeros(1, nbins);
theta_bin = zeros(1, nbins);

for i = 1:nbins
    idx = theta >= edges(i) & theta < edges(i+1);
    r_mean_bin(i) = mean(r(idx));
    theta_bin(i) = (edges(i)+edges(i+1))/2;
end

r_mean_bin = smoothdata(r_mean_bin, 'gaussian', 5);

if isempty(col)
    polarplot(ax, [theta_bin theta_bin(1)], [r_mean_bin r_mean_bin(1)], 'r-', 'LineWidth', 2)
else
    polarplot(ax, [theta_bin theta_bin(1)], [r_mean_bin r_mean_bin(1)], 'LineStyle','-', 'LineWidth', 2, 'Color', col)
end

set(gca, 'FontSize', 16)
end

function plot_spkjoy_correlation(state,unit_id,dest_dir)

FR                  = state.spk_n.(unit_id) ./ state.dur_s;
FR_idx              = FR > 0;
tlt                 = cellfun(@mean,state.js_monk_tlt);
snr                 = unique(state.rdp_coh);
col                 = cool(length(snr));

for iState = 1:length(state.dur_s)
    js_dev              = rad2deg(circ_dist(deg2rad(state.js_monk_dir{iState}),deg2rad(state.rdp_dir(iState)))); % Get circular distance to RDP direction
    state.js_acc{iState}= abs(1 - abs(js_dev / 180));                               % Calculate accuracy
end
acc                 = cellfun(@mean,state.js_acc);

figure; hold on
subplot(2,1,1); hold on
% sc = scatter(tlt(FR~=0),FR(FR~=0),'k.');
for iCoh = 1:length(snr)
    coh_idx = state.rdp_coh == snr(iCoh);
    sc = scatter(tlt(coh_idx & FR~=0),FR(coh_idx & FR~=0));
    sc.Marker = '.';
    sc.SizeData = 50;
    sc.MarkerEdgeColor = col(iCoh,:);
    sc.SizeData = 50;

    [pfit] = polyfit(tlt(coh_idx & FR~=0),FR(coh_idx & FR~=0),1);
    yfit = polyval(pfit,[0:.1:1]);
    plot([0:.1:1],yfit,'Color', col(iCoh,:),'LineWidth',2);
    % [r,p] = corrcoef(tlt(coh_idx & FR~=0),FR(coh_idx & FR~=0))
end

[r,p] = corrcoef(tlt(FR~=0),FR(FR~=0));
xlabel('JS tilt [norm]')
ylabel('FR [Hz]')
title(['r: ' num2str(r(2)) ' | p: ' num2str(p(2))])
set(gca,'fontsize',20)

subplot(2,1,2); hold on 
% sc = scatter(acc(FR~=0),FR(FR~=0),'k.');
for iCoh = 1:length(snr)
    coh_idx =  state.rdp_coh == snr(iCoh);
    sc = scatter(acc(coh_idx & FR~=0),FR(coh_idx & FR~=0));
    sc.Marker = '.';
    sc.SizeData = 50;
    sc.MarkerEdgeColor = col(iCoh,:);
    sc.SizeData = 50;
    [pfit] = polyfit(acc(coh_idx & FR~=0),FR(coh_idx & FR~=0),1);
    yfit = polyval(pfit,[0:.1:1]);
    plot([0:.1:1],yfit,'Color', col(iCoh,:),'LineWidth',2);
end

[r,p] = corrcoef(acc(FR~=0),FR(FR~=0));
xlabel('JS accuracy [norm]')
ylabel('FR [Hz]')
title(['r: ' num2str(r(2)) ' | p: ' num2str(p(2))])
set(gca,'fontsize',20)

% for iCoh = 1:length(snr)
%     lsl(iCoh).Color = col(iCoh,:);
%     lsl(iCoh).LineWidth = 2;
% end

print([dest_dir '/corr_spk_joy_' char(unit_id)],'-dpng')
print([dest_dir '/corr_spk_joy_' char(unit_id)],'-dsvg')
end

function plot_onset_raster(phy, unit_id, dest_dir)

str_chan = unit_id(1:end-6);
str_unit = unit_id(end-4:end);
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
ax.FontSize         = 20;

ax                  = subplot(3,3,7:9);
pl                  = plot(mean(sdf),'k', 'LineWidth', 1.5);
ax.XLim             = [0 1000];
ax.XLabel.String    = 'time [ms]';
ax.YLabel.String    = 'FR [Hz]';
ax.Box              = 'off';
ax.FontSize         = 20;

print([dest_dir '/cpr_onset_' unit_id],'-dpng')
print([dest_dir '/cpr_onset_' unit_id],'-dsvg')
end

function plot_cpr_tuning(state, unit_id, dest_dir)


dyad_idx            = contains(state.task,'dyad');
solo_idx            = contains(state.task,'solo');
incl                = state.include.(unit_id);

FR                  = state.spk_n.(unit_id) ./ state.dur_s;
snr                 = unique(state.rdp_coh);
stim_coh            = state.rdp_coh;
stim_dir            = deg2rad(state.rdp_dir);
col                 = cool(length(snr));

f1 = figure;
h = polaraxes; hold on
plot_CPRtuning(h,FR(solo_idx & incl),stim_dir(solo_idx & incl), false,[1 0 0])
plot_CPRtuning(h,FR(dyad_idx & incl),stim_dir(dyad_idx & incl), false,[0 0 1])
legend('solo','dyad')
print([dest_dir '/cpr_tuning_' unit_id],'-dpng')
print([dest_dir '/cpr_tuning_' unit_id],'-dsvg')

% Coherence-based tuning
f2                  = figure('Units','normalized','Position',[.2 .2 .6 .6]);
nrows = 2;
ncols = 2;
for iCoh = 1:length(snr)
    [row,clm] = ind2sub([nrows,ncols], iCoh);
    left = (clm-1)/ncols + 0.08;
    bottom = 1 - row/nrows + 0.05;
    w = 0.4; h = 0.4;
    pax = polaraxes('Position',[left bottom w h]);

    coh_idx         = stim_coh == snr(iCoh);
    plot_CPRtuning(pax,FR(coh_idx & solo_idx & incl),stim_dir(coh_idx & solo_idx & incl),true,col(iCoh,:));
    plot_CPRtuning(pax,FR(coh_idx & dyad_idx & incl),stim_dir(coh_idx & dyad_idx & incl),true,col(iCoh,:)./2);
    lg = legend('','solo','','dyad');

end
% lg = legend(cellfun(@num2str,{snr(1) snr(2) snr(3) snr(4)},'UniformOutput',false));

print([dest_dir '/cpr_coh_tuning_' unit_id],'-dpng')
print([dest_dir '/cpr_coh_tuning_' unit_id],'-dsvg')

end

function rho_smooth = circRunningMean(theta, rho, winWidth, isDegrees)
% CIRCRUNNINGMEAN  Circular running mean (Gaussian) for polar data.
%
%   rho_smooth = circRunningMean(theta, rho, winWidth)
%   rho_smooth = circRunningMean(theta, rho, winWidth, isDegrees)
%
%   Inputs:
%       theta      - angle values (radians or degrees)
%       rho        - corresponding data values
%       winWidth   - window width (same units as theta)
%       isDegrees  - optional, true if theta is in degrees (default: false)
%
%   Output:
%       rho_smooth - circularly smoothed rho (same size as input)
%
%   Example:
%       theta = linspace(0, 2*pi, 200);
%       rho = 1 + 0.3*sin(3*theta) + 0.1*randn(size(theta));
%       rho_smooth = circRunningMean(theta, rho, deg2rad(20));
%       polarplot(theta, rho, ':k'); hold on
%       polarplot(theta, rho_smooth, 'r', 'LineWidth', 1.5);

    if nargin < 4
        isDegrees = false;
    end

    % Convert degrees to radians if needed
    if isDegrees
        theta = deg2rad(theta);
        winWidth = deg2rad(winWidth);
    end

    % Ensure column vectors
    theta = theta(:);
    rho = rho(:);

    % Sort by theta just in case
    [theta, idx] = sort(theta);
    rho = rho(idx);

    % Estimate angular step
    dtheta = mean(diff(theta));

    % Create Gaussian kernel (±3σ range)
    sigma = (winWidth / 2.355) / dtheta; % FWHM to sigma in samples
    x = -3*sigma : 3*sigma;
    kernel = exp(-x.^2 / (2*sigma^2));
    kernel = kernel / sum(kernel);

    % Circular convolution using data wrapping
    rho_ext = [rho; rho; rho];  % triple to wrap
    rho_conv = conv(rho_ext, kernel, 'same');
    n = numel(rho);
    rho_smooth = rho_conv(n+1 : 2*n);  % extract central part

    % Restore original order
    [~, invIdx] = sort(idx);
    rho_smooth = rho_smooth(invIdx);
end
