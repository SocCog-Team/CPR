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

str_chan            = 'ch012_neg';
str_unit            = 'unit1';
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

%% CPR tuning curve based on state-wise firing rate

incl                = state.include.('ch012_neg_unit2');
nSpks               = cellfun(@(x) x< 1e6, state.spk_ts.('ch012_neg_unit2'), 'UniformOutput', false);
FR                  = cellfun(@sum,nSpks);
stim_dir            = deg2rad(state.rdp_dir);

f = figure;
h = polaraxes; hold on
plot_CPRtuning(h,FR(incl),stim_dir(incl), true,[0 0 0])

%% CPR tuning: coherence-wise and context-wise
plot_cpr_tuning(state, 'ch012_neg_unit2', [])

%% CPR cycle onset raster.

str_chan            = 'ch012_neg';
str_unit            = 'unit1';
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

%% Solo vs dyadic SDF - cycle onset
phy.brain.CPR.spks.include(phy.brain.CPR.spks.include.inclusion_flag,:)

str_chan            = 'ch043_neg';
str_unit            = 'unit3';
cpr_spk_times       = [];
for iCyc = 1:length(phy.stim.cpr_cyle)
    cpr_spk_times{iCyc} = double(phy.brain.CPR.spks.(str_chan).(str_unit){iCyc});% ./1e6;
    cpr_spk_times{iCyc}(cpr_spk_times{iCyc} > 1e6) = [];
end

cyc_idx             = phy.brain.CPR.spks.include.cyc_id{phy.brain.CPR.spks.include.unit_ID == [str_chan '_' str_unit]};
solo_idx            = cell2mat(phy.stim.cpr_solo);
tstep               = .001;
time                = [0:tstep:1];
[~, sdf_solo]       = FR_estimation(cpr_spk_times(cyc_idx' & solo_idx), time, false);
[~, sdf_dyad]       = FR_estimation(cpr_spk_times(cyc_idx' & ~solo_idx), time, false);

f                   = figure; hold on
pls                 = plot(mean(sdf_solo),'k', 'LineWidth', 1.5);
pld                 = plot(mean(sdf_dyad),'r', 'LineWidth', 1.5);
ax                	= gca;
ax.XLim             = [0 1000];
ax.XLabel.String    = 'time [ms]';
ax.YLabel.String    = 'FR [Hz]';
ax.Box              = 'off';
ax.FontSize         = 16;
ax.Title.String     = 'Example unit: Cycle onset';
lg                  = legend([pls pld], 'solo', 'dyad');

%% Solo vs dyadic SDF - state onset
close all
str_chan            = 'ch015_neg';
str_unit            = 'unit2';

state_idx           = state.include.([str_chan '_' str_unit]);
solo_idx            = state.cpr_solo;

bin_width           = 90;
[PD, ~,VS]          = preferredDirection(state.rdp_dir(state_idx), FR(state_idx));
roi                 = mod([PD-(bin_width/2) PD+(bin_width/2)],360);

if roi(1) > roi(2)
    PD_bin          = state.rdp_dir > roi(1) | state.rdp_dir < roi(2);
else
    PD_bin          = state.rdp_dir > roi(1) & state.rdp_dir < roi(2);
end

% Exclude cycle onset states
excl_cycOn              = [1 (find(diff(state.cIdx)))+1];
state_idx(excl_cycOn)   = 0;
solo_idx(excl_cycOn)    = 0;
% Exclude short states
excl_dur                = state.dur_s < 1.5;
state_idx(excl_dur)     = 0;
solo_idx(excl_dur)      = 0;

tstep               = .001;
time                = [0:tstep:1.5];

figure; hold on
[~, sdf_solo]       = FR_estimation(state.spk_ts.([str_chan '_' str_unit])(state_idx & solo_idx & PD_bin), time, true);
figure; hold on
[~, sdf_dyad]       = FR_estimation(state.spk_ts.([str_chan '_' str_unit])(state_idx & ~solo_idx & PD_bin), time, true);

f                   = figure; hold on
pls                 = plot(mean(sdf_solo),'k', 'LineWidth', 1.5);
pld                 = plot(mean(sdf_dyad),'r', 'LineWidth', 1.5);
ax                	= gca;
ax.XLim             = [0 1500];
ax.XLabel.String    = 'time [ms]';
ax.YLabel.String    = 'FR [Hz]';
ax.Box              = 'off';
ax.FontSize         = 16;
ax.Title.String     = 'Example unit: State onset';
lg                  = legend([pls pld], 'solo', 'dyad');

%% CPR spike - joystick correlation
% addpath /Users/cnl/Desktop/CPR/CircStat2012a
% 
% str_chan            = 'ch012_neg';
% str_unit            = 'unit2';
% snr                 = unique(state.rdp_coh);
% stim_coh            = state.rdp_coh;
% snr(snr==0)         = [];
% col                 = cool(3);
% FR                  = state.spk_n.([str_chan '_' str_unit]) ./ state.dur_s;
% tlt                 = cellfun(@mean,state.js_monk_tlt);
% 
% for iState = 1:length(state.dur_s)
%     js_dev              = rad2deg(circ_dist(deg2rad(state.js_monk_dir{iState}),deg2rad(state.rdp_dir(iState)))); % Get circular distance to RDP direction
%     state.js_acc{iState}= abs(1 - abs(js_dev / 180));                               % Calculate accuracy
% end
% acc                 = cellfun(@mean,state.js_acc);
% 
% figure; hold on
% sc = scatter(tlt(FR~=0),FR(FR~=0),'k.');
% % sc = scatter(acc(FR~=0),FR(FR~=0),'k.');
% sc.SizeData = 50;
% [r,p] = corrcoef(tlt(FR~=0),FR(FR~=0))
% % [r,p] = corrcoef(acc(FR~=0),FR(FR~=0))
% lsl = lsline;
% lsl.Color = 'r';
% lsl.LineWidth = 2;
% 
% % for iCoh = 1:length(snr)
% %     coh_idx         = stim_coh == snr(iCoh);
% %     sc = scatter(tlt(coh_idx & FR~=0),FR(coh_idx & FR~=0));
% %     sc.Marker = '.';
% %     sc.SizeData = 50;
% %     sc.MarkerEdgeColor = col(iCoh,:);
% %     [r,p] = corrcoef(tlt(coh_idx & FR~=0),FR(coh_idx & FR~=0))
% %     lsl = lsline;
% % end
% %
% % for iCoh = 1:length(snr)
% %     lsl(iCoh).Color = col(iCoh,:);
% %     lsl(iCoh).LineWidth = 2;
% % end
% 
% xlabel('JS tilt [norm]')
% ylabel('FR [Hz]')
% title(['r: ' num2str(r(2)) ' | p: ' num2str(p(2))])

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

function [prefDir, prefMag, vecStrength] = preferredDirection(angles, rates)
% preferredDirection computes the resultant vector and tuning strength from polar data.
%
%   [prefDir, prefMag, vecStrength] = preferredDirection(angles, rates)
%
%   INPUTS:
%       angles - vector of stimulus directions (in degrees or radians)
%       rates  - vector of corresponding firing rates
%
%   OUTPUTS:
%       prefDir     - preferred direction (same unit as input angle)
%       prefMag     - magnitude of the resultant vector
%       vecStrength - normalized vector magnitude (0–1), i.e. direction selectivity
%
%   Example:
%       angles = 0:45:315;
%       rates  = [5 8 12 9 4 3 2 6];
%       [dir, mag, vs] = preferredDirection(angles, rates)
%
%   See also: atan2, deg2rad, rad2deg

% Check input size
if numel(angles) ~= numel(rates)
    error('angles and rates must have the same length.');
end

% Detect whether angles are in degrees or radians
if max(abs(angles)) > 2*pi
    angRad = deg2rad(angles);
    useDegrees = true;
else
    angRad = angles;
    useDegrees = false;
end

% Compute resultant vector components
x = sum(rates .* cos(angRad));
y = sum(rates .* sin(angRad));

% Resultant vector magnitude
prefMag = sqrt(x^2 + y^2);

% Preferred direction (angle of resultant)
prefDir = atan2(y, x);

% Convert to degrees if needed
if useDegrees
    prefDir = rad2deg(prefDir);
    prefDir = mod(prefDir, 360);
end

% Normalized vector strength (0–1)
vecStrength = prefMag / sum(rates);
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
set(gca, 'GridAlpha', .2)
end

function plot_cpr_tuning(state, unit_id, dest_dir)

solo_idx            = state.cpr_solo;
incl                = state.include.(unit_id);

%%% full state %%%
FR                  = state.spk_n.(unit_id) ./ state.dur_s;
%%% first 500 ms %%%
% nSpks               = cellfun(@(x) x<500e3, state.spk_ts.(unit_id), 'UniformOutput', false);
% FR                  = cellfun(@sum,nSpks)./0.5;

%%% last 500ms %%%
% thresh = (state.dur_s-0.5).*1e6;
% for iS = 1:length(state.spk_ts.(unit_id))
%     if isempty(state.spk_ts.(unit_id){iS})
%         nSpikes(iS) = 0;
%     else
%         nSpks(iS) = sum(state.spk_ts.(unit_id){iS} > thresh(iS));
%     end
% end
% FR                  = nSpks./0.5;

snr                 = unique(state.rdp_coh);
stim_coh            = state.rdp_coh;
stim_dir            = deg2rad(state.rdp_dir);
col                 = cool(length(snr));

f1 = figure;
h = polaraxes; hold on
plot_CPRtuning(h,FR(solo_idx & incl),stim_dir(solo_idx & incl), false,[1 0 0])
plot_CPRtuning(h,FR(~solo_idx & incl),stim_dir(~solo_idx & incl), false,[0 0 1])
legend('solo','dyad')
% print([dest_dir '/cpr_tuning_' unit_id],'-dsvg')

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
    plot_CPRtuning(pax,FR(coh_idx & ~solo_idx & incl),stim_dir(coh_idx & ~solo_idx & incl),true,col(iCoh,:)./2);
    lg = legend('','solo','','dyad');

end
lg = legend(cellfun(@num2str,{snr(1) snr(2) snr(3) snr(4)},'UniformOutput',false));
% print([dest_dir '/cpr_coh_tuning_' unit_id],'-dsvg')

f3                  = figure('Units','normalized','Position',[.2 .2 .6 .6]);
h = polaraxes; hold on
for iCoh = 1:length(snr)
    coh_idx         = stim_coh == snr(iCoh);
    plot_CPRtuning(gca,FR(coh_idx & incl),stim_dir(coh_idx & incl),false,col(iCoh,:));
end
% print([dest_dir '/cpr_coh_tuning__comb' unit_id],'-dsvg')

end
