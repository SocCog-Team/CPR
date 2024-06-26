addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/

dyad = 51;
block = 1;
cycle = 15;

%% IMPORT %%%

xpth                = pwd;
dpth                = ['/Volumes/T7_Shield/CPR_psychophysics/Dyad' num2str(dyad) '/summary'];
dest_dir            = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG_task/raw/';

cd(dpth)
mat_files           = dir('*.mat');

p1                  = [];
for iFile = 1:2
    tbl             = load(mat_files(iFile).name);
    p1              = [p1; tbl.t];
end

p2                  = [];
for iFile = 3:4
    tbl          	= load(mat_files(iFile).name);
    p2              = [p2; tbl.t];
end

% Extract cycle data
idx                 = p1.cyc_no == cycle & p1.block == block;
tp1                 = p1(idx,:);
tp2                 = p2(idx,:);
ts1                 = tp1.frme_ts{1}(1);
cl                  = {[.9 .5 .3],[.3 .6 .75]};
fs                  = 18;
lw                  = 2;

% Initialise variables
ts                  = [];
trg_ts              = [];
trg_dir             = [];
rdp_coh             = [];
rdp_dir             = [];
p1_dir              = [];
p2_dir              = [];
p1_ecc              = [];
p2_ecc              = [];

% Append to vector
for iState = 1:size(tp1,1)
    ts              = [ts double((tp1.frme_ts{iState}-ts1))./1e6];                                      % Timestamps
    trg_ts          = [trg_ts (tp1.trg_ts{iState}-ts1)./1e6];                                   % Target timestamps
    trg_dir         = [trg_dir repmat(tp1.rdp_dir(iState),1,length(tp1.trg_ts{iState}))];       % Target direction
    
    rdp_coh         = [rdp_coh repmat(tp1.rdp_coh(iState),1,length(tp1.frme_ts{iState}))];
    rdp_dir         = [rdp_dir repmat(tp1.rdp_dir(iState),1,length(tp1.frme_ts{iState}))];
    p1_dir          = [p1_dir tp1.js_dir{iState}];
    p2_dir          = [p2_dir tp2.js_dir{iState}];
    p1_ecc          = [p1_ecc tp1.js_ecc{iState}];
    p2_ecc          = [p2_ecc tp2.js_ecc{iState}];
end

trg_dir(isnan(trg_ts)) = [];
trg_ts(isnan(trg_ts)) = [];

%% Partial response profile

f                	= figure;
ax                  = gca; hold on

time_bar                    = line([25 27],[130 130]); % 2s
time_bar.LineStyle          = '-';
time_bar.LineWidth          = 2;
time_bar.Color              = [0 0 0];

% Plot RDP stepfunction
pidx                        = 2450:6000; % Sample index
p                           = plot(ts(pidx),rdp_dir(pidx),'LineWidth',lw*1.5,'Color', [0 0 0]);
lum                         = 0;
alp                         = .75;

% Plot joystick response
sc1                         = scatter(ts(pidx), p1_dir(pidx),'filled');
sc1.CData                   = [p1_ecc(pidx)' zeros(length(ts(pidx)),1)+lum zeros(length(ts(pidx)),1)+lum];
sc1.MarkerFaceAlpha         = alp;
sc1.MarkerEdgeAlpha         = alp;

sc2                         = scatter(ts(pidx), p2_dir(pidx),'filled');
sc2.CData                   = [zeros(length(ts(pidx)),1)+lum p2_ecc(pidx)' zeros(length(ts(pidx)),1)+lum];
sc2.MarkerFaceAlpha         = alp;
sc2.MarkerEdgeAlpha         = alp;

ax.YAxis.Visible            =  'off';
ax.XAxis.Visible            =  'off';

print(f, [dest_dir '/js_example_response'], '-r500', '-dpng');
print(f, [dest_dir '/js_example_response'], '-r500', '-dsvg');

%% Coherence legend
 
% f               = figure;
% cb              = colorbar;
% cb.Ticks        = [];
% cb.Position     = [.4 .1 .2 .8];
% cb.Box          = 'off';
% cl              = linspace(0,.9,256)';
% 
% axis off
% colormap([cl cl cl])
% 
% print(f, [dest_dir '/lgnd_coh'], '-r300', '-dpng');

%% PLOT DIR+COH combined %%%

% trg_cl              = [.9 .5 .3];
% 
% f                	= figure('units','normalized','position',[0 0 .5 .5]);
% ax                  = gca; hold on
% p                   = plot(ts,rdp_dir,'LineWidth',lw*1.5,'Color', [0 0 0]);
% sc                  = scatter(trg_ts, trg_dir);
% sc.SizeData         = 50;
% sc.MarkerFaceColor  = trg_cl;
% sc.MarkerEdgeColor  = trg_cl;
% 
% ax.YLim             = [0 360];
% ax.XLim             = [1 ts(end)];
% ax.YTick            = [0 90 180 270 360];
% ax.YLabel.String    = {'JS direction';'raw [deg]'};
% ax.XLabel.String    = 'Time [s]';
% ax.FontSize         = fs;
% 
% yyaxis right
% coh_cl              = [.5 .5 .5];
% p                   = plot(ts,rdp_coh,'LineStyle', '--','LineWidth',lw*1.5,'Color', coh_cl);
% ax.YLabel.String    = {'RDP coherence';'[%]'};
% ax.YTick            = [0 .25 .5 .75 1];
% ax.YAxis(2).Color   = coh_cl;
% 
% print(f, [dest_dir '/rdp_example_line'], '-r500', '-dpng');
% print(f, [dest_dir '/rdp_example_line'], '-r500', '-dsvg');
