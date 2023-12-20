
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/

dyad = 51;
block = 1;
cycle = 15;

%% IMPORT %%%

xpth                = pwd;
dpth                = ['/Volumes/T7_Shield/CPR_psychophysics/Dyad' num2str(dyad) '/summary'];
dest_dir            = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/raw';

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
    ts              = [ts (tp1.frme_ts{iState}-ts1)./1e6];                                      % Timestamps
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

%% PLOT %%%

trg_cl              = [.9 .5 .3];

f                	= figure('units','normalized','position',[0 0 .5 .5]);
ax                  = gca; hold on

change_idx       	= [1 find(diff(rdp_coh))+1 length(rdp_coh)];
coh_id              = round(rdp_coh(change_idx(1:3)),2);
ccol                = [.85 .85 .85;.95 .95 .95; .9 .9 .9];
for iCoh = 1:3
%     ccol            = [rdp_coh(change_idx(iCoh)) .8 .8]
    pt(iCoh)           	= patch([ts(change_idx(iCoh)) ts(change_idx(iCoh)) ts(change_idx(iCoh+1)) ts(change_idx(iCoh+1))],[0 360 360 0],ccol(iCoh,:));
    pt(iCoh).EdgeColor  	= 'none';
end

p                   = plot(ts,rdp_dir,'LineWidth',lw*1.5,'Color', [0 0 0]);
sc                  = scatter(trg_ts, trg_dir);
sc.SizeData         = 50;
sc.MarkerFaceColor  = trg_cl;
sc.MarkerEdgeColor  = trg_cl;

ax.YLim             = [0 360];
ax.XLim             = [1 ts(end)];
ax.YTick            = [0 90 180 270 360];
ax.YLabel.String    = {'JS direction';'raw [deg]'};
ax.XLabel.String    = 'Time [s]';
ax.FontSize         = fs;

lg                  = legend(pt,{num2str(coh_id(1)),num2str(coh_id(2)),num2str(coh_id(3))});
lg.Location         = 'south';
lg.Box              = 'on';

print(f, [dest_dir '/rdp_example_patch'], '-r500', '-dpng');

%% Partial response profile

f                   = figure;
ax                  = gca; hold on

% Plot RDP stepfunction
pidx                = 2450:length(ts); % Sample index
p                   = plot(ts(pidx),rdp_dir(pidx),'LineWidth',lw*1.5,'Color', [0 0 0]);

% Use scatter to color-code 2D joystick response
for i = pidx
    if ~isnan(p1_ecc(i))
        sc                  = scatter(ts(i), p1_dir(i));
        sc.MarkerFaceColor  = [p1_ecc(i) 0 0];
        sc.MarkerEdgeColor  = [p1_ecc(i) 0 0];
        sc.SizeData         = 50;
    end
    
    if ~isnan(p2_ecc(i))
        sc                  = scatter(ts(i), p2_dir(i));
        sc.MarkerFaceColor  = [0 p2_ecc(i) 0];
        sc.MarkerEdgeColor  = [0 p2_ecc(i) 0];
        sc.SizeData         = 50;
    end
end

strt         	= 3600;
fs              = 1000/120;
len             = 2000; % 2sec
height          = 225;
ln              = line([ts(strt) ts(strt)+(len/fs)/100],[height height]);
ln.Color        = [0 0 0];
ln.LineWidth  	= 2;

ax.YAxis.Visible =  'off';
ax.XAxis.Visible =  'off';

print(f, [dest_dir '/js_example_response'], '-r500', '-dpng');
exportgraphics(f,[dest_dir '/js_example_response.pdf'],'ContentType','vector') 