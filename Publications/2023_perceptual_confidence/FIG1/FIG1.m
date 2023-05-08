%% FIG1 RDP ELEMENTS %%
% State sequence
f               = figure;
dest_dir        = '/Users/fschneider/ownCloud/Documents/Publications/CPR_psychophysics/Figures/FIG1/raw';

dur             = [0 randi([125 250],1,10)];
height          = 100;
fixcol          = 0;
maxcol          = .85;
cl              = linspace(maxcol,.1,10);
cl            	= cl(randperm(length(cl)));

% Direction
for j = 1:length(dur)-1
x               = [sum(dur(1:j)) sum(dur(1:j+1)) sum(dur(1:j+1)) sum(dur(1:j))];
y             	= [0 0 height height];
p               = patch(x,y,[cl(j) cl(j) cl(j)]);
p.EdgeColor     = 'none';
end

% Coherence
rnge            = [.24 .45 .31];
cdur            = [0 sum(dur)*rnge(1) sum(dur)*rnge(2) sum(dur)*rnge(3)];
m               = 5;

maxcol          = .75;
fixcol          = .2;
cl              = linspace(maxcol,0,10);
cl            	= cl(randperm(length(cl)));

for j = 1:length(cdur)-1
x               = [sum(cdur(1:j)) sum(cdur(1:j+1)) sum(cdur(1:j+1)) sum(cdur(1:j))];
y             	= [height*(m-1) height*(m-1) height*m height*m];
p               = patch(x,y,[fixcol cl(j) cl(j)]);
p.EdgeColor     = 'none';
end

% Cycle
mm              = 9;
xmat            = [0 sum(dur)*.2 sum(dur)*.2 0;...
                    sum(dur)*.25 sum(dur)*.75 sum(dur)*.75 sum(dur)*.25; ...
                    sum(dur)*.8 sum(dur) sum(dur) sum(dur)*.8];
y             	= [height*(mm-1) height*(mm-1) height*mm height*mm];

for i = 1:size(xmat,1)
p               = patch(xmat(i,:),y,[1 1 1]);
p.EdgeColor     = [0 0 0];
p.LineWidth     = 2;
end

% Connection lines
l1 = line([0 sum(dur)*rnge(1)],[height height*(m-1)]);
l1.LineStyle = '--';
l1.LineWidth = 2;
l1.Color = 'k';
l2 = line([sum(dur)*sum(rnge(1:2)) sum(dur)],[height*(m-1) height]);
l2.LineStyle = '--';
l2.LineWidth = 2;
l2.Color = 'k';

l3 = line([0 sum(dur)*.25],[height*m height*(mm-1)]);
l3.LineStyle = '--';
l3.LineWidth = 2;
l3.Color = 'k';
l4 = line([sum(dur)*.75 sum(dur)],[height*(mm-1) height*m]);
l4.LineStyle = '--';
l4.LineWidth = 2;
l4.Color = 'k';

box off
axis equal
set(gca,'visible','off')
print(f, [dest_dir '/state_sequence'], '-r300', '-dpng');

%% LEGENDS

% Direction legend
f               = figure;
maxcol          = .85;
fixcol          = 0;
cl              = linspace(maxcol,0,361);
scle            = 100;
ofs             = 0;
for iDir = 1:360
    theta       = linspace(iDir, iDir+1)/180*pi;
    p1          = patch([0 cos(theta).*scle 0]+ofs, [0 sin(theta).*scle 0]+height/2, [cl(iDir+1) cl(iDir+1) cl(iDir+1)]);
    p1.EdgeColor= 'none';
end

x               = 1.5;
theta           = linspace(0, 360)/180*pi;
p               = patch(([0 cos(theta).*scle 0]/x)+ofs, ([0 sin(theta).*scle 0]/x)+height/2, [.99999 .99999 .99999]);
p.EdgeColor     = 'none';

axis equal
axis off
print(f, [dest_dir '/lgnd_dir'], '-r300', '-dpng');

% Coherence legend
f               = figure;
cb              = colorbar;
cb.Ticks        = [];
cb.Position     = [.4 .1 .2 .8];
cb.Box          = 'off';
cl              = linspace(0,maxcol,256)';

axis off
colormap([repmat(fixcol,length(cl),1) cl cl])

print(f, [dest_dir '/lgnd_coh'], '-r300', '-dpng');

%% State sequence - actual data

load('/Users/fschneider/Documents/CPR_psychophysics/Dyad14/summary/20220325_LaS_CPRdyadic_block2_tbl.mat')
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a

trl             = [1 31 61];
lw              = 2;
col             = [0 0 0];
c_bkg           = .7;
fs              = 16;
ts1             = t.trl_frme_ts{1}(1);
ts_end          = t.trl_frme_ts{61}(end);
ax_ofs          = .04;

% Targets
f = figure;
f.Color = [1 1 1];
ax3 = subplot(3,1,3); hold on
maxh = 1;
minh = -.01;
for iTrl = 1:length(trl)
    ts_on(iTrl) = (t.trl_frme_ts{trl(iTrl)}(1) - ts1)  ./ 1e6;
    ts_off(iTrl) = (t.trl_frme_ts{trl(iTrl)}(end) - ts1)  ./ 1e6;
    
    p = patch([ts_on(iTrl) ts_on(iTrl) ts_off(iTrl) ts_off(iTrl)],[minh maxh maxh minh],col+c_bkg);
    p.EdgeColor = 'none';
end

tts = [];
for iState = 1:90
    if t.trg_shown(iState) == 1
        tts = [tts (t.trg_ts{iState}-ts1)./1e6];
    end
end

tl = [ones(1,length(tts)); zeros(1,length(tts)); nan(1,length(tts))];
p = plot(repmat(tts,3,1), tl,'Color', col, 'LineWidth', lw/2);
ax3.XLim = [0 (ts_end-ts1)./1e6];
ax3.XLabel.String = 'Time [s]';
ax3.YTick = [0 1];
ax3.YTickLabel = {'off','on'};
ax3.YLabel.String = 'Target';
ax3.FontSize = fs;
ax3.Position(2) = ax3.Position(2)+ax_ofs*2;

% Direction
ax2 = subplot(3,1,2); hold on
maxh = 361;
minh = -.01;
for iTrl = 1:length(trl)
    ts_on(iTrl) = (t.trl_frme_ts{trl(iTrl)}(1) - ts1)  ./ 1e6;
    ts_off(iTrl) = (t.trl_frme_ts{trl(iTrl)}(end) - ts1)  ./ 1e6;
    
    p = patch([ts_on(iTrl) ts_on(iTrl) ts_off(iTrl) ts_off(iTrl)],[minh maxh maxh minh],col+c_bkg);
    p.EdgeColor = 'none';
end

for iTrl = 1:length(trl)
    p = plot((t.trl_frme_ts{trl(iTrl)}-ts1)./1e6, t.trl_rdp_dir{trl(iTrl)});
    p.Color = col;
    p.LineWidth = lw;
end
ax2.XLim = [0 (ts_end-ts1)./1e6];
ax2.XAxis.Visible = 'off';
ax2.YTick = [0 180 360];
ax2.YLim = [0 361];
ax2.YLabel.String = {'Direction'; '[deg]'};
ax2.FontSize = fs;
ax2.Position(2) = ax2.Position(2)+ax_ofs;

% Cycle/Coherence
ax1 = subplot(3,1,1); hold on
maxh = 1;
minh = -.01;
for iTrl = 1:length(trl)
    ts_on(iTrl) = (t.trl_frme_ts{trl(iTrl)}(1) - ts1)  ./ 1e6;
    ts_off(iTrl) = (t.trl_frme_ts{trl(iTrl)}(end) - ts1)  ./ 1e6;
    
    p = patch([ts_on(iTrl) ts_on(iTrl) ts_off(iTrl) ts_off(iTrl)],[minh maxh maxh minh],col+c_bkg);
    p.EdgeColor = 'none';
end

for iTrl = 1:length(trl)
    p = plot((t.trl_frme_ts{trl(iTrl)}-ts1)  ./ 1e6, t.trl_rdp_coh{trl(iTrl)});
    p.Color = col;
    p.LineWidth = lw;
    
end
ax1.YLim = [minh .81];
ax1.XLim = [0 (ts_end-ts1) ./ 1e6];
ax1.YTick = [0 .4 .8];
ax1.YLabel.String = {'Coherence'; '[%]'};
ax1.FontSize = fs;
ax1.XAxis.Visible = 'off';


axb1 = axes('Position', [0.13 0.19 0.2405 0.73], 'Color',[col+c_bkg]);
axb1.YAxis.Visible = 'off';
axb1.XAxis.Visible = 'off';
uistack(axb1, 'bottom')

axb2 = axes('Position', [0.3835 0.19 0.2615 0.73], 'Color',[col+c_bkg]);
axb2.YAxis.Visible = 'off';
axb2.XAxis.Visible = 'off';
uistack(axb2, 'bottom')

axb3 = axes('Position', [0.6595 0.19 0.245 0.73], 'Color',[col+c_bkg]);
axb3.YAxis.Visible = 'off';
axb3.XAxis.Visible = 'off';
uistack(axb3, 'bottom')

addpath /Users/fschneider/Documents/MATLAB/altmany-export_fig-d7671fe
export_fig(f, [dest_dir '/state_sequence_data'], '-r300', '-dpng')
%print(f, [dest_dir '/state_sequence_data'], '-r300', '-dpng');

%% RDP - Joystick example

trl         = 31;
js          = [];
rdp         = [];
str         = [];

for i = trl:trl+10
js          = [js t.js_dir{i}];
str         = [str t.js_str{i}];
rdp         = [rdp repmat(t.rdp_dir(i),1,length(t.js_dir{i}))];
end

str(isnan(str)) = 0;

for j = 1:length(js)-1
d1          = deg2rad(js(j));
d2          = deg2rad(js(j+1));
df_js(j)    = rad2deg(circ_dist(d1,d2));

d1          = deg2rad(rdp(j));
d2          = deg2rad(rdp(j+1));
df_rdp(j)   = rad2deg(circ_dist(d1,d2));
end

df_js(isnan(df_js)) = 0;
cs_js       = cumsum(df_js) - js(1);
cs_rdp      = cumsum(df_rdp) - rdp(1);
idx         = 500:2500;

% figure
% plot(cs_rdp(idx), 'k')
% hold on
% plot(cs_js(idx),'r')
% box off
% axis off

f           = figure;
plt         = plot(cs_rdp(idx), 'k', 'LineWidth', 5);
hold on

for i = idx
sc                  = scatter(i-(idx(1)-1), cs_js(i));
sc.MarkerFaceColor  = [str(i) 0 0];
sc.MarkerEdgeColor  = [str(i) 0 0];
sc.SizeData         = 50;
end

box off
axis off

print(f, [dest_dir '/js_rdp_trial'], '-r300', '-dpng');

% Eccentricity Legend
f               = figure;
cb              = colorbar;
cb.FontSize     = 24;
cb.Ticks        = [];
cb.Position     = [.4 .1 .2 .8];
cb.Box          = 'off';
colormap([linspace(0,1,256)',zeros(256,1),zeros(256,1)])
axis off
box off

print(f, [dest_dir '/lgnd_ecc'], '-r300', '-dpng');
