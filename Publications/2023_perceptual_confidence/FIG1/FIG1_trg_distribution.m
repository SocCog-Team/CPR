f                           = figure('units','normalized','position',[0 0 .4 .635]);
lb_fs                       = 20;
lw                          = 3;
cmap                        = [0 0 0; gray(256)];
cmap                        = flipud(cmap);
pth                         = '/Volumes/T7_Shield/CPR_psychophysics/';

for iSubj = 1:2
    
    if iSubj == 1
        tbl             	= load([pth '/Dyad22/summary/20220510_anm_CPRdyadic_block1_tbl.mat']);
    else
        tbl             	= load([pth '/Dyad22/summary/20220510_akn_CPRdyadic_block1_tbl.mat']);
    end
    
    c                       = 0;
    t                    	= tbl.t;
    
    for iState = 1:size(t,1)
        if t.trg_shown(iState) == false
            continue
        end
        
        for iTarget = 1:length(t.trg_ts{iState})
            c                   = c+1;
            trg_acc{iSubj}(c)  	= t.trg_acc{iState}(iTarget);
            trg_conf{iSubj}(c)	= t.trg_ecc{iState}(iTarget);
        end
    end
end

% Calculate reward matrix
acc                         = 0:.001:1;
conf                        = 0:.001:1;
rew                         = acc' .* conf;

% Determine arc width for each confidence level
for j = 1:length(conf)
    arc(j)                  = 180 - (180 * conf(j));
end

% Cap arc width at target width (2dva == 12.7587deg at chosen position)
trg_idx                  	= arc < 12.7587;
arc(trg_idx)                = 12.7587;

% For each confidence level, calculate minimum accuracy required to hit
% the target at given arc width - normalised values
hit_width_acc               = 1 - ((arc/2) / 180);
hit_width_acc(trg_idx)   	= 1 - (12.7587/2)/180; % arc width fixed

% Remove position that cannot yield reward from reward matrix
for iAcc = 1:length(acc)
    indx                    = conf < hit_width_acc(iAcc);
    rew(iAcc,indx)          = nan;
end

% Plot reward matrix
ax0                         = gca; hold on
im                          = imagesc(acc,conf,rew);
ax0.XLabel.String           = 'Accuracy';
ax0.YLabel.String           = 'Eccentricity';
ax0.FontSize                = lb_fs;
ax0.XLim                    = [0 1];
ax0.YLim                    = [0 1];
ax0.XLim                    = [0 1];
ax0.XTick                   = [0:.2:1];
ax0.YTick                   = [0:.2:1];
ax0.XTickLabelRotation      = 0;

colormap(cmap)
cb                          = colorbar;
cb.Label.String             = 'Score';
cb.Location                 = 'northoutside';

% Add targets
sc1                         = scatter(trg_acc{1}, trg_conf{1}, 'filled');
sc1.CData                   = [.69 0 0];
sc1.SizeData                = 40;
sc1.MarkerFaceAlpha         = .75;

sc2                         = scatter(trg_acc{2}, trg_conf{2}, 'filled');
sc2.CData                   = [0 .69 0];
sc2.SizeData                = 40;
sc2.MarkerFaceAlpha         = .75;

ax0.Position                = [.25 .19 .65 .65];

% Add histograms on both axes
nBin                        = 25;
vofs                        = .175;
hofs                        = .175;

axh                        = axes('Position', [ax0.Position(1)-vofs ax0.Position(2) ax0.Position(3)/10 ax0.Position(4)]); hold on
for iSubj = 1:2
    [h, edg]            	= histcounts(trg_conf{iSubj},nBin);
    cntr                 	= edg(1:end-1) + diff(edg) ./ 2;
    st                     	= stairs(-h,cntr);
    st.LineWidth          	= lw;
    if iSubj == 1
        st.Color          	= [.65 0 0];
    else
        st.Color          	= [0 .65 0];
    end
    axh.YLim              	= [0 1];
    axh.XAxis.Visible     	= 'off';
    axh.YAxis.Visible       = 'off';
end

axv                         = axes('Position', [ax0.Position(1) ax0.Position(2)-vofs ax0.Position(3) ax0.Position(4)/10]); hold on
for iSubj = 1:2
    [v, edg]              	= histcounts(trg_acc{iSubj},nBin);
    cntr                   	= edg(1:end-1) + diff(edg) ./ 2;
    st                     	= stairs(cntr,-v);
    st.LineWidth          	= lw;
    if iSubj == 1
        st.Color          	= [.65 0 0];
    else
        st.Color          	= [0 .65 0];
    end
    axv.XLim              	= [0 1];
    axv.XAxis.Visible     	= 'off';
    axv.YAxis.Visible    	= 'off';
end

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/raw';

print(f, [dest_dir '/target_distribution_dyad'], '-r500', '-dpng');
print(f, [dest_dir '/target_distribution_dyad'], '-r500', '-dsvg');

%% Histograms: Target count

hist_col                    = [.3 .3 .3];
hist_alpha                  = 1;
nBin                        = 12;

f                           = figure;
tcount                      = cellfun(@length,t.trg_ts);
tcount(t.trg_shown == false) = 0;
[N,EDGES]                   = histcounts(tcount,5,'Normalization','probability');
h                           = histogram('BinEdges', EDGES,'BinCounts',N);
h.FaceColor                 = hist_col;
h.EdgeColor                 = hist_col;
h.FaceAlpha                 = hist_alpha;
h.EdgeAlpha                 = hist_alpha;
ax                          = gca;
ax.XTick                    = 0:4;
ax.XLabel.String            = 'Targets per state';
ax.YLabel.String            = 'Probability [%]';
ax.FontSize                 = lb_fs;
ax.Box                      = 'off';

print(f, [dest_dir '/histogram_target_count'], '-r500', '-dpng');
print(f, [dest_dir '/histogram_target_count'], '-r500', '-dsvg');

%% Histogram: Target interval
f                           = figure;
cnt                         = 0;
state_on                    = cellfun(@(x) x(1),t.frme_ts);

for iState = 1:length(state_on)
    for iTrg = 1:length(t.trg_ts{iState})
        cnt                 = cnt+1;
        trg_ts(cnt)      	= t.trg_ts{iState}(iTrg);
    end
end

trg_ts(isnan(trg_ts))       = [];
iti                         = diff(trg_ts)./1e3;

[N,EDGES]                   = histcounts(iti,nBin,'Normalization','probability');
h                           = histogram('BinEdges', EDGES,'BinCounts',N);
h.FaceColor                 = hist_col;
h.EdgeColor                 = hist_col;
h.FaceAlpha                 = hist_alpha;
h.EdgeAlpha                 = hist_alpha;
ax                          = gca;
ax.XLabel.String            = 'Inter-target interval [ms]';
ax.YLabel.String            = 'Probability [%]';
ax.FontSize                 = lb_fs;
ax.Box                      = 'off';

print(f, [dest_dir '/histogram_inter_target_interval'], '-r500', '-dpng');
print(f, [dest_dir '/histogram_inter_target_interval'], '-r500', '-dsvg');

% Histograms: Target time
f                           = figure;
cnt                         = 0;

for iState = 1:length(state_on)
    for iTrg = 1:length(t.trg_ts{iState})
        cnt             = cnt+1;
        tcount(cnt)       	= t.trg_ts{iState}(iTrg) - state_on(iState);
    end
end

[N,EDGES]                   = histcounts(tcount./1e3,nBin,'Normalization','probability');
h                           = histogram('BinEdges', EDGES,'BinCounts',N);
h.FaceColor                 = hist_col;
h.EdgeColor                 = hist_col;
h.FaceAlpha                 = hist_alpha;
h.EdgeAlpha                 = hist_alpha;
ax                          = gca;
ax.XLabel.String            = 'Time after direction change [ms]';
ax.YLabel.String            = 'Probability [%]';
ax.FontSize                 = lb_fs;
ax.Box                      = 'off';

print(f, [dest_dir '/histogram_target_time'], '-r500', '-dpng');
print(f, [dest_dir '/histogram_target_time'], '-r500', '-dsvg');

% Histograms: State duration
f                           = figure;
[N,EDGES]                   = histcounts(t.ss_dur,nBin,'Normalization','probability');
h                           = histogram('BinEdges', EDGES,'BinCounts',N);
h.FaceColor                 = hist_col;
h.EdgeColor                 = hist_col;
h.FaceAlpha                 = hist_alpha;
h.EdgeAlpha                 = hist_alpha;
ax                          = gca;
ax.XLabel.String            = 'State duration [ms]';
ax.YLabel.String            = 'Probability [%]';
ax.FontSize                 = lb_fs;
ax.Box                      = 'off';

print(f, [dest_dir '/histogram_state duration'], '-r500', '-dpng');
print(f, [dest_dir '/histogram_state_duration'], '-r500', '-dsvg');


% Histograms: Cycle duration
f                           = figure;
[N,EDGES]                   = histcounts(t.cyc_dur./1e3,nBin,'Normalization','probability');
h                           = histogram('BinEdges', EDGES,'BinCounts',N);
h.FaceColor                 = hist_col;
h.EdgeColor                 = hist_col;
h.FaceAlpha                 = hist_alpha;
h.EdgeAlpha                 = hist_alpha;
ax                          = gca;
ax.XLabel.String            = 'Cycle duration [s]';
ax.YLabel.String            = 'Probability [%]';
ax.FontSize                 = lb_fs;
ax.Box                      = 'off';

print(f, [dest_dir '/histogram_cyce duration'], '-r500', '-dpng');
print(f, [dest_dir '/histogram_cycle_duration'], '-r500', '-dsvg');