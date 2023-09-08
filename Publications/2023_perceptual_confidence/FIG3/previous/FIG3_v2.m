% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/WIP/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/comp_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/comp_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/dyad_human_comp_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/dyad_human_comp_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/dyad_human_human_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')

%% Compare conditions

nSample                                 = 30;
scnt                                    = 0;
snr                                     = solo_perf{1}.carr;

for iSubj = 1:size(dyad_pc_perf,2)
    if ~isempty(dyad_pc_perf{iSubj})
        id_dyad_pc{iSubj}           	= dyad_pc_perf{iSubj}.id;
    else
        id_dyad_pc{iSubj}           	= 'empty';
    end 
end

for iSubj = 1:size(dyad_perf,2)
    if ~isempty(dyad_perf{iSubj})
        id_dyad{iSubj}           	= dyad_perf{iSubj}.id;
    else
        id_dyad{iSubj}           	= 'empty';
    end 
end

for iSub = 1:length(solo_perf)

    if isempty(solo_perf{iSub})
        continue
    end
    
    idx_pc                          = cellfun(@(x) strcmp(x,solo_perf{iSub}.id),id_dyad_pc);
    
    if sum(idx_pc > 0)
        scnt                     	= scnt+1;
        hir_df(scnt,:)            	= dyad_pc_perf{idx_pc}.hir - solo_perf{iSub}.hir;
        
        for iCoh = 1:length(snr)
            
            auc_acc(scnt,iCoh)   	= getAUROC(solo_perf{iSub}.acc_trg{iCoh},dyad_pc_perf{idx_pc}.acc_trg{iCoh});
            auc_ecc(scnt,iCoh)     	= getAUROC(solo_perf{iSub}.ecc{iCoh},dyad_pc_perf{idx_pc}.ecc{iCoh});
            auc_score(scnt,iCoh)  	= getAUROC(solo_perf{iSub}.trg_score{iCoh},dyad_pc_perf{idx_pc}.trg_score{iCoh});
         
%             auc_acc(scnt,iCoh)          = f_auroc(solo_perf{iSub}.acc_trg{iCoh},dyad_pc_perf{idx_pc}.acc_trg{iCoh});
%             auc_ecc(scnt,iCoh)          = f_auroc(solo_perf{iSub}.ecc{iCoh},dyad_pc_perf{idx_pc}.ecc{iCoh});
%             auc_score(scnt,iCoh)        = f_auroc(solo_perf{iSub}.trg_score{iCoh},dyad_pc_perf{idx_pc}.trg_score{iCoh});
        end
    end
end

%% Load all agent sessions

sxc_agnt                    = [];
coh_agnt                    = [];
c                        	= 0;

for iComp = 1:length(agnt_perf)
    
    if isempty(agnt_perf{iComp})
        continue
    end
    
    c                       = c+1;
    macc_agnt(c,:)          = agnt_perf{iComp}.macc_trg;
    mecc_agnt(c,:)          = agnt_perf{iComp}.mecc;
    hir_agnt(c,:)           = agnt_perf{iComp}.hir;
    trg_score_agnt(c,:)     = agnt_perf{iComp}.trg_mscore;
    sxc_agnt                = [sxc_agnt; agnt_cr{iComp}.sxc];
    coh_agnt                = [coh_agnt; agnt_cr{iComp}.coh'];
end

%% PLOT

f                           = figure('units','normalized','position',[0 0 .5 1]);
height                    	= fliplr(linspace(.06,.85,4));
clmns                      	= linspace(.11,.75,3);
lb_fs                       = 14;
lg_fs                       = 10;
lw                          = 3;
frme_ms                     = 1000/120;
alp                         = .35;
avg_mult                    = 1.5;
nLag                        = 150;                                          % Cross-correlation lag       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Example agent joystick response %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [0.2 0.2]*1.5;
row                         = .65;
clm                         = .14; 
hofs                        = .1;
cmap                        = [0 0 0; gray(256)];
steps                       = .05;
bins                        = 0:steps:1;
c                           = 0;

ax0v                       	= axes('Position', [clm row-hofs dim(1) dim(2)/5]); hold on
ax0h                        = axes('Position', [clm-hofs row dim(1)/5 dim(2)]); hold on
ax0                       	= axes('Position', [clm row dim]); hold on

load('/Volumes/T7_Shield/CPR_psychophysics/RoH/summary/20220622_agnt_CPRagent_block1_tbl.mat')
t_agnt_plot                 = t;

% Extracte experimental data
clear trg_acc trg_conf trg_coh trg_hit
for iState = 1:size(t_agnt_plot,1)
    if t_agnt_plot.trg_shown(iState) == false
        continue
    end
    
    for iTarget = 1:length(t_agnt_plot.trg_ts{iState})
        c                   = c+1;
        trg_acc(c)          = t_agnt_plot.trg_acc{iState}(iTarget);
        trg_conf(c)         = t_agnt_plot.trg_ecc{iState}(iTarget);
        trg_coh(c)          = t_agnt_plot.rdp_coh(iState);
        trg_hit(c)          = t_agnt_plot.trg_hit{iState}(iTarget);
        
        if trg_acc(c) < .5 && trg_hit(c) == 1
            disp([num2str(iState) ' ' num2str(iTarget) ' ' num2str(c) ])
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
aidx                        = arc < 12.7587;
arc(aidx)                   = 12.7587;

% For each confidence level, calculate minimum accuracy required to hit
% the target at given arc width - normalised values
hit_width_acc               = 1 - ((arc/2) / 180);
hit_width_acc(aidx)         = 1 - (12.7587/2)/180; % arc width fixed

% Remove position that cannot yield reward from reward matrix
for iAcc = 1:length(acc)
    indx                    = conf < hit_width_acc(iAcc);
    rew(iAcc,indx)          = nan;
end

% Plot reward matrix
hold on
im                          = imagesc(acc,conf,rew);
ax0.XLabel.String           = 'Accuracy [norm]';
ax0.YLabel.String           = 'Eccentricity [norm]';
ax0.FontSize                = lb_fs;
ax0.XLim                    = [0 1];
ax0.YLim                    = [0 1];
cb                          = colorbar;
cb.Label.String             = '% Reward';
cb.Location                 = 'eastoutside';
ax0.XLim                    = [0 1];
ax0.XTick                   = [0:.2:1];
ax0.YTick                   = [0:.2:1];
ax0.XTickLabelRotation      = 0;
ax0.YLabel.Position(1)      = -.35;
ax0.XLabel.Position(2)      = -.22;

colormap(cmap)

cmap_coh                    = cool(size(snr,2));

for iCoh = 1:length(snr)
    cidx                    = trg_coh == snr(iCoh);
    sc                      = scatter(trg_acc(cidx), trg_conf(cidx), 'filled');
    sc.CData            	= cmap_coh(iCoh,:);
    sc.SizeData             = 20;
    sc.MarkerFaceAlpha      = .9;
end

ax0.Position                = [clm row dim];

nBin                        = 40;
axes(ax0h)
[h, edg]                    = histcounts(trg_conf,nBin);
cntr                        = edg(1:end-1) + diff(edg) ./ 2;
st                          = stairs(-h,cntr);
st.LineWidth                = lw/1.5;
st.Color                    = [0 0 0];
ax0h.YLim                   = [0 1];
ax0h.XAxis.Visible          = 'off';
ax0h.YAxis.Visible          = 'off';

axes(ax0v)
[v, edg]                    = histcounts(trg_acc,nBin);
cntr                        = edg(1:end-1) + diff(edg) ./ 2;
st                          = stairs(cntr,-v);
st.LineWidth                = lw/1.5;
st.Color                    = [0 0 0];
ax0v.XLim                   = [0 1];
ax0v.XAxis.Visible          = 'off';
ax0v.YAxis.Visible          = 'off';
uistack(ax0v,'bottom')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Performance AGNT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yof                         = .075;
xof                         = .075;
ax4                         = axes('Position', [clm-xof row-yof dim(1)/1.25 dim(2)/2]); hold on
cmap                        = cbrewer('div', 'PRGn', 4, 'PCHIP');

[pl1, ci1]                 	= plotAvgCI(macc_agnt, cmap(1,:), alp, lw);
[pl2, ci2]               	= plotAvgCI(hir_agnt, cmap(2,:), alp, lw);
[pl3, ci3]                	= plotAvgCI(mecc_agnt, cmap(3,:), alp, lw);
[pl4, ci4]                	= plotAvgCI(trg_score_agnt, cmap(4,:), alp, lw);

ax4.YLim                    = [40 100];
ax4.XLim                    = [1 size(hir_agnt,2)];
ax4.XLabel.String           = 'Coherence [%]';
ax4.YLabel.String           = 'Avg performance [%]';
ax4.XTick                   = 1:length(snr);
ax4.XTickLabel              = round(snr,2)*100;
ax4.FontSize                = lb_fs;
ax4.XTickLabelRotation      = 0;
ax4.Position                = [clmns(3)-xof row-yof dim(1)/1.25 dim(2)/2];
% ax4.XAxis.Visible           = 'off';

box off

lg                          = legend([pl1 pl2 pl3 pl4],{'Accuracy', 'Hit rate','Eccentricity' 'Score'}, 'NumColumns', 2)';
lg.Box                      = 'off';
lg.FontSize                 = lg_fs;
lg.Location                 = 'northeast';
lg.Position(1)              = .675;
lg.Position(2)              = .67;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Crosscorrelation AGNT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yof                         = .16;
ax3                         = axes('Position', [clm-xof row+yof dim(1)/1.25 dim(2)/2]); hold on
cmap                        = cool(size(snr,2));

for iCoh = 1:size(snr,2)
    avg_sxc{iCoh} = mean(sxc_agnt(coh_agnt == snr(iCoh),:));
end

clear pl
for iCoh = 1:size(avg_sxc,2)
    mmat(iCoh, :)           = avg_sxc{iCoh};
    pl(iCoh)               	= plot(mmat(iCoh, :),'LineWidth', lw, 'Color', [cmap(iCoh,:) .75]);
    lg_str{iCoh}            = num2str(round(snr(iCoh)*100));
end

avg_peak_pos                = find(mean(mmat) == max(mean(mmat)));
px                          = plot(avg_peak_pos, max(mean(mmat)), 'kx');
px.MarkerSize               = 8;
px.LineWidth                = lw/2;

ln                          = line([avg_peak_pos avg_peak_pos],[0 max(mean(mmat))]);
ln.LineStyle                = ':';
ln.LineWidth                = lw/2;
ln.Color                    = [.5 .5 .5];

% ln                          = line([150 150],[0 .2]);
% ln.LineStyle                = ':';
% ln.LineWidth                = lw;
% ln.Color                    = [0 0 0];

ax3.YLabel.String           = 'XC Coef';
ax3.XLabel.String           = 'Lag [ms]';
ax3.XLim                    = [150 301];
ax3.YLim                    = [0 .16];
ax3.XTick                   = [150 avg_peak_pos 300];
ax3.FontSize              	= lb_fs;
ax3.XTickLabel              = round((cellfun(@str2num, ax3.XTickLabel)-nLag) * frme_ms);
ax3.Position                = [clmns(3)-xof row+yof dim(1)/1.25 dim(2)/2];
ax3.XTickLabelRotation      = 0;

lg                          = legend(pl,lg_str,'Location','northwest','NumColumns', 2);
lg.Box                      = 'off';
lg.TextColor                = [.99 .99 .99];
lg.FontSize                 = lg_fs;
lg.Position(1)              = .17;
lg.Position(2)              = .67;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBPLOT: Lag difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax6                         = axes('Position', [clmns(1) height(3) dim]); hold on
% 
% for iL = 1:length(lag)
%     pl(iL)               	= plot([1 2],[lag(iL,1) lag(iL,2)]);
%     pl(iL).Color          	= [.5 .5 .5 alp];
%     pl(iL).LineWidth      	= lw;
% end
% 
% bx                          = boxplot(lag, 'Colors', 'k');
% set(bx,'MarkerEdgeColor','k')
% set(bx, {'linew'},{lw})
% 
% ax6.YLabel.String           = 'Time [ms]';
% ax6.XTick                   = [1 2];
% ax6.XLim                    = [.5 2.5];
% ax6.XTickLabel              = {'Solo','Agnt'};
% ax6.FontSize                = lb_fs;
% ax6.XTickLabelRotation      = 0;
% ax6.Position                = [clmns(1) height(3) dim];
% ax6.Box                     = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Score difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [.2 .2];

ax6                         = axes('Position', [clmns(1) height(3) dim]); hold on
ln                          = line([0 1], [0 1]);
ln.LineStyle                = '--';
ln.LineWidth                = lw/2;
ln.Color                    = [0 0 0];

for iSubj = 1:size(solo_perf,2)
    idx_pc              = cellfun(@(x) strcmp(x,solo_perf{iSubj}.id),id_dyad_pc);
    idx_dy            	= cellfun(@(x) strcmp(x,solo_perf{iSubj}.id),id_dyad);
    
    if sum(idx_pc>0) && sum(idx_dy>0)
        sc_agnt       	= scatter(mean(solo_perf{iSubj}.score_norm),mean(dyad_pc_perf{idx_pc}.score_norm), 'MarkerFaceColor', [.5 .5 .5],'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
        sc_dyad       	= scatter(mean(solo_perf{iSubj}.score_norm), mean(dyad_perf{idx_dy}.score_norm), 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
    end
end

ax6.YLabel.String           = {'Dyadic'};
ax6.XLabel.String           = {'Solo'};
ax6.FontSize                = lb_fs;
ax6.YLim                    = [0 .5];
ax6.XLim                    = [0 .5];

lg                          = legend([sc_agnt(1) sc_dyad(1)], 'COMP','HUMAN' );
lg.Location                 = 'southeast';
lg.FontSize                 = lg_fs;
lg.Box                      = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Hit rate comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSub = 1:size(solo_perf,2)
    hir.solo(iSub)   	= solo_perf{iSub}.hir_pool;
    idx_pc              = cellfun(@(x) strcmp(x,solo_perf{iSub}.id),id_dyad_pc);
    idx_dy            	= cellfun(@(x) strcmp(x,solo_perf{iSub}.id),id_dyad);

    if sum(idx_pc > 0)
        hir.dyad_pc(iSub)	= dyad_pc_perf{idx_pc}.hir_pool;
    else
        hir.dyad_pc(iSub) 	= nan;
    end
    if sum(idx_dy > 0)
        hir.dyad(iSub) 	= dyad_perf{idx_dy}.hir_pool;
    else
        hir.dyad(iSub) 	= nan;
  	end
end

nRep                        = 500;
e.solo                      = bootci(nRep, {@mean, hir.solo},'alpha', .001);
e.dyad_pc               	= bootci(nRep, {@nanmean, hir.dyad_pc},'alpha', .001);
e.dyad                      = bootci(nRep, {@nanmean, hir.dyad},'alpha', .001);

ax1                      	= axes('Position', [clmns(1) height(4) dim]); hold on
bp                          = bar([1:3],[mean(hir.solo) nanmean(hir.dyad_pc) nanmean(hir.dyad)]);
bp.FaceColor                = [.5 .5 .5];
bp.EdgeColor                = 'none';

er                          = errorbar([1:3],[mean(hir.solo) nanmean(hir.dyad_pc) nanmean(hir.dyad)],[nanmean(hir.solo)-e.solo(2) nanmean(hir.dyad_pc)-e.dyad_pc(2) nanmean(hir.dyad)-e.dyad(2)],[e.solo(1)-mean(hir.solo) e.dyad_pc(1)-nanmean(hir.dyad_pc) e.dyad(1)-nanmean(hir.dyad)]);    
er.Color                    = [0 0 0];                            
er.LineStyle                = 'none';
er.LineWidth                = lw/1.5;

ax1.YLabel.String           = 'Avg hit rate';
ax1.XLim                    = [.5 3.5];
ax1.YLim                    = [.3 .6];
ax1.YTick                 	= .3:.1:.6;
ax1.XTickLabel              = {'Solo','Comp','Dyad'};
ax1.FontSize                = lb_fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Hit rate difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax10                      	= axes('Position', [clmns(2) height(3) dim]); hold on

ln                          = line([0 8],[0 0], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax10                        = plotData(ax10, hir_df, 'Hit rate diff', snr, alp, lw, lb_fs, avg_mult,true);
ax10.YLim                 	= [-.1 .4];
ax10.YLim                 	= [-.2 .4];
ax10.YTick                  = [-.2 0 .2 .4];
ax10.Position             	= [clmns(2) height(3) dim];
ax10.XAxis.Visible          = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Score difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax11                    	= axes('Position', [clmns(3) height(3) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax11                        = plotData(ax11, auc_score, {'Score diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax11.Position           	= [clmns(3) height(3) dim];
ax11.XAxis.Visible          = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Eccentricity difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax12                        = axes('Position', [clmns(2) height(4) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax12                        = plotData(ax12, auc_ecc, {'Eccentricity diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax12.Position           	= [clmns(2) height(4) dim];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Accuracy difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax13                     	= axes('Position', [clmns(3) height(4) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax13                        = plotData(ax13, auc_acc, {'Accuracy diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax13.Position             	= [clmns(3) height(4) dim];

text(3.5,0.825, 'Comp > Solo', 'Parent', ax13, 'FontSize', lb_fs, 'Color', [0 0 0])
text(3.5,0.175, 'Comp < Solo', 'Parent', ax13, 'FontSize', lb_fs, 'Color', [0 0 0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax0                         = axes('Position',[0 0 1 1],'Visible','off');
lofs                        = .225;
% text(clmns(1)-.1,height(1)+lofs, 'A', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
% text(clmns(3)-.175,height(1)+lofs, 'B', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
% text(clmns(3)-.175,height(2)+lofs, 'C', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
% text(clmns(1)-.1,height(3)+lofs, 'D', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
% text(clmns(2)-.1,height(3)+lofs, 'E', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')

% text(0.2,.75, 'Coherence', 'Parent', ax0, 'FontSize', lg_fs, 'Color', [.99 .99 .99])

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG3/FIG3', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG3/FIG3', '-r500', '-dsvg', '-painters');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG3/FIG3', '-r500', '-depsc2', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pl, ci] = plotAvgCI(dat, col, alp, lw)

ci                	= bootci(500,{@median, dat.*100}, 'alpha', .01);
pt                  = patch([1:size(dat,2) fliplr(1:size(dat,2))], [ci(1,:) fliplr(ci(2,:))], col, 'FaceAlpha',alp, 'EdgeColor','none');

pl                 	= plot(median(dat).*100);
pl.Color           	= col;
pl.Marker       	= '.';
pl.MarkerSize       = 15;
pl.LineWidth       	= lw/2;
end

function [ax] = plotData(ax, dat, ylab, snr, alp, lw, lb_fs, avg_mult, flag)

if nargin < 9
    flag = false;
end

% if flag == false
%     pt                      = patch([0 8 8 0], [.5 .5 1 1], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
%     pt                      = patch([0 8 8 0], [0 0 .5 .5], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
% else
%     pt                      = patch([0 8 8 0], [0 0 .5 .5], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
%     pt                      = patch([0 8 8 0], [-.5 -.5 0 0], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
% end

for iL = 1:size(dat,1)
    pl(iL)               	= plot(dat(iL,:));
    pl(iL).Color          	= [.5 .5 .5 alp];
    pl(iL).LineWidth      	= lw/lw;
end

% Boostrap confidence intervals
nRep                        = 1000;
[CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.05);

% Prepare filled area
vec                         = 1:length(CI);
x_spacing                   = [vec fliplr(vec)];
ci                          = [CI(1,:) fliplr(CI(2,:))];

% Overlay confidence intervals
fl                          = fill(x_spacing,ci,[.3 0 0],'EdgeColor','none', 'FaceAlpha', alp);

% Plot mean curve
pm                          = plot(nanmean(dat),'LineWidth', lw/1.5, 'Color', [0 0 0]);

ax.XLim                     = [1 size(dat,2)];
ax.XLabel.String            = 'Coherence [%]';
ax.YLabel.String            = ylab;
ax.XTick                    = [1 4 length(snr)];
ax.XTickLabel               = round([snr(1) snr(4) snr(end)],2)*100;
ax.FontSize                 = lb_fs;
ax.YLim                     = [.1 .9];
ax.YTick                    = [.25 .5 .75];
ax.XTickLabelRotation       = 0;

box off
end

function [out] = getAUROC(in1, in2)

if size(in1,1) == 1
    in1         = in1';
    in2         = in2';
end

lab          	= [zeros(length(in1),1); ones(length(in2),1)];
[~,~,~,out]     = perfcurve(lab,[in1; in2],1);

end
