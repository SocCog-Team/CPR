% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

source_pth = '/Users/fschneider/ownCloud/var_plot/';
load([source_pth '/solo_performance.mat'])
load([source_pth '/hh_dyad_performance.mat'])
load([source_pth '/hc_dyad_performance.mat'])
load([source_pth '/hc_dyad_correlation.mat'])
load([source_pth '/comp_performance.mat'])
load([source_pth '/comp_correlation.mat'])

%% Compare conditions

nSample                                 = 30;
scnt                                    = 0;
snr                                     = solo_perf{1}.carr;

% Extract subject sequence
for iSubj = 1:size(hc_dyad_perf,2)
    if ~isempty(hc_dyad_perf{iSubj})
        id_hc_dyad{iSubj}           	= hc_dyad_perf{iSubj}.id;
    else
        id_hc_dyad{iSubj}           	= 'empty';
    end 
end

for iSubj = 1:size(dyad_perf,2)
    if ~isempty(dyad_perf{iSubj})
        id_dyad{iSubj}           	= dyad_perf{iSubj}.id;
    else
        id_dyad{iSubj}           	= 'empty';
    end 
end

% Contrast data
scnt = 0;
for iSub = 1:length(solo_perf)

    if isempty(solo_perf{iSub})
        continue
    end
    
    idx_pc                          = cellfun(@(x) strcmp(x,solo_perf{iSub}.id),id_hc_dyad);
    idx_dy                          = cellfun(@(x) strcmp(x,solo_perf{iSub}.id),id_dyad);
    
    if sum(idx_pc) == 0 || sum(idx_dy) == 0
        continue
    end
    
    scnt = scnt+1;
    hir_df_SC(scnt,:)               = hc_dyad_perf{idx_pc}.hir - solo_perf{iSub}.hir;
    hir_df_SH(scnt,:)               = dyad_perf{idx_pc}.hir - solo_perf{iSub}.hir;
    hir_df_CH(scnt,:)               = hc_dyad_perf{idx_pc}.hir - dyad_perf{idx_dy}.hir;

    auc_acc_pooled_SC(scnt)         = getAUROC(cell2mat(solo_perf{iSub}.acc_trg),cell2mat(hc_dyad_perf{idx_pc}.acc_trg));
    auc_acc_pooled_SH(scnt)         = getAUROC(cell2mat(solo_perf{iSub}.acc_trg),cell2mat(dyad_perf{idx_dy}.acc_trg));
    auc_acc_pooled_CH(scnt)         = getAUROC(cell2mat(hc_dyad_perf{idx_pc}.acc_trg),cell2mat(dyad_perf{idx_dy}.acc_trg));
    
    auc_ecc_pooled_SC(scnt)         = getAUROC(cell2mat(solo_perf{iSub}.ecc_state'),cell2mat(hc_dyad_perf{idx_pc}.ecc_state'));
    auc_ecc_pooled_SH(scnt)         = getAUROC(cell2mat(solo_perf{iSub}.ecc_state'),cell2mat(dyad_perf{idx_dy}.ecc_state'));
    auc_ecc_pooled_CH(scnt)         = getAUROC(cell2mat(hc_dyad_perf{idx_pc}.ecc_state'),cell2mat(dyad_perf{idx_dy}.ecc_state'));

    for iCoh = 1:length(snr)
        % Solo - Computer dyad
%         auc_acc_SC(scnt,iCoh)   	= getAUROC(solo_perf{iSub}.acc_trg{iCoh},hc_dyad_perf{idx_pc}.acc_trg{iCoh});
%         auc_ecc_SC(scnt,iCoh)     	= getAUROC(solo_perf{iSub}.ecc_state{iCoh},hc_dyad_perf{idx_pc}.ecc_state{iCoh});
%         auc_score_SC(scnt,iCoh)  	= getAUROC(solo_perf{iSub}.trg_score{iCoh},hc_dyad_perf{idx_pc}.trg_score{iCoh});
            
        auc_acc_SC(scnt,iCoh)   	= getAUROC(hc_dyad_perf{idx_pc}.acc_trg{iCoh},solo_perf{iSub}.acc_trg{iCoh});
        auc_ecc_SC(scnt,iCoh)     	= getAUROC(hc_dyad_perf{idx_pc}.ecc_state{iCoh},solo_perf{iSub}.ecc_state{iCoh});
        auc_score_SC(scnt,iCoh)  	= getAUROC(hc_dyad_perf{idx_pc}.trg_score{iCoh},solo_perf{iSub}.trg_score{iCoh});
        
        % Solo - Human dyad
        auc_acc_SH(scnt,iCoh)   	= getAUROC(solo_perf{iSub}.acc_trg{iCoh},dyad_perf{idx_dy}.acc_trg{iCoh});
        auc_ecc_SH(scnt,iCoh)     	= getAUROC(solo_perf{iSub}.ecc_state{iCoh},dyad_perf{idx_dy}.ecc_state{iCoh});
        auc_score_SH(scnt,iCoh)  	= getAUROC(solo_perf{iSub}.trg_score{iCoh},dyad_perf{idx_dy}.trg_score{iCoh});
        
        % Human dyad - Computer dyad
%         auc_acc_CH(scnt,iCoh)       = getAUROC(hc_dyad_perf{idx_pc}.acc_trg{iCoh},dyad_perf{idx_dy}.acc_trg{iCoh});
%         auc_ecc_CH(scnt,iCoh)     	= getAUROC(hc_dyad_perf{idx_pc}.ecc_state{iCoh},dyad_perf{idx_dy}.ecc_state{iCoh});
%         auc_score_CH(scnt,iCoh)  	= getAUROC(hc_dyad_perf{idx_pc}.trg_score{iCoh},dyad_perf{idx_dy}.trg_score{iCoh});
                
        auc_acc_CH(scnt,iCoh)       = getAUROC(dyad_perf{idx_dy}.acc_trg{iCoh},hc_dyad_perf{idx_pc}.acc_trg{iCoh});
        auc_ecc_CH(scnt,iCoh)     	= getAUROC(dyad_perf{idx_dy}.ecc_state{iCoh},hc_dyad_perf{idx_pc}.ecc_state{iCoh});
        auc_score_CH(scnt,iCoh)  	= getAUROC(dyad_perf{idx_dy}.trg_score{iCoh},hc_dyad_perf{idx_pc}.trg_score{iCoh});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE SETTINGS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                           = figure('units','centimeters','position',[0 0 16 16]);
dim                         = [.2 .2];
height                      = linspace(.79, .075,4);
clmns                       = linspace(.1,.78,3);
lb_fs                       = 8;
lg_fs                       = 8;
lw                          = 3;
frme_ms                     = 1000/120;
alp                         = .4;
coh_col                     = cool(length(snr));
col_dat                     = [0 0 0];
col_ci                      = [.3 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Score/hit rate difference Dyad - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax0                         = axes('Position', [clmns(1) height(1) dim]); hold on
ln                          = line([0 1], [0 1]);
ln.LineStyle                = ':';
ln.LineWidth                = lw/2;
ln.Color                    = [0 0 0];

for iSubj = 1:size(solo_perf,2)
    idx_pc                  = cellfun(@(x) strcmp(x,solo_perf{iSubj}.id),id_hc_dyad);
    idx_dy                  = cellfun(@(x) strcmp(x,solo_perf{iSubj}.id),id_dyad);
    
    if sum(idx_pc>0) && sum(idx_dy>0)
        sc_agnt             = scatter(mean(solo_perf{iSubj}.trg_all_score),mean(hc_dyad_perf{idx_pc}.trg_all_score), 'MarkerFaceColor', [.6 .6 .6],'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
        sc_dyad             = scatter(mean(solo_perf{iSubj}.trg_all_score), mean(dyad_perf{idx_dy}.trg_all_score), 'MarkerFaceColor', [.3 .3 .3],'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);

    end
end

ax0.YLabel.String           = {'Dyadic'};
ax0.XLabel.String           = {'Solo'};
ax0.FontSize                = lb_fs;
ax0.YLim                    = [0 .5];
ax0.XLim                    = [0 .5];

lg                          = legend([sc_agnt(1) sc_dyad(1)], 'HC dyad','HH dyad');
lg.Location                 = 'southeast';
lg.Interpreter              = 'none';
lg.FontSize                 = lg_fs;
lg.Box                      = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Average hit rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSub = 1:size(solo_perf,2)
    hir.solo(iSub)   	= solo_perf{iSub}.hir_pool;
    scr.solo(iSub)   	= mean(solo_perf{iSub}.trg_all_score);
    idx_pc              = cellfun(@(x) strcmp(x,solo_perf{iSub}.id),id_hc_dyad);
    idx_dy            	= cellfun(@(x) strcmp(x,solo_perf{iSub}.id),id_dyad);

    if sum(idx_pc > 0)
        hir.dyad_pc(iSub)	= hc_dyad_perf{idx_pc}.hir_pool;
        scr.dyad_pc(iSub)   = mean(hc_dyad_perf{idx_pc}.trg_all_score);
    else
        hir.dyad_pc(iSub) 	= nan;
        scr.dyad_pc(iSub) 	= nan;
    end
    if sum(idx_dy > 0)
        hir.dyad(iSub)      = dyad_perf{idx_dy}.hir_pool;
        scr.dyad(iSub)      = mean(dyad_perf{idx_dy}.trg_all_score);
    else
        hir.dyad(iSub)      = nan;
        scr.dyad(iSub)      = nan;
  	end
end

ax11                      	= axes('Position', [clmns(2) height(1) dim]); hold on
ax11                     	= plotBar(ax11,scr,'Mean score',lb_fs,lw);
ax11.YLim               	= [.2 .35];
ax11.YTick                 	= .2:.05:.35;

ax12                       	= axes('Position', [clmns(3) height(1) dim]); hold on
ax12                     	= plotBar(ax12,hir,'Mean hit rate',lb_fs,lw);
ax12.YLim                	= [.3 .6];
ax12.YTick                 	= .3:.1:.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Hit rate difference Computer-Human %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax2                         = axes('Position', [clmns(2) height(2) dim]); hold on
[ax2,pl]                    = plotData(ax2,hir_df_CH,true,lw,alp,col_dat,col_ci);
lm                          = line([1 7],[0 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth',lw/2);
ax2.XLim                    = [1 size(hir_df_CH,2)];
ax2.FontSize                = lb_fs;
ax2.YLabel.String           = 'Difference';
ax2.XAxis.Visible           = 'off';
ax2.YLim                    = [-15 40];
tx                        	= text(2.25,-10, 'HH dyad > HC dyad', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(2.25,35, 'HC dyad > HH dyad', 'FontSize', lb_fs, 'Color', 'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: AUROC Accuracy Computer-Human %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3                         = axes('Position', [clmns(2) height(3) dim]); hold on
ax3                         = plotAUROC(ax3,auc_acc_CH,'AUC',lb_fs,snr,alp,lw,col_dat,col_ci);
tx                        	= text(2.25,.1, 'HH dyad > HC dyad', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(2.25,.9, 'HC dyad > HH dyad', 'FontSize', lb_fs, 'Color', 'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: AUROC eccentricity Computer-Human %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax4                         = axes('Position', [clmns(2) height(4) dim]); hold on
ax4                         = plotAUROC(ax4,auc_ecc_CH,'AUC',lb_fs,snr,alp,lw,col_dat,col_ci);
ax4.XAxis.Visible           = 'on';
tx                        	= text(2.25,.1, 'HH dyad > HC dyad', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(2.25,.9, 'HC dyad > HH dyad', 'FontSize', lb_fs, 'Color', 'k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Scatter accuracy average Human-Computer %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax7                     	= axes('Position', [clmns(1) height(2) dim]); hold on
ax7                        	= plotScatter(ax7,solo_perf,hc_dyad_perf, dyad_perf,coh_col,'hir',lb_fs,id_hc_dyad,id_dyad,snr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Scatter accuracy average Human-Computer %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax5                        	= axes('Position', [clmns(1) height(3) dim]); hold on
ax5                       	= plotScatter(ax5,solo_perf,hc_dyad_perf, dyad_perf,coh_col,'macc_trg',lb_fs,id_hc_dyad,id_dyad,snr);
legend off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Scatter eccentricity average Human-Computer %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax6                        	= axes('Position', [clmns(1) height(4) dim]); hold on
ax6                      	= plotScatter(ax6,solo_perf,hc_dyad_perf, dyad_perf,coh_col,'mecc_state',lb_fs,id_hc_dyad,id_dyad,snr);
ax6.XAxis.Visible       	= 'on';
legend off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Scatter subjetc-wise accuracy AUROC  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax8                        	= axes('Position', [clmns(3) height(3) dim]); hold on
% ax8                         = scatterAUROC(ax8,auc_acc_CH, auc_acc_SC,lb_fs,coh_col);

ax10                      	= axes('Position', [clmns(3) height(3) dim(1) dim(2)/2.1]); hold on
ax10                        = histogramAUROC(ax10,auc_acc_CH,lb_fs,coh_col,lw);
ax11                      	= axes('Position', [clmns(3) height(3)+dim(2)/2 dim(1) dim(2)/2.1]); hold on
ax11                        = histogramAUROC(ax11,1-auc_acc_SC,lb_fs,coh_col,lw); % mirror data to bring CH data on same side of plot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Scatter subjetc-wise eccentricity AUROC  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax9                        	= axes('Position', [clmns(3) height(4) dim]); hold on
% ax9                         = scatterAUROC(ax9,auc_ecc_CH, auc_ecc_SC,lb_fs,coh_col);
% ax9.XAxis.Visible         	= 'on';

ax10                      	= axes('Position', [clmns(3) height(4) dim(1) dim(2)/2.1]); hold on
ax10                        = histogramAUROC(ax10,auc_ecc_CH,lb_fs,coh_col,lw);
ax10.XAxis.Visible          = 'on';
tx                        	= text(.6,12, 'HC > HH', 'FontSize', lb_fs, 'Color', 'k');
tx                        	= text(.1,12, 'HH > HC', 'FontSize', lb_fs, 'Color', 'k');

ax11                      	= axes('Position', [clmns(3) height(4)+dim(2)/2 dim(1) dim(2)/2.1]); hold on
ax11                        = histogramAUROC(ax11,1-auc_ecc_SC,lb_fs,coh_col,lw); % mirror data to bring CH data on same side of plot
tx                        	= text(.6,12, 'HC > S', 'FontSize', lb_fs, 'Color', 'k');
tx                        	= text(.1,12, 'S > HC', 'FontSize', lb_fs, 'Color', 'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Scatter subjetc-wise hit rate difference  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax10                      	= axes('Position', [clmns(3) height(2) dim]); hold on
% ax10                    	= scatterAUROC(ax10,hir_df_CH, hir_df_SC,lb_fs,coh_col);

ax10                      	= axes('Position', [clmns(3) height(2) dim(1) dim(2)/2.1]); hold on
ax10                        = histogramAUROC(ax10,hir_df_CH.*100,lb_fs,coh_col,lw,false);
ax10.XAxis.Visible          = 'on';

ax11                      	= axes('Position', [clmns(3) height(2)+dim(2)/2 dim(1) dim(2)/2.1]); hold on
ax11                        = histogramAUROC(ax11,hir_df_SC.*100,lb_fs,coh_col,lw,false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Export
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG3/FIG3', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG3/FIG3', '-r500', '-dsvg', '-painters');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG3/FIG3', '-r500', '-depsc2', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reported stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accuracy effect
better_all_coh  = sum(sum(auc_acc_CH < .5,2) == size(auc_acc_CH,2)) / size(auc_acc_CH,1);
better_zero_coh = sum(auc_acc_CH(:,1) < .5) / size(auc_acc_CH,1);

% Eccentricity effect
worse_all_coh  = sum(sum(auc_ecc_CH > .5,2) == size(auc_ecc_CH,2)) / size(auc_ecc_CH,1);
worse_98_coh = sum(auc_ecc_CH(:,end) > .5) / size(auc_ecc_CH,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ax = plotBar(ax,in,lab_str,lb_fs,lw)
axes(ax); hold on

nRep                        = 500;
e.solo                      = bootci(nRep, {@mean, in.solo},'alpha', .01);
e.dyad_pc               	= bootci(nRep, {@nanmean, in.dyad_pc},'alpha', .01);
e.dyad                      = bootci(nRep, {@nanmean, in.dyad},'alpha', .01);

bp                          = bar([1:3],[mean(in.solo) nanmean(in.dyad_pc) nanmean(in.dyad)]);
bp.FaceColor                = [.5 .5 .5];
bp.EdgeColor                = 'none';

er                          = errorbar([1:3],[mean(in.solo) nanmean(in.dyad_pc) nanmean(in.dyad)],[nanmean(in.solo)-e.solo(2) nanmean(in.dyad_pc)-e.dyad_pc(2) nanmean(in.dyad)-e.dyad(2)],[e.solo(1)-mean(in.solo) e.dyad_pc(1)-nanmean(in.dyad_pc) e.dyad(1)-nanmean(in.dyad)]);
er.Color                    = [0 0 0];
er.LineStyle                = 'none';
er.LineWidth                = lw/1.5;

ax.YLabel.String            = lab_str;
ax.XLim                     = [.5 3.5];
ax.XTick                    = 1:3;
ax.XTickLabel               = {'Solo','HC dyad','HH dyad'};
ax.XTickLabelRotation      	= 15;
ax.FontSize                 = lb_fs;
end

function [ax,pl] = plotData(ax,dat,scale_flag,lw,alp,col_dat,col_ci)

% Remove NaNs and zero rows
dat(sum(isnan(dat),2)>0,:)  = [];
dat(sum(dat==0,2)>0,:)      = [];

if scale_flag
    dat                     = dat.*100;
end

% Boostrap confidence intervals
nRep                        = 1000;
[CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.01);

% Prepare filled area
vec                         = 1:length(CI);
x_spacing                   = [vec fliplr(vec)];
ci                          = [CI(1,:) fliplr(CI(2,:))];

% Plot subject-wise data
hold on
for iL = 1:size(dat,1)
    pl(iL)                  = plot(dat(iL,:), 'LineWidth', lw/lw, 'Color', [.5 .5 .5 alp]);
end

% Overlay confidence intervals
fl                          = fill(x_spacing,ci,col_ci,'EdgeColor','none', 'FaceAlpha', alp);

% Plot mean curve
pl                          = plot(mean(dat), 'LineWidth', lw/1.5, 'Color', col_dat);
end

function ax = plotAUROC(ax,dat,str,lb_fs,snr,alp,lw,col_dat,col_ci)

axes(ax); hold on

% Remove zero rows
dat(sum(dat == 0,2) == size(dat,2),:) = [];

% Plot subject-wise data
for iL = 1:size(dat,1)
    pl(iL)               	= plot(dat(iL,:));
    pl(iL).Color          	= [.5 .5 .5 alp];
    pl(iL).LineWidth      	= lw/lw;
end

% Boostrap confidence intervals
nRep                        = 1000;
[CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.01);

% Prepare filled area
vec                         = 1:length(CI);
x_spacing                   = [vec fliplr(vec)];
ci                          = [CI(1,:) fliplr(CI(2,:))];

% Overlay confidence intervals
fl                          = fill(x_spacing,ci,col_ci,'EdgeColor','none', 'FaceAlpha', alp);

% Plot mean curve
pm                          = plot(mean(dat),'LineWidth', lw/1.5, 'Color', col_dat);

ax.XLim                  	= [1 length(snr)];
ax.YLim                  	= [0 1];
ax.FontSize               	= lb_fs;
ax.YLabel.String           	= str;
ax.XLabel.String           	= 'Coherence [%]';
ax.XTick                 	= 1:7;
ax.YTick                  	= [.25 .5 .75];
ax.XTickLabel            	= round(snr,2).*100;
ax.XTickLabelRotation      	= 0;
ax.XAxis.Visible            = 'off';

% n                           = .05;
% crit                        = [.5-n .5+n];
lm                          = line([1 7],[.5 .5], 'Color', 'k', 'LineStyle', ':', 'LineWidth',lw/2);
% ll                          = line([1 7],[crit(1) crit(1)], 'Color', 'k', 'LineStyle', '-.', 'LineWidth',lw/2);
% lh                          = line([1 7],[crit(2) crit(2)], 'Color', 'k', 'LineStyle', '-.', 'LineWidth',lw/2);

end

function ax = plotScatter(ax,in_solo,in_pc, in_hum,coh_col,cond_str,lb_fs,id_dyad_pc,id_dyad,snr)

axes(ax); hold on

for iSubj = 1:size(in_solo,2)
    idx_pc                  = cellfun(@(x) strcmp(x,in_solo{iSubj}.id),id_dyad_pc);
    idx_dy                  = cellfun(@(x) strcmp(x,in_solo{iSubj}.id),id_dyad);
    
    if sum(idx_pc>0) && sum(idx_dy>0)
        for iCoh = 1:length(snr)
            
            y_mat(iSubj,iCoh)   = in_pc{idx_pc}.(cond_str)(iCoh);
            x_mat(iSubj,iCoh)   = in_hum{idx_dy}.(cond_str)(iCoh);
            sc(iCoh)            = scatter(in_hum{idx_dy}.(cond_str)(iCoh),in_pc{idx_pc}.(cond_str)(iCoh), 'MarkerFaceColor', coh_col(iCoh,:)./2,'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .3);
            lab{iCoh}           = num2str(round(snr(iCoh),2));
        end
    else
        x_mat(iSubj,:)      = nan([1 length(snr)]);
        y_mat(iSubj,:)      = nan([1 length(snr)]);
    end
end

for iCoh = 1:length(snr)
    xx              = x_mat(:,iCoh);
    yy              = y_mat(:,iCoh);
    x_ci            = (bootci(500, {@nanmedian,  xx},'alpha', .001));
    y_ci            = (bootci(500, {@nanmedian,  yy},'alpha', .001));
    lny             = line([nanmedian(xx) nanmedian(xx)],[y_ci(1) y_ci(2)], 'Color', coh_col(iCoh,:),'LineWidth',1);
    lnx             = line([x_ci(1) x_ci(2)],[nanmedian(yy) nanmedian(yy)], 'Color', coh_col(iCoh,:),'LineWidth',1);
    sc              = scatter(nanmedian(xx),nanmedian(yy), 'MarkerFaceColor', coh_col(iCoh,:),'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1,'SizeData', 20);
end

ax.FontSize                 = lb_fs;
ax.YLabel.String            = 'HC dyad';
ax.XLabel.String            = 'HH dyad';
ax.XLim                 	= [0 1];
ax.YLim                     = [0 1];
ax.XTick                    = [.1 .5 .9];
ax.YTick                    = [.1 .5 .9];
ax.XAxis.Visible         	= 'off';

ln                          = line([0 1],[0 1]);
ln.LineWidth                = 1.5;
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];

end

function ax = scatterAUROC(ax,inCH, inSC,lb_fs,coh_col)

axes(ax); hold on

xdat                        = abs(inCH - .5);
ydat                        = abs(inSC - .5);

for iCoh = 1:7
    sc(iCoh)                = scatter(xdat(:,iCoh),ydat(:,iCoh), 'MarkerFaceColor', coh_col(iCoh,:),'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
end

ax.FontSize                 = lb_fs;
ax.XLabel.String            = {'Effect size'; '[Computer vs Human]'};
ax.YLabel.String            = {'Effect size'; '[Solo vs Computer]'};
ax.XLim                 	= [0 .6];
ax.YLim                     = [0 .6];
ax.XTick                    = [0 .3 .6];
ax.YTick                    = [0 .3 .6];
ax.XAxis.Visible         	= 'off';

ln                          = line([0 1],[0 1]);
ln.LineWidth                = 1.5;
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];

end

function ax = histogramAUROC(ax,in,lb_fs,coh_col,lw,flag)

axes(ax); hold on

if nargin < 6
    flag = true;
end

if flag
    edges = 0:.05:1;
else
    edges = 0:2.5:50;
end

for iCoh = 1:size(in,2)
    h = histogram(in(:,iCoh));
    h.BinEdges = edges;
    h.FaceColor = coh_col(iCoh,:);
    h.EdgeColor = 'none';
    h.FaceAlpha = .5;
    h.EdgeAlpha = .5;
end

ax.FontSize                 = lb_fs;
ax.YLabel.String            = '#';
ax.XAxis.Visible            = 'off';

if flag
    lm                  	= line([.5 .5],[0 15], 'Color', 'k', 'LineStyle', ':', 'LineWidth',lw/2);
    ax.XLabel.String     	= 'AUC';
    ax.XTick               	= [0 .5 1];
    ax.YTick               	= [0 10];
    ax.YLim               	= [0 15];
else
    ax.YLim             	= [0 12];
    ax.XLim               	= [0 50];
    ax.XTick              	= [0 25 50];
    ax.XLabel.String     	= 'Difference';
end

end

function [out] = getAUROC(in1, in2)

if size(in1,1) == 1
    in1         = in1';
    in2         = in2';
end

lab          	= [zeros(length(in1),1); ones(length(in2),1)];
[~,~,~,out]     = perfcurve(lab,[in1; in2],1);

end