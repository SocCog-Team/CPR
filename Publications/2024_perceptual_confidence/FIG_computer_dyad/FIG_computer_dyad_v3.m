% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

% Adjust path
source_pth = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/var_plot/';
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

%     auc_acc_pooled_SC(scnt)         = getAUROC(cell2mat(solo_perf{iSub}.acc_trg),cell2mat(hc_dyad_perf{idx_pc}.acc_trg));
%     auc_acc_pooled_SH(scnt)         = getAUROC(cell2mat(solo_perf{iSub}.acc_trg),cell2mat(dyad_perf{idx_dy}.acc_trg));
%     auc_acc_pooled_CH(scnt)         = getAUROC(cell2mat(dyad_perf{idx_dy}.acc_trg),cell2mat(hc_dyad_perf{idx_pc}.acc_trg)); 
    
    auc_acc_pooled_SC(scnt)         = getAUROC(cell2mat(solo_perf{iSub}.acc_state'),cell2mat(hc_dyad_perf{idx_pc}.acc_state'));
    auc_acc_pooled_SH(scnt)         = getAUROC(cell2mat(solo_perf{iSub}.acc_state'),cell2mat(dyad_perf{idx_dy}.acc_state'));
    auc_acc_pooled_CH(scnt)         = getAUROC(cell2mat(dyad_perf{idx_dy}.acc_state'),cell2mat(hc_dyad_perf{idx_pc}.acc_state'));
    
    auc_ecc_pooled_SC(scnt)         = getAUROC(cell2mat(solo_perf{iSub}.ecc_state'),cell2mat(hc_dyad_perf{idx_pc}.ecc_state'));
    auc_ecc_pooled_SH(scnt)         = getAUROC(cell2mat(solo_perf{iSub}.ecc_state'),cell2mat(dyad_perf{idx_dy}.ecc_state'));
    auc_ecc_pooled_CH(scnt)         = getAUROC(cell2mat(dyad_perf{idx_dy}.ecc_state'),cell2mat(hc_dyad_perf{idx_pc}.ecc_state'));

    p_ecc_pooled_CH(scnt)        	= ranksum(cell2mat(dyad_perf{idx_dy}.ecc_state'),cell2mat(hc_dyad_perf{idx_pc}.ecc_state'));
    p_acc_pooled_CH(scnt)       	= ranksum(cell2mat(dyad_perf{idx_dy}.acc_state'),cell2mat(hc_dyad_perf{idx_pc}.acc_state'));

    for iCoh = 1:length(snr)
        % Solo - Computer dyad
%         auc_acc_SC(scnt,iCoh)   	= getAUROC(solo_perf{iSub}.acc_trg{iCoh},hc_dyad_perf{idx_pc}.acc_trg{iCoh});
%         auc_ecc_SC(scnt,iCoh)     	= getAUROC(solo_perf{iSub}.ecc_state{iCoh},hc_dyad_perf{idx_pc}.ecc_state{iCoh});
%         auc_score_SC(scnt,iCoh)  	= getAUROC(solo_perf{iSub}.trg_score{iCoh},hc_dyad_perf{idx_pc}.trg_score{iCoh});
            
%         auc_acc_SC(scnt,iCoh)   	= getAUROC(hc_dyad_perf{idx_pc}.acc_trg{iCoh},solo_perf{iSub}.acc_trg{iCoh});
        auc_acc_SC(scnt,iCoh)   	= getAUROC(hc_dyad_perf{idx_pc}.acc_state{iCoh},solo_perf{iSub}.acc_state{iCoh});
        auc_ecc_SC(scnt,iCoh)     	= getAUROC(hc_dyad_perf{idx_pc}.ecc_state{iCoh},solo_perf{iSub}.ecc_state{iCoh});
        auc_score_SC(scnt,iCoh)  	= getAUROC(hc_dyad_perf{idx_pc}.trg_score{iCoh},solo_perf{iSub}.trg_score{iCoh});
        
        % Solo - Human dyad
%         auc_acc_SH(scnt,iCoh)   	= getAUROC(solo_perf{iSub}.acc_trg{iCoh},dyad_perf{idx_dy}.acc_trg{iCoh});
        auc_acc_SH(scnt,iCoh)   	= getAUROC(solo_perf{iSub}.acc_state{iCoh},dyad_perf{idx_dy}.acc_state{iCoh});
        auc_ecc_SH(scnt,iCoh)     	= getAUROC(solo_perf{iSub}.ecc_state{iCoh},dyad_perf{idx_dy}.ecc_state{iCoh});
        auc_score_SH(scnt,iCoh)  	= getAUROC(solo_perf{iSub}.trg_score{iCoh},dyad_perf{idx_dy}.trg_score{iCoh});
        
        % Human dyad - Computer dyad
%         auc_acc_CH(scnt,iCoh)       = getAUROC(hc_dyad_perf{idx_pc}.acc_trg{iCoh},dyad_perf{idx_dy}.acc_trg{iCoh});
%         auc_ecc_CH(scnt,iCoh)     	= getAUROC(hc_dyad_perf{idx_pc}.ecc_state{iCoh},dyad_perf{idx_dy}.ecc_state{iCoh});
%         auc_score_CH(scnt,iCoh)  	= getAUROC(hc_dyad_perf{idx_pc}.trg_score{iCoh},dyad_perf{idx_dy}.trg_score{iCoh});
                
%         auc_acc_CH(scnt,iCoh)       = getAUROC(dyad_perf{idx_dy}.acc_trg{iCoh},hc_dyad_perf{idx_pc}.acc_trg{iCoh});
        auc_acc_CH(scnt,iCoh)       = getAUROC(dyad_perf{idx_dy}.acc_state{iCoh},hc_dyad_perf{idx_pc}.acc_state{iCoh});
        auc_ecc_CH(scnt,iCoh)     	= getAUROC(dyad_perf{idx_dy}.ecc_state{iCoh},hc_dyad_perf{idx_pc}.ecc_state{iCoh});
        auc_score_CH(scnt,iCoh)  	= getAUROC(dyad_perf{idx_dy}.trg_score{iCoh},hc_dyad_perf{idx_pc}.trg_score{iCoh});
        
%         p_acc_CH(scnt,iCoh)            = ranksum(dyad_perf{idx_dy}.acc_trg{iCoh},hc_dyad_perf{idx_pc}.acc_trg{iCoh});
        p_acc_CH(scnt,iCoh)            = ranksum(dyad_perf{idx_dy}.acc_state{iCoh},hc_dyad_perf{idx_pc}.acc_state{iCoh});
        p_ecc_CH(scnt,iCoh)            = ranksum(dyad_perf{idx_dy}.ecc_state{iCoh},hc_dyad_perf{idx_pc}.ecc_state{iCoh});
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
% [ax2,pl]                    = plotData(ax2,hir_df_CH,true,lw,alp,col_dat,col_ci);
ax2                         = add_mean_to_bar(ax2,hir_df_CH.*100,[-10 30],[-10:10:30],lw,alp,col_ci);

lm                          = line([0 8],[0 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth',lw/2);
ax2.XLim                    = [0 size(hir_df_CH,2)+1];
ax2.FontSize                = lb_fs;
ax2.YLabel.String           = 'Difference';
ax2.XAxis.Visible           = 'off';
ax2.YAxis(1).Color           = 'none';
tx                        	= text(2,-5, 'HH dyad > HC dyad', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(2,25, 'HC dyad > HH dyad', 'FontSize', lb_fs, 'Color', 'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBPLOT: AUROC Accuracy Computer-Human %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax3                         = axes('Position', [clmns(2) height(3) dim]); hold on
% ax3                         = plotAUROC(ax3,auc_acc_CH,'AUC',lb_fs,snr,alp,lw,col_dat,col_ci);
% tx                        	= text(2.25,.1, 'HH dyad > HC dyad', 'FontSize', lb_fs, 'Color', 'k');
% tx                       	= text(2.25,.9, 'HC dyad > HH dyad', 'FontSize', lb_fs, 'Color', 'k');

sig_boundary                = .05 / (size(p_acc_CH,1) * size(p_acc_CH,2)); % Bonferroni correction

ax3                         = axes('Position', [clmns(2) height(3) dim]); hold on
ax3                         = plotBar_coherence(ax3, auc_acc_CH, p_acc_CH, sig_boundary, lb_fs, snr);
ax3.YLim                    = [-40 80];
ax3                         = add_mean_to_bar(ax3,auc_acc_CH,[.35 .8],[.4:.1:8],lw,alp,col_ci);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: AUROC eccentricity Computer-Human %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax4                         = plotAUROC(ax4,auc_ecc_CH,'AUC',lb_fs,snr,alp,lw,col_dat,col_ci);
% tx                        	= text(2.25,.1, 'HH dyad > HC dyad', 'FontSize', lb_fs, 'Color', 'k');
% tx                       	= text(2.25,.9, 'HC dyad > HH dyad', 'FontSize', lb_fs, 'Color', 'k');

sig_boundary                = .05 / (size(p_ecc_CH,1) * size(p_ecc_CH,2)); % Bonferroni correction

ax4                         = axes('Position', [clmns(2) height(4) dim]); hold on
ax4                         = plotBar_coherence(ax4, auc_ecc_CH, p_ecc_CH, sig_boundary, lb_fs, snr);
ax4.YLim                    = [-80 40];
ax4                         = add_mean_to_bar(ax4,auc_ecc_CH,[.2 .65],[.2:.1:.6],lw,alp,col_ci);
ax4.XAxis.Visible           = 'on';

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

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/FIG_computer_dyad/raw/FIG_computer_dyad', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/FIG_computer_dyad/raw/FIG_computer_dyad', '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reported stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Accuracy [coherence pooled, within subject]
acc_better_sign_bool  	= auc_acc_pooled_CH > .5 & p_acc_pooled_CH < ( .05 / length(p_acc_pooled_CH));
acc_better_perc_sign    = sum(acc_better_sign_bool) / length(acc_better_sign_bool);
acc_worse_sign_bool     = auc_acc_pooled_CH < .5 & p_acc_pooled_CH < ( .05 / length(p_acc_pooled_CH));
acc_worse_perc_sign     = sum(acc_worse_sign_bool) / length(acc_worse_sign_bool);

% Eccentricity [coherence pooled, within subject]
ecc_better_sign_bool  	= auc_ecc_pooled_CH > .5 & p_ecc_pooled_CH < ( .05 / length(p_ecc_pooled_CH));
ecc_better_perc_sign  	= sum(ecc_better_sign_bool) / length(ecc_better_sign_bool);
ecc_worse_sign_bool  	= auc_ecc_pooled_CH < .5 & p_ecc_pooled_CH < ( .05 / length(p_ecc_pooled_CH));
ecc_worse_perc_sign    	= sum(ecc_worse_sign_bool) / length(ecc_worse_sign_bool);

% Test across subjects
median(auc_ecc_pooled_CH)
[p,h,stats] = signrank(auc_ecc_pooled_CH,.5)
median(auc_acc_pooled_CH)
[p,h,stats] = signrank(auc_acc_pooled_CH,.5)

% Accuracy coherence effect
sign_better_acc = auc_acc_CH > .5 & p_acc_CH < ( .05 / (size(p_acc_CH,1)*size(p_acc_CH,2)));
acc_better_perc_sign = sum(sign_better_acc) / size(sign_better_acc,1);
sign_worse_acc = auc_acc_CH < .5 & p_acc_CH < ( .05 / (size(p_acc_CH,1)*size(p_acc_CH,2)));
acc_worse_perc_sign = sum(sign_worse_acc) / size(sign_worse_acc,1);

n_better_zero = sum(sign_better_acc(:,1));
acc_median_zero = median(auc_acc_CH(:,1));

% Eccentricity coherence effect
sign_better_ecc = auc_ecc_CH > .5 & p_ecc_CH < ( .05 / (size(p_ecc_CH,1)*size(p_ecc_CH,2)));
ecc_better_perc_sign = sum(sign_better_ecc) / length(sign_better_ecc(:,1));
sign_worse_ecc = auc_ecc_CH < .5 & p_ecc_CH < ( .05 / (size(p_ecc_CH,1)*size(p_ecc_CH,2)));
ecc_worse_perc_sign = sum(sign_worse_ecc) / length(sign_worse_ecc(:,1));

acc_median_98 = median(auc_ecc_CH(:,end));

[round(min(ecc_better_perc_sign),3) round(max(ecc_better_perc_sign),3)]
[round(min(ecc_worse_perc_sign),2) round(max(ecc_worse_perc_sign),2)]

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

function ax = plotBar_coherence(ax,dat,p_mat,p_val,lb_fs,snr)

axes(ax); hold on

pos = (dat > .5) & (p_mat < p_val);
neg = (dat < .5) & (p_mat < p_val);
all = size(dat,1);

bp                          = bar( (sum(pos)/all) .*100);
bp.FaceColor                = [.4 .4 .4];

bn                          = bar( -(sum(neg)/all) .*100);
bn.FaceColor                = [.1 .1 .1];

ax.YLim                     = [-60 60];
ax.YLabel.String            = '[%]';
ax.FontSize                 = lb_fs;
ax.XTick                 	= [1 4 7];
ax.XTickLabel            	= round(snr([1 4 7]),2).*100;
ax.XTickLabelRotation      	= 0;
ax.XAxis.Visible         	= 'off';
ax.YLabel.String           	= '[%]';
ax.XLabel.String           	= 'Coherence [%]';
end


function ax = add_mean_to_bar(ax,dat,ylim,ytick,lw,alp,col_ci)

yyaxis right
[CI,~]                      = bootci(1000,{@mean,dat},'Alpha',0.01); % Boostrap confidence intervals
vec                         = 1:length(CI); 
x_spacing                   = [vec fliplr(vec)];
ci                          = [CI(1,:) fliplr(CI(2,:))];
fl                          = fill(x_spacing,ci,col_ci,'EdgeColor','none', 'FaceAlpha', alp); % Overlay confidence intervals
mpl                         = plot(1:size(dat,2),mean(dat), 'LineWidth', lw/1.5, 'LineStyle', '-', 'Color', col_ci.*2); % Plot mean curve
ax.YLim                     = ylim;
ax.YTick                    = ytick;
ax.YLabel.String            = 'AUC';
ax.YAxis(1).Color           = 'k';
ax.YAxis(2).Color           = col_ci.*2;

end

% function ax = plotAUROC(ax,dat,str,lb_fs,snr,alp,lw,col_dat,col_ci)
% 
% axes(ax); hold on
% 
% % Remove zero rows
% dat(sum(dat == 0,2) == size(dat,2),:) = [];
% 
% % Plot subject-wise data
% for iL = 1:size(dat,1)
%     pl(iL)               	= plot(dat(iL,:));
%     pl(iL).Color          	= [.5 .5 .5 alp];
%     pl(iL).LineWidth      	= lw/lw;
% end
% 
% % Boostrap confidence intervals
% nRep                        = 1000;
% [CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.01);
% 
% % Prepare filled area
% vec                         = 1:length(CI);
% x_spacing                   = [vec fliplr(vec)];
% ci                          = [CI(1,:) fliplr(CI(2,:))];
% 
% % Overlay confidence intervals
% fl                          = fill(x_spacing,ci,col_ci,'EdgeColor','none', 'FaceAlpha', alp);
% 
% % Plot mean curve
% pm                          = plot(mean(dat),'LineWidth', lw/1.5, 'Color', col_dat);
% 
% ax.XLim                  	= [1 length(snr)];
% ax.YLim                  	= [0 1];
% ax.FontSize               	= lb_fs;
% ax.YLabel.String           	= str;
% ax.XLabel.String           	= 'Coherence [%]';
% ax.XTick                 	= 1:7;
% ax.YTick                  	= [.25 .5 .75];
% ax.XTickLabel            	= round(snr,2).*100;
% ax.XTickLabelRotation      	= 0;
% ax.XAxis.Visible            = 'off';
% 
% % n                           = .05;
% % crit                        = [.5-n .5+n];
% lm                          = line([1 7],[.5 .5], 'Color', 'k', 'LineStyle', ':', 'LineWidth',lw/2);
% % ll                          = line([1 7],[crit(1) crit(1)], 'Color', 'k', 'LineStyle', '-.', 'LineWidth',lw/2);
% % lh                          = line([1 7],[crit(2) crit(2)], 'Color', 'k', 'LineStyle', '-.', 'LineWidth',lw/2);
% 
% end

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