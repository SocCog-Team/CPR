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
        auc_acc_SC(scnt,iCoh)   	= getAUROC(solo_perf{iSub}.acc_trg{iCoh},hc_dyad_perf{idx_pc}.acc_trg{iCoh});
        auc_ecc_SC(scnt,iCoh)     	= getAUROC(solo_perf{iSub}.ecc_state{iCoh},hc_dyad_perf{idx_pc}.ecc_state{iCoh});
        auc_score_SC(scnt,iCoh)  	= getAUROC(solo_perf{iSub}.trg_score{iCoh},hc_dyad_perf{idx_pc}.trg_score{iCoh});
        
        % Solo - Human dyad
        auc_acc_SH(scnt,iCoh)   	= getAUROC(solo_perf{iSub}.acc_trg{iCoh},dyad_perf{idx_dy}.acc_trg{iCoh});
        auc_ecc_SH(scnt,iCoh)     	= getAUROC(solo_perf{iSub}.ecc_state{iCoh},dyad_perf{idx_dy}.ecc_state{iCoh});
        auc_score_SH(scnt,iCoh)  	= getAUROC(solo_perf{iSub}.trg_score{iCoh},dyad_perf{idx_dy}.trg_score{iCoh});
        
        % Human dyad - Computer dyad
        auc_acc_CH(scnt,iCoh)       = getAUROC(hc_dyad_perf{idx_pc}.acc_trg{iCoh},dyad_perf{idx_dy}.acc_trg{iCoh});
        auc_ecc_CH(scnt,iCoh)     	= getAUROC(hc_dyad_perf{idx_pc}.ecc_state{iCoh},dyad_perf{idx_dy}.ecc_state{iCoh});
        auc_score_CH(scnt,iCoh)  	= getAUROC(hc_dyad_perf{idx_pc}.trg_score{iCoh},dyad_perf{idx_dy}.trg_score{iCoh});
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
%%% SUBPLOT: Hit rate difference Computer-Human %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax2                         = axes('Position', [clmns(2) height(2) dim]); hold on
[ax2,pl]                    = plotData(ax2,hir_df_SC,true,lw,alp,col_dat,col_ci);
lm                          = line([1 7],[0 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth',lw/2);
ax2.XLim                    = [1 size(hir_df_SC,2)];
ax2.FontSize                = lb_fs;
ax2.YLabel.String           = 'Difference';
ax2.XAxis.Visible           = 'off';
ax2.YLim                    = [-20 40];
tx                        	= text(2.25,35, 'HC dyad > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(2.25,-15, 'Solo > HC dyad', 'FontSize', lb_fs, 'Color', 'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: AUROC Accuracy Computer-Human %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3                         = axes('Position', [clmns(2) height(3) dim]); hold on
ax3                         = plotAUROC(ax3,auc_acc_SC,'AUC',lb_fs,snr,alp,lw,col_dat,col_ci);
tx                        	= text(2.25,.9, 'HC dyad > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(2.25,.1, 'Solo > HC dyad', 'FontSize', lb_fs, 'Color', 'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: AUROC eccentricity Computer-Human %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax4                         = axes('Position', [clmns(2) height(4) dim]); hold on
ax4                         = plotAUROC(ax4,auc_ecc_SC,'AUC',lb_fs,snr,alp,lw,col_dat,col_ci);
ax4.XAxis.Visible           = 'on';
tx                        	= text(2.25,.9, 'HC dyad > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(2.25,.1, 'Solo > HC dyad', 'FontSize', lb_fs, 'Color', 'k');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Export
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG_computer_dyad/SFIG_computer_dyad/raw/SFIGraw_part2', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG_computer_dyad/SFIG_computer_dyad/raw/SFIGraw_part2', '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Igor - control plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([source_pth 'hh_dyad_pairwise_performance.mat'])
tmp_dyad = [dyad_pw_perf(:,1); dyad_pw_perf(:,2)];
for iSubj = 1:length(solo_perf)
        id_solo{iSubj}           	= solo_perf{iSubj}.id;
end

for iSubj = 1:size(tmp_dyad,1)
    for iCoh = 1:length(snr)
        idx_solo = cellfun(@(x) strcmp(x,tmp_dyad{iSubj}.id), id_solo);

        % Solo - Computer dyad
        auc_acc_SH_pw(iSubj,iCoh) = getAUROC(solo_perf{idx_solo}.acc_trg{iCoh},tmp_dyad{iSubj}.acc_trg{iCoh});
        auc_ecc_SH_pw(iSubj,iCoh) = getAUROC(solo_perf{idx_solo}.ecc_trg{iCoh},tmp_dyad{iSubj}.ecc_trg{iCoh});
      
    end
end

figure; hold on
cl = [0 0 1;.8 0 0];
for ii = 1:2
    clear xdat ydat
    
    if ii == 1
        xdat = mean(auc_ecc_SC,2);
        ydat = mean(auc_acc_SC,2);
        str = 'Solo vs HC Dyad';
        ofs = .1;
    else
        xdat = mean(auc_acc_SH_pw,2);
        ydat = mean(auc_ecc_SH_pw,2);
        str = 'Solo vs HH Dyad [pairwise]';
        ofs = .6;
    end
    
    size(xdat)
    sc = scatter(xdat, ydat, 'MarkerFaceColor', cl(ii,:),'MarkerFaceAlpha',.3, 'MarkerEdgeColor','none'); lsline
    [r,p] = corrcoef(xdat, ydat);
    text(ofs,.35, str, 'FontSize', 12, 'Color', cl(ii,:));
    text(ofs,.25, ['r = ' num2str(r(2))], 'FontSize', 12, 'Color', cl(ii,:));
    text(ofs,.3, ['p = ' num2str(p(2))], 'FontSize', 12, 'Color', cl(ii,:));
end

line([.5 .5], [0 1],'Color', 'k', 'LineStyle', ':', 'LineWidth',1.5)
line([0 1], [.5 .5],'Color', 'k', 'LineStyle', ':', 'LineWidth',1.5);
xlabel('AUC: Eccentricity')
ylabel('AUC: Accuracy')
ylim([0 1])
xlim([0 1])
set(gca,'fontsize', 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    if sum(idx_pc>0)
        for iCoh = 1:length(snr)
            
            y_mat(iSubj,iCoh)   = in_pc{idx_pc}.(cond_str)(iCoh);
            x_mat(iSubj,iCoh)   = in_solo{iSubj}.(cond_str)(iCoh);
            sc(iCoh)            = scatter(in_solo{iSubj}.(cond_str)(iCoh), in_pc{idx_pc}.(cond_str)(iCoh),'MarkerFaceColor', coh_col(iCoh,:)./2,'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .3);
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
    x_ci            = (bootci(500, {@nanmedian,  xx},'alpha', .01));
    y_ci            = (bootci(500, {@nanmedian,  yy},'alpha', .01));
    lny             = line([nanmedian(xx) nanmedian(xx)],[y_ci(1) y_ci(2)], 'Color', coh_col(iCoh,:),'LineWidth',1);
    lnx             = line([x_ci(1) x_ci(2)],[nanmedian(yy) nanmedian(yy)], 'Color', coh_col(iCoh,:),'LineWidth',1);
    sc              = scatter(nanmedian(xx),nanmedian(yy), 'MarkerFaceColor', coh_col(iCoh,:),'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1,'SizeData', 20);
end

ax.FontSize                 = lb_fs;
ax.YLabel.String            = 'HC dyad';
ax.XLabel.String            = 'Solo';
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

function [out] = getAUROC(in1, in2)

if size(in1,1) == 1
    in1         = in1';
    in2         = in2';
end

lab          	= [zeros(length(in1),1); ones(length(in2),1)];
[~,~,~,out]     = perfcurve(lab,[in1; in2],1);

end