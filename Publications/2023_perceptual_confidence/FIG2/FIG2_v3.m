% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/MATLAB/cbrewer/
addpath /Users/fschneider/Documents/MATLAB/palamedes1_10_9/Palamedes

close all
clear all

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_performance.mat')

% Import subject summary spreadsheet
pth                         = '/Volumes/T7_Shield/CPR_psychophysics/';      % Local hard drive
x                           = readtable([pth 'Subjects_summary.xlsx']);     % Spreadsheet
sbj_lst                     = x.Abbreviation;                               % Subject ID list
sbj_lst(cellfun(@isempty,sbj_lst)) = [];

nSample                     = 30;                                           % Time window size [samples]
nLag                        = 150;                                          % Cross-correlation lag                               

%% Convert relevant data to matrix for simplicity

cnt                         = 0;
for iSubj = 1:length(sbj_lst)
    if isempty(solo_perf{iSubj}.hir)
        continue
    end
    
    cnt                   	= cnt + 1;
    solo_hir(cnt,:)        	= solo_perf{iSubj}.hir;          % Hit rate
    solo_trg_score(cnt,:)  	= solo_perf{iSubj}.trg_mscore;   % Avg target scores
    solo_macc(cnt,:)       	= solo_perf{iSubj}.macc_trg;     % Avg accuracy
    solo_mecc(cnt,:)       	= solo_perf{iSubj}.mecc;         % Avg eccentricity
    
    try
        dyad_hir(cnt,:)       	= dyad_perf{iSubj}.hir;          % Hit rate
        dyad_trg_score(cnt,:)	= dyad_perf{iSubj}.trg_mscore;   % Avg target scores
        dyad_macc(cnt,:)     	= dyad_perf{iSubj}.macc_trg;     % Avg accuracy
        dyad_mecc(cnt,:)    	= dyad_perf{iSubj}.mecc;         % Avg eccentricity
    catch
        dyad_hir(cnt,:)       	= nan;
        dyad_trg_score(cnt,:)	= nan;
        dyad_macc(cnt,:)     	= nan;
        dyad_mecc(cnt,:)    	= nan;
    end
end

%% PLOT

% f                           = figure('units','normalized','position',[0 0 .5 1]);
f                           = figure('units','centimeters','position',[0 0 14 16]);
height                      = [linspace(.555, .056,4) .1];
colmn                       = linspace(.075, .85,4);
lb_fs                       = 8;
lg_fs                       = 8;
lw                          = 3;
frme_ms                     = 1000/120;
alp                         = .4;
dim                         = [.18 .15];
col_dat                     = [0 0 0];
col_ci                      = [.3 0 0];
snr                         = solo_perf{end}.carr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Hit rate raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1                         = axes('Position', [colmn(1) height(1) dim]); hold on
[ax1,pl]                    = plotData(ax1,solo_hir,true,lw,alp,col_dat,col_ci);
ax1.YLim                    = [10 70];
ax1.XLim                    = [1 size(solo_hir,2)];
ax1.XLabel.String           = 'Coherence [%]';
ax1.YLabel.String           = 'Hit rate [%]';
ax1.YLabel.Position(1)      = -.4;
ax1.XTick                   = 1:length(snr);
ax1.XTickLabel              = round(snr,2)*100;
ax1.FontSize                = lb_fs;
ax1.XTickLabelRotation      = 0;
ax1.XAxis.Visible           = 'off';
% ax1.XGrid                   = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Lag %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSubj = 1:size(solo_cr,2)
    solo_mlag_id{iSubj} = solo_cr{iSubj}.id;

    for iCoh = 1:length(snr)
        clear cIdx
        cIdx             	= solo_cr{iSubj}.coh == snr(iCoh);
        solo_mlag(iSubj,iCoh) = mean(solo_cr{iSubj}.lag(cIdx));
    end
end

ax2                         = axes('Position', [colmn(1) height(4) dim]); hold on
[ax2,pl]                    = plotData(ax2,solo_mlag ./ 1e3,false,lw,alp,col_dat,col_ci);
ax2.XLim                    = [1 size(solo_mlag,2)];
ax2.YLim                    = [.3 1];
ax2.XLabel.String           = 'Coherence [%]';
ax2.YLabel.String           = 'Lag [s]';
ax2.YLabel.Position(1)      = -.4;
ax2.XTick                   = 1:length(snr);
ax2.XTickLabel              = round(snr,2)*100;
ax2.FontSize                = lb_fs;
ax2.XTickLabelRotation      = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Avg accuracy raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3                         = axes('Position', [colmn(1) height(2) dim]); hold on
[ax3,pl]                    = plotData(ax3,solo_macc,true,lw,alp,col_dat,col_ci);
ax3.XLim                    = [1 size(solo_macc,2)];
ax3.YLim                    = [30 100];
ax3.XLabel.String           = 'Coherence [%]';
ax3.YLabel.String           = 'Accuracy [%]';
ax3.YLabel.Position(1)      = -.4;
ax3.XTick                   = 1:length(snr);
ax3.XTickLabel              = round(snr,2)*100;
ax3.FontSize                = lb_fs;
ax3.XTickLabelRotation      = 0;
ax3.XAxis.Visible           = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Avg eccentricity raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax4                         = axes('Position', [colmn(1) height(3) dim]); hold on
[ax4,pl]                    = plotData(ax4,solo_mecc,true,lw,alp,col_dat,col_ci);
ax4.XLim                    = [1 size(solo_mecc,2)];
ax4.YLim                    = [20 100];
ax4.XLabel.String           = 'Coherence [%]';
ax4.YLabel.String           = 'Eccentricity [%]';
ax4.YLabel.Position(1)      = -.4;
ax4.XTick                   = 1:length(snr);
ax4.XTickLabel              = round(snr,2)*100;
ax4.FontSize                = lb_fs;
ax4.XTickLabelRotation      = 0;
ax4.XAxis.Visible           = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Hit rate raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax5                         = axes('Position', [colmn(2) height(1) dim]); hold on
[ax5]                       = plotScatter(ax5, solo_hir, dyad_hir, snr, lb_fs, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Dyadic - Lag %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub_cnt = 0;
for iSubj = 1:size(dyad_cr,2)
    if isempty(dyad_cr{iSubj})
        continue
    end
    
    sub_cnt                     = sub_cnt+1;
    dyad_mlag_id{sub_cnt}     	= dyad_cr{iSubj}.id;
        
    for iCoh = 1:length(snr)
        clear cIdx
        cIdx                	= dyad_cr{iSubj}.coh == snr(iCoh);
        dyad_mlag(sub_cnt,iCoh)	= mean(dyad_cr{iSubj}.lag(cIdx));
    end
end

for iSubj = 1:size(dyad_mlag_id,2)
    idx                         = cellfun(@(x) strcmp(x,dyad_mlag_id{iSubj}),solo_mlag_id);
    solo_mlag_sorted(iSubj,:)   = solo_mlag(idx,:);
end

ax6                         = axes('Position', [colmn(2) height(4) dim]); hold on
ln                          = line([0 1000],[0 1000]);
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];
[ax6]                       = plotScatter(ax6, solo_mlag_sorted ./ 1e3, dyad_mlag ./ 1e3, snr, lb_fs, false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Avg accuracy raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax7                         = axes('Position', [colmn(2) height(2) dim]); hold on
[ax7]                       = plotScatter(ax7, solo_macc, dyad_macc, snr, lb_fs, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Avg eccentricity raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax8                         = axes('Position', [colmn(2) height(3) dim]); hold on
[ax8]                       = plotScatter(ax8, solo_mecc, dyad_mecc, snr, lb_fs, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBPLOT: Hit rate comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 0;
for iSubj = 1:size(dyad_perf,2)
    if ~isempty(dyad_perf{iSubj})
        idx                         = cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),lower(sbj_lst));        
        cnt                         = cnt+1;
        df_hir(cnt,:)               = solo_perf{idx}.hir - dyad_perf{iSubj}.hir;
        df_hir_pool(cnt,:)          = solo_perf{idx}.hir_pool - dyad_perf{iSubj}.hir_pool;
    end
end

ax13                      	= axes('Position', [colmn(3) height(1) dim]); hold on
% pt                        = patch([0 8 8 0], [0 0 25 25], [.7 .7 .7], 'FaceAlpha',.2, 'EdgeColor','none');
% pt                       	= patch([0 8 8 0], [-25 -25 0 0], [.3 .3 .3], 'FaceAlpha',.2, 'EdgeColor','none');
tx                        	= text(1.25,18, 'Dyadic > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(1.25,-18, 'Solo > Dyadic', 'FontSize', lb_fs, 'Color', 'k');

[ax13,pl]                   = plotData(ax13,df_hir,true,lw,alp,col_dat,col_ci);
lm                          = line([1 7],[0 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth',lw/2);
ax13.XLim                   = [1 size(df_hir,2)];
ax13.FontSize               = lb_fs;
ax13.YLabel.String          = 'Difference';
ax13.XAxis.Visible          = 'off';
ax13.YLim                   = [-25 25];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUROC solo - dyadic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:size(dyad_perf,2)
    if isempty(dyad_perf{j})
        ply_id{j} = nan;
    else
        ply_id{j} = dyad_perf{j}.id;
    end
end

cnt = 0;

for iSubj = 1:size(dyad_perf,2)
    
    d.sp                  	= solo_perf{iSubj};                             % Performance data
    d.sc                  	= solo_cr{iSubj};                               % Correlation data
    
    sIdx                    = cellfun(@(x) strcmp(x,solo_perf{iSubj}.id),ply_id);
    
    if sum(sIdx) == 0
        continue
    end
    
    d.dp                  	= dyad_perf{sIdx};
    d.dc                  	= dyad_cr{sIdx};
         
    if isempty(d.dc)
        continue
    end
    
    cnt = cnt+1;

    for iCoh = 1:length(snr)
        bl_ecc(cnt,iCoh)     	= d.sp.mecc(iCoh);
        bl_acc(cnt,iCoh)     	= d.sp.macc_trg(iCoh);
        bl_scr(cnt,iCoh)     	= d.sp.trg_mscore(iCoh);
        bl_hir(cnt,iCoh)     	= d.sp.hir(iCoh);

        auc_acc(cnt,iCoh)     	= getAUROC(d.sp.acc_trg{iCoh},d.dp.acc_trg{iCoh});
        auc_ecc(cnt,iCoh)     	= getAUROC(d.sp.ecc{iCoh},d.dp.ecc{iCoh});
        auc_score(cnt,iCoh)    	= getAUROC(d.sp.trg_score{iCoh},d.dp.trg_score{iCoh});
        
        auc_cc(cnt,iCoh)    	= getAUROC(d.sc.cc(snr(iCoh) == d.sc.coh),d.dc.cc(snr(iCoh) == d.dc.coh));
        auc_xcp(cnt,iCoh)    	= getAUROC(d.sc.posPk(snr(iCoh) == d.sc.coh),d.dc.posPk(snr(iCoh) == d.dc.coh));
        auc_xc(cnt,iCoh)    	= getAUROC(d.sc.maxR(snr(iCoh) == d.sc.coh),d.dc.maxR(snr(iCoh) == d.dc.coh));
        auc_lag(cnt,iCoh)    	= getAUROC(d.sc.lag(snr(iCoh) == d.sc.coh),d.dc.lag(snr(iCoh) == d.dc.coh));
      
        p_acc(cnt,iCoh)      	= ranksum(d.sp.acc_trg{iCoh},d.dp.acc_trg{iCoh});
%         p_acc(cnt,iCoh)      	= ranksum(d.sp.acc_state{iCoh},d.dp.acc_state{iCoh});
        p_ecc(cnt,iCoh)      	= ranksum(d.sp.ecc{iCoh},d.dp.ecc{iCoh});
        p_score(cnt,iCoh)       = ranksum(d.sp.trg_score{iCoh},d.dp.trg_score{iCoh});
        p_cc(cnt,iCoh)      	= ranksum(d.sc.cc(snr(iCoh) == d.sc.coh),d.dc.cc(snr(iCoh) == d.dc.coh));
        p_xcp(cnt,iCoh)      	= ranksum(d.sc.posPk(snr(iCoh) == d.sc.coh),d.dc.posPk(snr(iCoh) == d.dc.coh));
        p_xc(cnt,iCoh)       	= ranksum(d.sc.maxR(snr(iCoh) == d.sc.coh),d.dc.maxR(snr(iCoh) == d.dc.coh));
        p_lag(cnt,iCoh)      	= ranksum(d.sc.lag(snr(iCoh) == d.sc.coh),d.dc.lag(snr(iCoh) == d.dc.coh));
    end
    
    p_ecc_pooled(cnt)           = ranksum(cell2mat(d.sp.ecc'),cell2mat(d.dp.ecc'));
    auc_ecc_pooled(cnt)         = getAUROC(cell2mat(d.sp.ecc'),cell2mat(d.dp.ecc'));
    
    p_acc_pooled(cnt)           = ranksum(cell2mat(d.sp.acc_trg),cell2mat(d.dp.acc_trg));
    auc_acc_pooled(cnt)         = getAUROC(cell2mat(d.sp.acc_trg),cell2mat(d.dp.acc_trg));
end

ax9                       	= axes('Position', [colmn(3) height(4) dim]); hold on

ax10                     	= axes('Position', [colmn(3) height(2) dim]); hold on
tx                        	= text(1.25,.8, 'Dyadic > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(1.25,.2, 'Solo > Dyadic', 'FontSize', lb_fs, 'Color', 'k');

ax11                       	= axes('Position', [colmn(3) height(3) dim]); hold on
% pt                         	= patch([0 8 8 0], [.5 .5 1 1], [.75 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
% pt                         	= patch([0 8 8 0], [0 0 .5 .5], [.5 .5 .5], 'FaceAlpha',.1, 'EdgeColor','none');

ax9                     	= plotAUROC(ax9,auc_lag,'AUROC',lb_fs,snr,alp,lw,col_dat,col_ci);
ax9.XAxis.Visible         	= 'on';
ax9.YLim(2)                 = .93;
ax10                     	= plotAUROC(ax10,auc_acc,'AUROC',lb_fs,snr,alp,lw,col_dat,col_ci);
ax11                    	= plotAUROC(ax11,auc_ecc,'AUROC',lb_fs,snr,alp,lw,col_dat,col_ci);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUROC Significance testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_boundary                = .05 / (size(p_ecc,1) * size(p_ecc,2)); % Bonferroni correction

ax15                       	= axes('Position', [colmn(4) height(2) .1 dim(2)]); hold on
ax15                      	= plotBar(ax15, auc_acc, p_acc, sig_boundary, lb_fs, snr);

ax16                       	= axes('Position', [colmn(4) height(3) .1 dim(2)]); hold on
ax16                      	= plotBar(ax16, auc_ecc, p_ecc, sig_boundary, lb_fs, snr);

ax17                       	= axes('Position', [colmn(4) height(4) .1 dim(2)]); hold on
ax17                      	= plotBar(ax17, auc_lag, p_lag, sig_boundary, lb_fs, snr);
ax17.XAxis.Visible         	= 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG2/';
ax00                         = axes('Position',[0 0 1 1],'Visible','off');

titl_fs = 10;
ofs = 0;
text(colmn(1)+ofs,.71, 'Solo condition', 'Parent', ax00, 'FontSize', titl_fs, 'Color', 'k')
text(colmn(2)+ofs,.71, 'Solo vs Dyadic', 'Parent', ax00, 'FontSize', titl_fs, 'Color', 'k')
text(colmn(3)+ofs,.71, 'Effect size', 'Parent', ax00, 'FontSize', titl_fs, 'Color', 'k')
text(colmn(4)+ofs,.55, 'Stats', 'Parent', ax00, 'FontSize', titl_fs, 'Color', 'k')

print(f, [dest_dir '/FIG2'], '-r500', '-dpng');
print(f, [dest_dir '/FIG2'], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG2'], '-r500', '-depsc2', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBPLOT: Score comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                           = figure('units','centimeters','position',[0 0 17.2 5.2]);
sbjs                        = cellfun(@lower, sbj_lst, 'UniformOutput', false);

% pl_dim                      = dim*1.33;
row                         = .79;
cmap_coh                    = cool(size(snr,2));

% ax20                    	= axes('Position', [colmn(1) row pl_dim]); hold on
ax20                        = subplot(1,3,1); hold on

for iSubj = 1:length(dyad_perf)
    if isempty(dyad_perf{iSubj})
        continue
    end
    sidx = cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),sbjs);
    for iCoh = 1:length(snr)
        sc(iCoh)                    = scatter(solo_perf{sidx}.trg_mscore(iCoh),dyad_perf{iSubj}.trg_mscore(iCoh));
        sc(iCoh) .MarkerFaceColor   = cmap_coh(iCoh,:);
        sc(iCoh) .MarkerFaceAlpha   = alp;
        sc(iCoh) .MarkerEdgeColor   = 'none';
        lg_str{iCoh}               	= num2str(round(snr(iCoh),2)*100);
    end
end

lg0                         = legend(sc,lg_str,'Location','northwest','NumColumns', 1);
lg0.Box                     = 'off';
lg0.Position(1)             = .25;
lg0.Position(2)             = .1;

ln                          = line([0 1],[0 1]);
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];

ax20.FontSize               = lb_fs;
ax20.YLabel.String          = 'Dyad score [a.u.]';
ax20.XLabel.String          = 'Solo score [a.u.]';
ax20.XTickLabelRotation     = 0;
axis tight

% Implement gramm solution later!
% clear g
% g=gramm('x',xsolo,'y',ydyad,'color',[0 0 0]);
% g.facet_grid([],b);
% g.geom_point();
% g.stat_cornerhist('aspect',0.5);
% g.geom_abline();
% g.set_title('Visualize x-y with stat_cornerhist()');
% figure('Position',[100 100 800 600]);
% g.draw();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Cumulative score comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cnt = 0;
for iSubj = 1:length(dyad_perf)
    if isempty(dyad_perf{iSubj})
        continue
    end
    
    cnt                     = cnt+1;
    sidx                    = cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),sbjs);        
    mscore_solo(cnt)    	= mean(solo_perf{sidx}.score_norm);
    mscore_dyadic(cnt)   	= mean(dyad_perf{iSubj}.score_norm);
end

ax21                        = subplot(1,3,2); hold on
axis tight

ln                          = line([0 1],[0 1]);
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];

sc                          = scatter(mscore_solo,mscore_dyadic);
sc.MarkerFaceColor          = [.3 .3 .3];
sc.MarkerFaceAlpha          = .75;
sc.MarkerEdgeColor          = 'none';

ax21.FontSize               = lb_fs;
ax21.YLim                   = [0 .4];
ax21.XLim                   = [0 .4];
ax21.YTick                  = [0 .2 .4];
ax21.XTick                  = [0 .2 .4];
ax21.YLabel.String          = 'Dyad score [a.u.]';
ax21.XLabel.String          = 'Solo score [a.u.]';
ax21.XTickLabelRotation     = 0;

[pv,h]                      = signrank(mscore_solo,mscore_dyadic); % paired, two-sided test for the null hypothesis that x – y comes from a distribution with zero median
tx                          = text(.05,.35,{['p = ' num2str(round(pv,2))]});
tx.Color                    = [0 0 0];
tx.FontSize                 = 8;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUROC Significant subjects histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nBin                        = 15;
ax12                        = subplot(2,3,6); hold on
hs1                         = histogram(auc_ecc_pooled(p_ecc_pooled < .05/length(p_ecc_pooled)),nBin);
hs1.FaceColor               = [.3 .3 .3];
hs1.EdgeColor               = 'none';
hs1.FaceAlpha               = .5;
hs2                         = histogram(auc_ecc_pooled(p_ecc_pooled >= .05/length(p_ecc_pooled)),nBin);
edges                       = hs1.BinEdges;
hs2.BinEdges                = edges;
hs2.FaceColor               = [1 .3 .3];
hs2.EdgeColor               = 'none';
hs2.FaceAlpha               = .5;
ax12.FontSize               = lb_fs;
ax12.XTick                  = [.3 .5 .7];
ax12.XLim                   = [.25 .75];
ax12.YLim                   = [0 5];
ax12.YLabel.String          = '# Subjects';
ax12.XLabel.String          = 'AUROC';
ln                          = line([.5 .5],[0 6]);
ln.LineWidth                = 1.5;
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];
ax12.Position               = [ax20.Position(1)+.6 ax20.Position(2) .15 .3];

ax13                        = subplot(2,3,3); hold on
hs2                         = histogram(auc_acc_pooled(p_acc_pooled >= .05/length(p_acc_pooled)),nBin);
hs2.BinEdges                = edges;
hs2.FaceColor               = [1 .3 .3];
hs2.EdgeColor               = 'none';
hs2.FaceAlpha               = .5;
hs1                         = histogram(auc_acc_pooled(p_acc_pooled < .05/length(p_acc_pooled)),nBin);
hs1.BinEdges                = edges;
hs1.FaceColor               = [.3 .3 .3];
hs1.EdgeColor               = 'none';
hs1.FaceAlpha               = .5;
ax13.FontSize               = lb_fs;
ax13.XTick                  = [.3 .5 .7];
ax13.XLim                   = [.25 .75];
ax13.YLim                   = [0 5];
ax13.YLabel.String          = '# Subjects';
ax13.XLabel.String          = 'AUROC';
ln                          = line([.5 .5],[0 6]);
ln.LineWidth                = 1.5;
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];
ax13.Position               = [ax20.Position(1)+.6 ax20.Position(2)+.5 .15 .3];

% PRINT
print(f, [dest_dir '/FIG2_a'], '-r500', '-dpng');
print(f, [dest_dir '/FIG2_a'], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG2_a'], '-r500', '-depsc2', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBPLOT: AUC - Solo scatter %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

figure(1);hold on
figure(2);hold on
coh_col                     = cool(length(snr));

for iPlot = 1:4
    clear xdat
    if iPlot == 1
        xdat = bl_ecc;
        str = 'Eccentricity';
    elseif iPlot == 2
        xdat = bl_acc;
        str = 'Accuracy'; 
     elseif iPlot == 3
        xdat = bl_scr;
        str = 'Hit score';  
            elseif iPlot == 4
        xdat = bl_hir;
        str = 'Hit rate'; 
    end
    
    figure(1)
    ax = subplot(2,2,iPlot); hold on
    ln = line([0 1],[.5 .5], 'Color','k', 'LineStyle', ':', 'LineWidth',1);
    
    for iCoh = 1:7
        scatter(xdat(:,iCoh),auc_ecc(:,iCoh),'MarkerFaceColor', coh_col(iCoh,:),'MarkerEdgeColor','none')
    end
    
    xlabel({str ' [Solo]'});
    ylabel({'Eccentricity AUC'; '[Solo vs Dyadic]'});
    xlim([0 1])
    ylim([0 1])
    set(ax,'fontsize', 14)
    
    figure(2)
    ax = subplot(2,2,iPlot); hold on
    
    for iCoh = 1:7
        [~,idx_solo] = sort(xdat(:,iCoh));
        [~,idx_auc] = sort(auc_ecc(:,iCoh));
        mat = [(1:34)' idx_solo idx_auc];
        for j = 1:34
            xx(j) = find(mat(:,2) == j);
            yy(j) = find(mat(:,3) == j);
        end
        scatter(xx,yy,'MarkerFaceColor',coh_col(iCoh,:),'MarkerEdgeColor','none')
    end
    
    xlabel({str ' [Solo rank]'});
    ylabel('Rank AUC');
    set(ax,'fontsize', 14)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reported stats in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hit rate drop for highest coherence level
for iSubj = 1:length(solo_cr)
    hir             = solo_perf{iSubj}.hir;
    hir_drop(iSubj)	= hir(end-1) > hir(end);
end
n                   = sum(hir_drop);
rate                = n / length(hir_drop);

% Dyad vs solo scores comparison
for iSubj = 1:size(solo_perf,2)
    solo_id{iSubj}  = solo_perf{iSubj}.id;
end

cnt = 0;
for iSubj = 1:length(dyad_perf)
    if ~isempty(dyad_perf{iSubj})
        cnt      	= cnt+1;
        idx        	= cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),solo_id);
        d_score     = mean(dyad_perf{iSubj}.score_norm);
        s_score     = mean(solo_perf{idx}.score_norm);
        
        d_better(cnt) 	= (d_score - s_score) > 0;
    end
end
n                   = sum(d_better);
rate                = n / length(d_better);
[pv,h]            	= signrank(mscore_solo,mscore_dyadic); % paired, two-sided test for the null hypothesis that x – y comes from a distribution with zero median

% Accuracy [coherence pooled, within subject]
sign_bool           = p_acc_pooled < ( .05 / length(p_acc_pooled));
perc_sign           = sum(sign_bool) / length(sign_bool);

% Eccentricity [coherence pooled, within subject]
sign_bool           = p_ecc_pooled < ( .05 / length(p_ecc_pooled));
perc_sign           = sum(sign_bool) / length(sign_bool);

% Coherence-wise: Eccentricity
sig_boundary      	= .05 / (size(p_ecc,1) * size(p_ecc,2)); % Bonferroni correction
pos = (auc_ecc > .5) & (p_ecc < sig_boundary);
neg = (auc_ecc < .5) & (p_ecc < sig_boundary);

percent_pos = (sum(pos)./size(auc_ecc,1)) .*100;
[round(min(percent_pos)) round(max(percent_pos))]
percent_neg = (sum(neg)./size(auc_ecc,1)) .*100;
[round(min(percent_neg)) round(max(percent_neg))]

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
[CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.05);

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

function [ax] = plotScatter(ax, solo_dat, dyad_dat, snr, lb_fs, ax_flag)

hold on
coh_col                     = cool(length(snr));

for iCoh = 1:size(solo_dat,2)
    sc                      = scatter(solo_dat(:,iCoh),dyad_dat(:,iCoh));
    sc.MarkerFaceColor      = coh_col(iCoh,:);
    sc.MarkerEdgeColor      = 'none';
    sc.MarkerFaceAlpha      = .4;
    sc.SizeData             = 20;
end

ax.FontSize                 = lb_fs;
ax.XLabel.String            = 'Solo';
ax.YLabel.String            = 'Dyad';
ax.XLim                 	= [0 1];
ax.YLim                     = [0 1];
ax.XTick                    = [.1 .5 .9];
ax.YTick                    = [.1 .5 .9];

if ax_flag 
    ax.XAxis.Visible        = 'off';
end

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
[CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.05);

% Prepare filled area
vec                         = 1:length(CI);
x_spacing                   = [vec fliplr(vec)];
ci                          = [CI(1,:) fliplr(CI(2,:))];

% Overlay confidence intervals
fl                          = fill(x_spacing,ci,col_ci,'EdgeColor','none', 'FaceAlpha', alp);

% Plot mean curve
pm                          = plot(mean(dat),'LineWidth', lw/1.5, 'Color', col_dat);

ax.XLim                  	= [1 length(snr)];
ax.YLim                  	= [.1 .9];
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

function ax = plotBar(ax,dat,p_mat,p_val,lb_fs,snr)

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