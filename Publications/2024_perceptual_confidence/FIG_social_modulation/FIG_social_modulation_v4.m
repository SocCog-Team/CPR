% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

% close all
clear all

source_dir = '/Users/fschneider/ownCloud/var_plot/';
load([ source_dir '/solo_correlation.mat'])
load([ source_dir '/solo_performance.mat'])
load([ source_dir '/hh_dyad_correlation.mat'])
load([ source_dir '/hh_dyad_performance.mat'])
load([ source_dir '/hh_dyad_pairwise_correlation.mat'])
load([ source_dir '/hh_dyad_pairwise_performance.mat'])

for iSubj = 1:size(solo_perf,2)
    sbj_lst{iSubj}          = solo_perf{iSubj}.id;                          % Extract subject ID
end

nSample                     = 30;                                           % Time window size [samples]
nLag                        = 150;                                          % Cross-correlation lag                               

state_alignment           	= true;                                         % Target-alignment == False

%% Convert relevant data to matrix for simplicity

cnt                         = 0;
for iSubj = 1:length(sbj_lst)
    if isempty(dyad_perf{iSubj})
        continue
    end
    
    cnt                   	= cnt + 1;
    dyad_hir(cnt,:)        	= dyad_perf{iSubj}.hir;                         % Hit rate
    
    if state_alignment == true
        dyad_macc(cnt,:)       	= dyad_perf{iSubj}.macc_state;              % Avg accuracy
        dyad_mecc(cnt,:)       	= dyad_perf{iSubj}.mecc_state;              % Avg eccentricity
    else
        dyad_mecc(cnt,:)       	= dyad_perf{iSubj}.mecc_trg;                % Avg eccentricity
        dyad_macc(cnt,:)       	= dyad_perf{iSubj}.macc_trg;               	% Avg accuracy
    end
end

%% PLOT

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
%%% SUBPLOT: Hit rate comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 0;
for iSubj = 1:size(dyad_perf,2)
    if ~isempty(dyad_perf{iSubj})
        idx                         = cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),lower(sbj_lst));        
        cnt                         = cnt+1;
        df_hir(cnt,:)               = dyad_perf{iSubj}.hir - solo_perf{idx}.hir;
        df_hir_pool(cnt,:)          = dyad_perf{iSubj}.hir_pool - solo_perf{idx}.hir_pool;
    end
end

ax13                      	= axes('Position', [colmn(2) height(1) dim]); hold on
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
%%% AUROC solo - dyadic
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
        if state_alignment == true
            bl_ecc(cnt,iCoh)     	= d.sp.mecc_state(iCoh);
            bl_acc(cnt,iCoh)     	= d.sp.macc_state(iCoh);
            [~,~,auc_acc(cnt,iCoh)] = getAUROC(d.sp.acc_state{iCoh},d.dp.acc_state{iCoh});
            [~,~,auc_ecc(cnt,iCoh)] = getAUROC(d.sp.ecc_state{iCoh},d.dp.ecc_state{iCoh});
            p_acc(cnt,iCoh)      	= ranksum(d.sp.acc_state{iCoh},d.dp.acc_state{iCoh});
            p_ecc(cnt,iCoh)      	= ranksum(d.sp.ecc_state{iCoh},d.dp.ecc_state{iCoh});
        else
            bl_ecc(cnt,iCoh)     	= d.sp.mecc_trg(iCoh);
            bl_acc(cnt,iCoh)     	= d.sp.macc_trg(iCoh);
            [~,~,auc_acc(cnt,iCoh)] = getAUROC(d.sp.acc_trg{iCoh},d.dp.acc_trg{iCoh});
            [~,~,auc_ecc(cnt,iCoh)] = getAUROC(d.sp.ecc_trg{iCoh},d.dp.ecc_trg{iCoh});
            p_acc(cnt,iCoh)      	= ranksum(d.sp.acc_trg{iCoh},d.dp.acc_trg{iCoh});
            p_ecc(cnt,iCoh)      	= ranksum(d.sp.ecc_trg{iCoh},d.dp.ecc_trg{iCoh});
        end
            
        bl_scr(cnt,iCoh)            = d.sp.trg_mscore(iCoh);
        bl_hir(cnt,iCoh)            = d.sp.hir(iCoh);
        bl_lag(cnt,iCoh)            = mean(d.sc.lag(d.sc.coh == snr(iCoh)));

        [~,~,auc_score(cnt,iCoh)]	= getAUROC(d.sp.trg_score{iCoh},d.dp.trg_score{iCoh});
        [~,~,auc_cc(cnt,iCoh)]    	= getAUROC(d.sc.cc(snr(iCoh) == d.sc.coh),d.dc.cc(snr(iCoh) == d.dc.coh));
        [~,~,auc_xcp(cnt,iCoh)]    	= getAUROC(d.sc.posPk(snr(iCoh) == d.sc.coh),d.dc.posPk(snr(iCoh) == d.dc.coh));
        [~,~,auc_xc(cnt,iCoh)]    	= getAUROC(d.sc.maxR(snr(iCoh) == d.sc.coh),d.dc.maxR(snr(iCoh) == d.dc.coh));
        [~,~,auc_lag(cnt,iCoh)]    	= getAUROC(d.sc.lag(snr(iCoh) == d.sc.coh),d.dc.lag(snr(iCoh) == d.dc.coh));

        p_score(cnt,iCoh)           = ranksum(d.sp.trg_score{iCoh},d.dp.trg_score{iCoh});
        p_cc(cnt,iCoh)              = ranksum(d.sc.cc(snr(iCoh) == d.sc.coh),d.dc.cc(snr(iCoh) == d.dc.coh));
        p_xcp(cnt,iCoh)             = ranksum(d.sc.posPk(snr(iCoh) == d.sc.coh),d.dc.posPk(snr(iCoh) == d.dc.coh));
        p_xc(cnt,iCoh)              = ranksum(d.sc.maxR(snr(iCoh) == d.sc.coh),d.dc.maxR(snr(iCoh) == d.dc.coh));
        p_lag(cnt,iCoh)             = ranksum(double(d.sc.lag(snr(iCoh) == d.sc.coh)),double(d.dc.lag(snr(iCoh) == d.dc.coh)));
    end
    
    if state_alignment == true
        p_ecc_pooled(cnt)           = ranksum(cell2mat(d.sp.ecc_state'),cell2mat(d.dp.ecc_state'));
        [~,~,auc_ecc_pooled(cnt)]  	= getAUROC(cell2mat(d.sp.ecc_state'),cell2mat(d.dp.ecc_state'));
        
        p_acc_pooled(cnt)           = ranksum(cell2mat(d.sp.acc_state'),cell2mat(d.dp.acc_state'));
        [~,~,auc_acc_pooled(cnt)]  	= getAUROC(cell2mat(d.sp.acc_state'),cell2mat(d.dp.acc_state'));
    else
        p_ecc_pooled(cnt)           = ranksum(cell2mat(d.sp.ecc_trg),cell2mat(d.dp.ecc_trg));
        [~,~,auc_ecc_pooled(cnt)]  	= getAUROC(cell2mat(d.sp.ecc_trg),cell2mat(d.dp.ecc_trg));
        
        p_acc_pooled(cnt)           = ranksum(cell2mat(d.sp.acc_trg),cell2mat(d.dp.acc_trg));
        [~,~,auc_acc_pooled(cnt)]  	= getAUROC(cell2mat(d.sp.acc_trg),cell2mat(d.dp.acc_trg));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: AUROC per coherence %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax10                     	= axes('Position', [colmn(2) height(2) dim]); hold on
ax10                     	= plotAUROC(ax10,auc_acc,'AUC',lb_fs,snr,alp,lw,col_dat,col_ci);

tx                        	= text(1.25,.8, 'Dyadic > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(1.25,.2, 'Solo > Dyadic', 'FontSize', lb_fs, 'Color', 'k');

ax11                       	= axes('Position', [colmn(2) height(3) dim]); hold on
ax11                    	= plotAUROC(ax11,auc_ecc,'AUC',lb_fs,snr,alp,lw,col_dat,col_ci);
ax11.XAxis.Visible         	= 'on';

tx                        	= text(1.25,.8, 'Dyadic > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(1.25,.2, 'Solo > Dyadic', 'FontSize', lb_fs, 'Color', 'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC AUC - Solo hit rate  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax5                         = axes('Position', [colmn(3) height(1) dim]); hold on
[ax5,mdl,R,PV]              = plotQuartiles(ax5,bl_hir,df_hir,'Hit rate dff',lw,lb_fs,false);
ax5.YLim                    = [-.1 .1];
ax5.YTick                   = [-.1 0 .1];
ax5.XAxis.Visible         	= 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC AUC - Solo accuracy %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax7                         = axes('Position', [colmn(3) height(2) dim]); hold on
[ax7,mdl,R,PV]              = plotQuartiles(ax7,bl_acc,auc_acc,'Accuracy [AUC]',lw,lb_fs, true);
ax7.YLim                    = [.4 .6];
ax7.YTick                   = [.4 .5 .6];
ax7.XAxis.Visible         	= 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC AUC - Solo eccentricity %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax8                         = axes('Position', [colmn(3) height(3) dim]); hold on
[ax8,mdl,R,PV]              = plotQuartiles(ax8,bl_ecc,auc_ecc,'Eccentricity [AUC]',lw,lb_fs, true);
ax8.YLim                    = [.25 .75];
ax8.YTick                   = [.3 .5 .7];
ax8.XAxis.Visible         	= 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AUROC Significance testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_boundary                = .05 / (size(p_ecc,1) * size(p_ecc,2)); % Bonferroni correction

ax15                       	= axes('Position', [colmn(4) height(2) .1 dim(2)]); hold on
ax15                      	= plotBar(ax15, auc_acc, p_acc, sig_boundary, lb_fs, snr);
ax15                        = add_mean_to_bar(ax15,auc_acc,[.4 .6],[.4 .5 .6],lw,alp,col_ci);

ax16                       	= axes('Position', [colmn(4) height(3) .1 dim(2)]); hold on
ax16                      	= plotBar(ax16, auc_ecc, p_ecc, sig_boundary, lb_fs, snr);
ax16.XAxis.Visible         	= 'on';
ax16                        = add_mean_to_bar(ax16,auc_ecc,[.25 .75],[.3 .5 .7],lw,alp,col_ci);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/FIG_social_modulation/raw/';
ax00                     	= axes('Position',[0 0 1 1],'Visible','off');
titl_fs                     = 10;
ofs                         = 0;

text(colmn(2)+ofs,.71, 'Solo vs Dyadic', 'Parent', ax00, 'FontSize', titl_fs, 'Color', 'k')
text(colmn(3)+ofs,.71, 'Effect size vs Solo', 'Parent', ax00, 'FontSize', titl_fs, 'Color', 'k')
text(colmn(4)+ofs,.55, 'Stats', 'Parent', ax00, 'FontSize', titl_fs, 'Color', 'k')

if state_alignment == true
    print(f, [dest_dir '/FIG_social_modulation_state'], '-r500', '-dpng');
    print(f, [dest_dir '/FIG_social_modulation_state'], '-r500', '-dsvg', '-painters');
else
    print(f, [dest_dir '/FIG_social_modulation_target'], '-r500', '-dpng');
    print(f, [dest_dir '/FIG_social_modulation_target'], '-r500', '-dsvg', '-painters');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UPPER PANEL: Score comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                           = figure('units','centimeters','position',[0 0 17.2 5.2]);
sbjs                        = cellfun(@lower, sbj_lst, 'UniformOutput', false);
row                         = .79;
cmap_coh                    = cool(size(snr,2));
ax20                        = subplot(1,3,1); hold on

% Get scores
for iSubj = 1:length(dyad_perf)
    if isempty(dyad_perf{iSubj})
        x_mat(iSubj,:) = nan([1 7]);
        y_mat(iSubj,:) = nan([1 7]);
        continue
    end
    
    sidx                    = cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),sbjs);
    
    for iCoh = 1:length(snr)
        x_mat(iSubj,iCoh)           = mean(solo_perf{sidx}.trg_all_score(solo_perf{sidx}.trg_all_coh == snr(iCoh)));
        y_mat(iSubj,iCoh)           = mean(dyad_perf{sidx}.trg_all_score(dyad_perf{sidx}.trg_all_coh == snr(iCoh)));
    end
end

% Plot scatter plot
for iCoh = 1:length(snr)
    sc(iCoh)                    = scatter(x_mat(:,iCoh),y_mat(:,iCoh) );  
    sc(iCoh) .MarkerFaceColor   = cmap_coh(iCoh,:)/2;
    sc(iCoh) .MarkerFaceAlpha   = .3;
    sc(iCoh) .MarkerEdgeColor   = 'none';
    lg_str{iCoh}               	= num2str(round(snr(iCoh),2)*100);
end

for iCoh = 1:size(x_mat,2)
    xx                      = x_mat(:,iCoh);
    yy                      = y_mat(:,iCoh);
    x_ci                    = (bootci(500, {@nanmedian,  xx},'alpha', .001));
    y_ci                    = (bootci(500, {@nanmedian,  yy},'alpha', .001));
    lny                     = line([nanmedian(xx) nanmedian(xx)],[y_ci(1) y_ci(2)], 'Color', cmap_coh(iCoh,:),'LineWidth',1);
    lnx                     = line([x_ci(1) x_ci(2)],[nanmedian(yy) nanmedian(yy)], 'Color', cmap_coh(iCoh,:),'LineWidth',1);
    sc                      = scatter(nanmedian(xx),nanmedian(yy), 'MarkerFaceColor', cmap_coh(iCoh,:),'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1,'SizeData', 20);
end

ln                          = line([0 .5],[0 .5]);
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];

ax20.FontSize               = lb_fs;
ax20.YLabel.String          = 'Dyadic';
ax20.XLabel.String          = 'Solo';
ax20.XTickLabelRotation     = 0;
ax20.XTick                  = 0:.1:.5;
ax20.YTick                  = 0:.1:.5;
axis tight

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
    
    mscore_solo(cnt)    	= mean(solo_perf{sidx}.trg_all_score);
    mscore_dyadic(cnt)   	= mean(dyad_perf{iSubj}.trg_all_score);    
%     mscore_solo(cnt)    	= mean(solo_perf{sidx}.score_norm);
%     mscore_dyadic(cnt)   	= mean(dyad_perf{iSubj}.score_norm);
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

linearCoefficients          = polyfit(mscore_solo, mscore_dyadic, 1);
x_fit                       = linspace(0, 1, 50);
y_fit                       = polyval(linearCoefficients, x_fit);
pl                          = plot(x_fit, y_fit, '-', 'LineWidth', 1,'Color',[.6 .1 .1 alp]);

ax21.FontSize               = lb_fs;
ax21.YLim                   = [0 .4];
ax21.XLim                   = [0 .4];
ax21.YTick                  = [0:.1:.4];
ax21.XTick                  = [0:.1:.4];
ax21.YLabel.String          = 'Dyadic';
ax21.XLabel.String          = 'Solo';
ax21.XTickLabelRotation     = 0;

[pv,h,z]                    = signrank(mscore_solo,mscore_dyadic); % paired, two-sided test for the null hypothesis that x – y comes from a distribution with zero median
tx                          = text(.05,.35,{['p = ' num2str(round(pv,2))]});
tx.Color                    = [0 0 0];
tx.FontSize                 = 8;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUROC Significant subjects histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hist_sig                    = [.6 .1 .1];
hist_nonsig                 = [.1 .1 .1];
nBin                        = 17;

ax12                        = subplot(2,3,6); hold on
hs1                         = histogram(auc_ecc_pooled(p_ecc_pooled < .05/length(p_ecc_pooled)),nBin);
hs1.FaceColor               = hist_sig;
hs1.EdgeColor               = 'none';
hs1.FaceAlpha               = .5;
hs2                         = histogram(auc_ecc_pooled(p_ecc_pooled >= .05/length(p_ecc_pooled)),nBin);
edges                       = hs1.BinEdges;
hs2.BinEdges                = edges;
hs2.FaceColor               = hist_nonsig;
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
hs2.FaceColor               = hist_nonsig;
hs2.EdgeColor               = 'none';
hs2.FaceAlpha               = .5;
hs1                         = histogram(auc_acc_pooled(p_acc_pooled < .05/length(p_acc_pooled)),nBin);
hs1.BinEdges                = edges;
hs1.FaceColor               = hist_sig;
hs1.EdgeColor               = 'none';
hs1.FaceAlpha               = .5;
ax13.FontSize               = lb_fs;
ax13.XTick                  = [.3 .5 .7];
ax13.XLim                   = [.25 .75];
ax13.YLim                   = [0 12];
ax13.YLabel.String          = '# Subjects';
ax13.XLabel.String          = 'AUROC';
ln                          = line([.5 .5],[0 12]);
ln.LineWidth                = 1.5;
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];
ax13.Position               = [ax20.Position(1)+.6 ax20.Position(2)+.5 .15 .3];

% PRINT
if state_alignment == true
    print(f, [dest_dir '/FIG_social_modulation_top_state'], '-r500', '-dpng');
    print(f, [dest_dir '/FIG_social_modulation_top_state'], '-r500', '-dsvg', '-painters');
else
    print(f, [dest_dir '/FIG_social_modulation_top_target'], '-r500', '-dpng');
    print(f, [dest_dir '/FIG_social_modulation_top_target'], '-r500', '-dsvg', '-painters');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scatter plot for supplementary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                           = figure('units','centimeters','position',[0 0 22.5 5]); hold on
ax5                         = subplot(1,4,1); hold on
[ax5]                       = plotScatter(ax5, bl_hir, df_hir, snr, lb_fs, false, false);
ax7                         = subplot(1,4,2); hold on
[ax7]                       = plotScatter(ax7, bl_acc, auc_acc, snr, lb_fs, false, true);
ax8                         = subplot(1,4,3); hold on
[ax8]                       = plotScatter(ax8, bl_ecc, auc_ecc, snr, lb_fs, false, true);

if state_alignment == true 
    print(f, [dest_dir '/SFIG_social_modulation_scatter_state'], '-r500', '-dpng');
    print(f, [dest_dir '/SFIG_social_modulation_scatter_state'], '-r500', '-dsvg', '-painters');
else
    print(f, [dest_dir '/SFIG_social_modulation_scatter_target'], '-r500', '-dpng');
    print(f, [dest_dir '/SFIG_social_modulation_scatter_target'], '-r500', '-dsvg', '-painters');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Different grouping for supplementary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure;

% X: hit rate
ax = subplot(3,3,1); hold on
plotQuartiles(ax,bl_hir,df_hir,'Hit rate [Dff]',lw,lb_fs, false);
ax.YLim                     = [-.1 .1];
ax.Title.String             = 'X: Hit rate';
ax = subplot(3,3,4); hold on
plotQuartiles(ax,bl_hir,auc_acc,'Accuracy [AUC]',lw,lb_fs, true);
ax.YLim                    = [.4 .6];
ax = subplot(3,3,7); hold on
plotQuartiles(ax,bl_hir,auc_ecc,'Eccentricity [AUC]',lw,lb_fs, true);
ax.YLim                    = [.3 .7];

% X: accuracy
ax = subplot(3,3,2); hold on
plotQuartiles(ax,bl_acc,df_hir,'Hit rate [Dff]',lw,lb_fs, false);
ax.YLim                    = [-.1 .1];
ax.Title.String             = 'X: Accuracy';
ax = subplot(3,3,5); hold on

plotQuartiles(ax,bl_acc,auc_acc,'Accuracy [AUC]',lw,lb_fs, true);
ax.YLim                    = [.4 .6];
%%
ax = subplot(3,3,8); hold on
plotQuartiles(ax,bl_acc,auc_ecc,'Eccentricity [AUC]',lw,lb_fs, true);
ax.YLim                    = [.3 .7];

% X: eccentricity
ax = subplot(3,3,3); hold on
plotQuartiles(ax,bl_ecc,df_hir,'Hit rate [Dff]',lw,lb_fs, false);
ax.YLim                    = [-.1 .1];
ax.Title.String             = 'X: Eccentricity';
ax = subplot(3,3,6); hold on
plotQuartiles(ax,bl_ecc,auc_acc,'Accuracy [AUC]',lw,lb_fs, true);
ax.YLim                    = [.4 .6];
ax = subplot(3,3,9); hold on
plotQuartiles(ax,bl_ecc,auc_ecc,'Eccentricity [AUC]',lw,lb_fs, true);
ax.YLim                    = [.3 .7];

if state_alignment == true 
print(f, [dest_dir '/SFIG_social_modulation_matrix_state'], '-r500', '-dpng');
print(f, [dest_dir '/SFIG_social_modulation_matrix_state'], '-r500', '-dsvg', '-painters');
else
print(f, [dest_dir '/SFIG_social_modulation_matrix_target'], '-r500', '-dpng');
print(f, [dest_dir '/SFIG_social_modulation_matrix_target'], '-r500', '-dsvg', '-painters'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hit rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;hold on
% cnt = 0;
% for iSubj = 1:length(dyad_perf)
%     if isempty(dyad_perf{iSubj})
%         continue
%     end
%     
%     cnt                     = cnt+1;
%     sidx                    = cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),sbjs);  
%     
%     hir_solo(cnt)    	= mean(solo_perf{sidx}.hir);
%     hir_dyadic(cnt)   	= mean(dyad_perf{iSubj}.hir);    
% end
% 
% axis tight
% 
% ln                          = line([0 1],[0 1]);
% ln.LineStyle                = ':';
% ln.Color                    = [0 0 0];
% 
% sc                          = scatter(hir_solo,hir_dyadic);
% sc.MarkerFaceColor          = [.3 .3 .3];
% sc.MarkerFaceAlpha          = .75;
% sc.MarkerEdgeColor          = 'none';
% 
% 
% linearCoefficients          = polyfit(hir_solo, hir_dyadic, 1);
% x_fit                       = linspace(0, 1, 50);
% y_fit                       = polyval(linearCoefficients, x_fit);
% pl                          = plot(x_fit, y_fit, '-', 'LineWidth', 1,'Color',[.6 .1 .1 alp]);
% 
% ax21 = gca
% ax21.FontSize               = 16;
% ax21.YLim                   = [.2 .6];
% ax21.XLim                   = [.2 .6];
% ax21.YTick                  = [.2:.1:.6];
% ax21.XTick                  = [.2:.1:.6];
% ax21.YLabel.String          = 'Dyadic';
% ax21.XLabel.String          = 'Solo';
% ax21.XTickLabelRotation     = 0;
% 
% [pv,h,z]                    = signrank(hir_solo,hir_dyadic); % paired, two-sided test for the null hypothesis that x – y comes from a distribution with zero median
% tx                          = text(.25,.55,{['p = ' num2str(round(pv,2))]});
% tx.Color                    = [0 0 0];
% tx.FontSize                 = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Performance quartile - social modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Population %%%
% clear bl dat quartile_boundaries out
% ax                          = subplot(1,3,3); hold on
% ccol                        = [.5 .5 .5; 0 0 0; .8 .5 .5; .5 .1 .1];
% 
% for iDat = 1:4
%     if iDat == 1
%         dat                 = mean(auc_acc,2);
%         bl                  = mean(bl_acc,2);
%     elseif iDat == 2
%         dat                 = mean(auc_ecc,2);
%         bl                  = mean(bl_ecc,2);
%     elseif iDat == 3
%         dat                 = mean(auc_acc,2);
%         bl                  = mean(bl_ecc,2);
%     elseif iDat == 4
%         dat                 = mean(auc_ecc,2);
%         bl                  = mean(bl_acc,2);
%     end
%     
%     mdl_comb{iDat}                  = fitlm(bl,dat);
%     [R_comb{iDat},PV_comb{iDat}]    = corrcoef(bl,dat); 
%     
%     quartile_boundaries   	= prctile(bl, [25 50 75]);
%     
%     for iQuart = 1:4
%         if iQuart == 1
%             idx             = bl < quartile_boundaries(1);
%         elseif iQuart == 2 || iQuart == 3
%             idx             = bl  >= quartile_boundaries(iQuart-1) & bl  < quartile_boundaries(iQuart);
%         elseif iQuart == 4
%             idx             = bl  > quartile_boundaries(3);
%         end
%         out(iQuart)         = mean(dat(idx));
%     end
%     
%     ln                      = line([0 5],[.5 .5], 'Color', [0 0 0],'LineWidth',1,'LineStyle',':');
% 
%     plt(iDat)               = plot(out,'Color',[ccol(iDat,:)],'LineWidth', 2);
%     ax                      = gca;
%     ax.XLim                 = [0 5];
%     ax.XTick                = [1:4];
%     ax.XTickLabel           = {'Bottom','2nd','3rd','Top'};
%     ax.XLabel.String      	= 'Solo performance [Quartiles]';
%     ax.YLabel.String      	= 'Avg. social modulation [AUC]';
%     ax.FontSize             = lb_fs;
%     ax.XTickLabelRotation   = 0;
% end
% 
% legend(plt,'x:Acc y:Acc', 'x:Ecc y:Ecc', 'x:Ecc y:Acc','x:Acc y:Ecc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bl                      = mean(bl_acc,2);
% mconf                   = mean(bl_ecc,2);
% quartile_boundaries   	= prctile(bl, [25 50 75]);
% 
% for iQuart = 1:4
%     if iQuart == 1
%         idx             = bl < quartile_boundaries(1);
%     elseif iQuart == 2 || iQuart == 3
%         idx             = bl  >= quartile_boundaries(iQuart-1) & bl  < quartile_boundaries(iQuart);
%     elseif iQuart == 4
%         idx             = bl  > quartile_boundaries(3);
%     end
%     conf(iQuart)     	= mean(mconf(idx));
%     cpmt(iQuart)      	= mean(bl(idx));
% end
% 
% figure; hold on
% pl(1)               	= plot(conf,'Color',[.6 .1 .1],'LineWidth', 2);
% pl(2)                 	= plot(cpmt,'Color',[.1 .1 .1],'LineWidth', 2);
% ax                      = gca;
% ax.XLim                 = [.75 4.25];
% ax.XTick                = [1:4];
% ax.XTickLabel           = {'Bottom','2nd','3rd','Top'};
% ax.FontSize             = 24;
% ax.YLabel.String      	= 'Solo [norm]';
% ax.XLabel.String      	= 'Accuracy';
% legend(pl,{'Confidence','Accuracy'},'Location', 'southeast')

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
        d_score(cnt)= mean(dyad_perf{iSubj}.score_norm);
        s_score(cnt)= mean(solo_perf{idx}.score_norm);
        
        d_better(cnt) 	= (d_score(cnt) - s_score(cnt)) > 0;
    end
end
n                   = sum(d_better);
rate                = n / length(d_better);
[pv,h,stats]     	= signrank(d_score,s_score); % paired, two-sided test for the null hypothesis that x – y comes from a distribution with zero median
[mean(s_score) iqr(s_score)]
[mean(d_score) iqr(d_score)]

% Dyad vs solo hit rate comparison
cnt = 0;
for iSubj = 1:length(dyad_perf)
    if ~isempty(dyad_perf{iSubj})
        cnt      	= cnt+1;
        idx        	= cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),solo_id);
        d_hir(cnt)     = mean(dyad_perf{iSubj}.hir);
        s_hir(cnt)     = mean(solo_perf{idx}.hir);
        
        d_better(cnt) 	= (d_hir(cnt) - s_hir(cnt)) > 0;
    end
end
n                   = sum(d_better);
rate                = n / length(d_better);
[pv,h,stats]     	= signrank(d_hir,s_hir); % paired, two-sided test for the null hypothesis that x – y comes from a distribution with zero median
[mean(s_hir) iqr(s_hir)]
[mean(d_hir) iqr(d_hir)]

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

% Test across subjects
median(auc_ecc_pooled)
[p,h,stats] = signrank(auc_ecc_pooled,.5)
median(auc_acc_pooled)
[p,h,stats] = signrank(auc_acc_pooled,.5)

%%% Lag difference
snr = unique(dyad_cr{1}.coh);
cnt = 0;
for iSubj = 1:length(dyad_cr)
    if ~isempty(dyad_perf{iSubj})
        cnt                         = cnt+1;
        for iCoh = 1:length(snr)
            avg_dyad_lag(cnt,iCoh)	= mean(dyad_cr{iSubj}.lag(dyad_cr{iSubj}.coh == snr(iCoh)));
%             iqr_dyad_lag(cnt,iCoh)	= iqr(dyad_cr{iSubj}.lag(dyad_cr{iSubj}.coh == snr(iCoh)));
            idx                   	= cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),solo_id);
            avg_solo_lag(cnt,iCoh)	= mean(solo_cr{idx}.lag(solo_cr{idx}.coh == snr(iCoh)));
%             iqr_solo_lag(cnt,iCoh)	= iqr(solo_cr{idx}.lag(solo_cr{idx}.coh == snr(iCoh)));
        end
    end
end

[mean(mean(avg_solo_lag,2)) iqr(mean(avg_solo_lag,2))]
[mean(mean(avg_dyad_lag,2)) iqr(mean(avg_dyad_lag,2))]

[p_lag,h_lag,stats_lag] = signrank(mean(avg_solo_lag,2), mean(avg_dyad_lag,2));

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

function [ax] = plotScatter(ax, solo_dat, dyad_dat, snr, lb_fs, ax_flag, auc_flag)


hold on
coh_col                     = cool(length(snr));

for iCoh = 1:size(solo_dat,2)
    xx                      = solo_dat(:,iCoh);
    yy                      = dyad_dat(:,iCoh);
    
    sc                      = scatter(xx,yy);
    sc.MarkerFaceColor      = coh_col(iCoh,:)/2;
    sc.MarkerEdgeColor      = 'none';
    sc.MarkerFaceAlpha      = .3;
    sc.SizeData             = 20;
end

for iCoh = 1:size(solo_dat,2)
    xx                      = solo_dat(:,iCoh);
    yy                      = dyad_dat(:,iCoh);
    x_ci                    = (bootci(500, {@nanmedian,  xx},'alpha', .01));
    y_ci                    = (bootci(500, {@nanmedian,  yy},'alpha', .01));
    lny                     = line([nanmedian(xx) nanmedian(xx)],[y_ci(1) y_ci(2)], 'Color', coh_col(iCoh,:),'LineWidth',1);
    lnx                     = line([x_ci(1) x_ci(2)],[nanmedian(yy) nanmedian(yy)], 'Color', coh_col(iCoh,:),'LineWidth',1);
    sc                      = scatter(nanmedian(xx),nanmedian(yy), 'MarkerFaceColor', coh_col(iCoh,:),'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1,'SizeData', 20);
end

ax.FontSize                 = lb_fs;
ax.XLabel.String            = 'Solo [norm]';
ax.XLim                 	= [0 1];
ax.XTick                    = [.1 .5 .9];

if ax_flag
    ax.XAxis.Visible        = 'off';
end

if auc_flag
    ln                    	= line([0 1],[.5 .5]);
    ax.YLim                 = [.1 .9];
    ax.YTick               	= [.25 .5 .75];
    ax.YLabel.String     	= 'AUC';
else
    ln                   	= line([0 1],[0 0]);
    ax.YTick                = [-.2 0 .2];
    ax.YLim                 = [-.25 .25];
    ax.YLabel.String      	= 'Difference';
end

ln.LineWidth                = 1.5;
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];

end

function [x,y,auc] = getAUROC(in1, in2)

if size(in1,1) == 1
    in1         = in1';
    in2         = in2';
end

lab          	= [zeros(length(in1),1); ones(length(in2),1)];
[x,y,~,auc]     = perfcurve(lab,double([in1; in2]),1);
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

function [ax,mdl,R,PV] = plotQuartiles(ax,x,y,str,lw,lb_fs,auc_flag)

addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
coh_col                 = cool(size(y,2));
quartile_boundaries   	= prctile(x, [25 50 75]);
out                     = nan([size(y,2),4]);

for iCoh = 1:size(y,2)
    
    mdl{iCoh}           = fitlm(x(:,iCoh),y(:,iCoh));
    [R{iCoh},PV{iCoh}]  = corrcoef(x(:,iCoh),y(:,iCoh));
    
    for iQuart = 1:4
        if iQuart == 1
            idx         = x(:,iCoh) < quartile_boundaries(1,iCoh);
        elseif iQuart == 2 || iQuart == 3
            idx         = x(:,iCoh)  >= quartile_boundaries(iQuart-1,iCoh) & x(:,iCoh)  <= quartile_boundaries(iQuart,iCoh);
        elseif iQuart == 4
            idx         = x(:,iCoh)  > quartile_boundaries(3,iCoh);
        end
        
        dat(iCoh,iQuart) = mean(y(idx,iCoh));
    end
    
    plot(dat(iCoh,:),'Color',coh_col(iCoh,:),'LineWidth', lw/2)
end

ax.XLim                 = [0.75 4.25];
ax.XTick                = [1:4];
ax.XTickLabel           = {'Bottom','2nd','3rd','Top'};
ax.XLabel.String      	= 'Solo performance [Quartiles]';
ax.FontSize             = lb_fs;
ax.XTickLabelRotation   = 0;

if auc_flag
    ln                      = line([0 5],[.5 .5]);
    ax.YLabel.String     	= str;
else
    ln                   	= line([0 5],[0 0]);
    ax.YLabel.String      	= str;
end

ln.LineWidth            = 1.5;
ln.LineStyle            = ':';
ln.Color                = [0 0 0];
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