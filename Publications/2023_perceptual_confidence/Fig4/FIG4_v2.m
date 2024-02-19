% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hc_dyad_pairwise_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hc_dyad_pairwise_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_performance.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Within-dyad effect size %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[acc_df, ecc_df, auc, score, hir, ply1_flag, ply2_flag, raw] = ...
    dyad_effect_size(solo_perf, dyad_pw_perf, dyad_perf);


%%% PLOT PARAMS %%%
dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG4/';
rcol                       	= [.6 .1 .1];
scol                        = [.1 .1 .1];
hcol                        = [.1 .1 .1];
lw                         	= 1;
lb_fs                    	= 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Scatter dyadic modulation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                           = figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
ofs                         = .025;
label_neg                   = false;
win_size                    = 12;
nBin                        = 40;
h_ofs                       = .3;

for iPlot = 1:4
    
    ax = subplot(2,2,iPlot); hold on
    
    ln                     	= line([0 0],[0 1]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    if iPlot < 3
        ln                     	= line([-.35 .35],[.5 .5]);
        ln.Color              	= [0 0 0];
        ln.LineWidth           	= lw;
        ln.LineStyle           	= ':';
    end
    
    % Raw dyadic AUC data
    if iPlot == 1
        plot_scatter(ecc_df.solo,auc.ecc1,auc.ecc2,label_neg,ply1_flag,ply2_flag)
        
        ax0v                    = axes('Position', [ax.Position(1)+h_ofs ax.Position(2) ax.Position(3)/5 ax.Position(4)]); hold on
        [v, edg]                = histcounts([auc.ecc1,auc.ecc2],nBin);
        cntr                    = edg(1:end-1) + diff(edg) ./ 2;
        st                      = stairs(v,cntr);
        st.LineWidth            = lw;
        st.Color                = [scol];
        ln                      = line([0 max(v)],[median([auc.ecc1,auc.ecc2]) median([auc.ecc1 auc.ecc2])]);
        ln.LineWidth            = lw;
        ln.Color                = rcol;
        ax0v.XAxis.Visible      = 'off';
        ax0v.YAxis.Visible      = 'off';
        
        ax.YLabel.String        = 'Eccentricity [AUC: Solo vs Dyadic]';
        str                     = 'Eccentricity';
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [.1 .9];
        ax0v.YLim               = [.1 .9];
        
    elseif iPlot == 2
        plot_scatter(acc_df.solo,auc.acc1,auc.acc2,label_neg,ply1_flag,ply2_flag)
        
        ax0v                    = axes('Position', [ax.Position(1)+h_ofs ax.Position(2) ax.Position(3)/5 ax.Position(4)]); hold on
        [v, edg]                = histcounts([auc.acc1,auc.acc2],nBin);
        cntr                    = edg(1:end-1) + diff(edg) ./ 2;
        st                      = stairs(v,cntr);
        st.LineWidth            = lw;
        st.Color                = [scol];
        ln                      = line([0 max(v)],[median([auc.acc1,auc.acc2]) median([auc.acc1,auc.acc2])]);
        ln.LineWidth            = lw;
        ln.Color                = rcol;    
        ax0v.XAxis.Visible      = 'off';
        ax0v.YAxis.Visible      = 'off';
        
        ax.YLabel.String        = 'Accuracy [AUC: Solo vs Dyadic]';
        str                     = 'Accuracy';
        ax.XLim                 = [-.1 .1];
        ax.YLim                 = [.3 .7];
        ax0v.YLim               = [.3 .7];

    % AUC distance between players        
    elseif iPlot == 3
        x                       = ecc_df.solo';
        y                       = abs(auc.ecc1-auc.ecc2)';
        sc                      = scatter(x,y,'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc.SizeData             = 10;
        [ft1,gof1]              = fit(x,y,'poly1');
        [ft2,gof2]              = fit(x,y,'poly2');
        ecc_rsquared            = [gof1.adjrsquare gof2.adjrsquare];
        x_vec                   = -.35:.01:.35;
        p_fit                   = plot(x_vec,ft2(x_vec), 'Color', rcol, 'LineWidth',lw);
        tx                      = text(.001, .4,['adjR^{2}: ' num2str(round(gof2.adjrsquare,3))]);
        tx.FontSize             = lb_fs;
        ax.YLabel.String        = 'AUC difference';
        str                     = 'Eccentricity';
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [0 .5];
        [x_avg, y_avg] = avg_win_value(x,y, win_size);
    elseif iPlot == 4
        x                       = acc_df.solo';
        y                       = abs(auc.acc1-auc.acc2)';
        sc                      = scatter(x,y,'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc.SizeData             = 10;
        [ft1,gof1]              = fit(x,y,'poly1');
        [ft2,gof2]              = fit(x,y,'poly2');
        acc_rsquared            = [gof1.adjrsquare gof2.adjrsquare];
        x_vec                   = -.1:.01:.1;
        p_fit                   = plot(x_vec,ft2(x_vec), 'Color', rcol, 'LineWidth',lw);
        tx                      = text(.001, .175,['adjR^{2}: ' num2str(round(gof2.adjrsquare,3))]);
        tx.FontSize             = lb_fs;
        ax.YLabel.String        = 'AUC difference';
        str                     = 'Accuracy';
        ax.XLim                 = [-.1 .1];
        ax.YLim                 = [0 .2];  
        [x_avg, y_avg]          = avg_win_value(x,y,win_size);
    end
    
    ax.Position(1)              = ax.Position(1) - ofs;
    ax.FontSize                 = lb_fs;
    ax.XLabel.String            = {'Within-dyad difference [P1 - P2]',str};
end

if label_neg
    print(f, [dest_dir '/FIG4_auc_labelled'], '-r500', '-dpng');
    print(f, [dest_dir '/FIG4_auc_labelled'], '-r500', '-dsvg', '-painters');
else
    print(f, [dest_dir '/FIG4_auc'], '-r500', '-dpng');
    print(f, [dest_dir '/FIG4_auc'], '-r500', '-dsvg', '-painters');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Dyadic modulation - histogram %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG4/';
f                           = figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
ofs                         = .025;
nbin                        = 12;

for iPlot = 1:4
    ax                      = subplot(2,2,iPlot); hold on
    
    if iPlot == 1
        auc1                    = auc.ecc1 - .5;
        auc2                    = auc.ecc2 - .5;
        h                       = histogram([auc1 auc2],nbin,'FaceColor',hcol, 'EdgeColor','none');
        h.BinLimits             = [-.5 .5];
        ax.XLabel.String        = 'Eccentricity change';
        ax.YLabel.String        = '# Players';
        ax.XLim                 = h.BinLimits;
        ax.YLim                 = [0 max(h.BinCounts)];

        ln                     	= line([median([auc1 auc2]) median([auc1 auc2])],[0 max(h.BinCounts)]);
        ln.Color              	= rcol;
        ln.LineWidth           	= lw;
        ln.LineStyle           	= ':';
    
    elseif iPlot == 2
        auc1                    = auc.acc1 - .5;
        auc2                    = auc.acc2 - .5;
        h                       = histogram([auc1 auc2],nbin,'FaceColor',hcol, 'EdgeColor','none');
        h.BinLimits             = [-.2 .2];
        ax.XLabel.String        = 'Accuracy change';
        ax.YLabel.String        = '# Players';
        ax.XLim                 = h.BinLimits;
        ax.YLim                 = [0 max(h.BinCounts)];
        %%%  Significance test             signrank([[auc1 auc2]])   
        ln                     	= line([median([auc1 auc2]) median([auc1 auc2])],[0 max(h.BinCounts)]);
        ln.Color              	= rcol;
        ln.LineWidth           	= lw;
        ln.LineStyle           	= ':';

    elseif iPlot == 3
        
        ln                     	= line([0 1],[.5 .5]);
        ln.Color              	= [0 0 0];
        ln.LineWidth           	= lw;
        ln.LineStyle           	= ':';
        
        ln                     	= line([.5 .5],[0 1]);
        ln.Color              	= [0 0 0];
        ln.LineWidth           	= lw;
        ln.LineStyle           	= ':';
        
        sc                      = scatter([auc.ecc1 auc.ecc2], [auc.acc1 auc.acc2]);
        sc.MarkerFaceColor      = [hcol];
        sc.MarkerEdgeColor      = 'none';
        sc.MarkerFaceAlpha      = .5;
        ax.XLabel.String        = 'Eccentricity change [AUC]';
        ax.YLabel.String        = 'Accuracy change [AUC]';
        ax.XLim                 = [0 1];
        ax.YLim                 = [0 1];
        
        [r,pv]               	= corrcoef([auc.ecc1 auc.ecc2],[auc.acc1 auc.acc2]);
        r                       = r(2);
        pv                      = pv(2);
        
        linearCoefficients      = polyfit([auc.ecc1 auc.ecc2], [auc.acc1 auc.acc2], 1);
        x_fit                   = linspace(-1, 1, 50);
        y_fit                   = polyval(linearCoefficients, x_fit);
        pl                      = plot(x_fit, y_fit, '-', 'LineWidth', 1);

        if pv < .05
            pl.Color         	= rcol;
            tx                 	= text(.4,.15,{['r=' num2str(round(r,2))];['p=' num2str(round(pv,2))]});
            tx.Color          	= hcol;
            tx.FontSize      	= 8;
        end
    end
        
    ax.Position(1)              = ax.Position(1) - ofs;
    ax.FontSize                 = lb_fs;
end

if label_neg
    print(f, [dest_dir '/FIG4_auc_labelled_hist'], '-r500', '-dpng');
    print(f, [dest_dir '/FIG4_auc_labelled_hist'], '-r500', '-dsvg', '-painters');
else
    print(f, [dest_dir '/FIG4_auc_hist'], '-r500', '-dpng');
    print(f, [dest_dir '/FIG4_auc_hist'], '-r500', '-dsvg', '-painters');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Effect size difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                          	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on

for i = 1:2
    ax = subplot(2,2,i); hold on
    
    ln                     	= line([0 0],[-1 1]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    ln                     	= line([-1 1],[0 0]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    if i == 1
        df_subset               = auc.ecc1 - auc.ecc2;
        sc1                     = scatter(ecc_df.solo,df_subset,'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc1.SizeData            = 10;
        
        [x_fit, y_fit, r(i,1), pv(i,1)]	= regr_line(ax,ecc_df.solo,df_subset,rcol, .3);
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [-.55 .55];
        ax.YTick                = [-.5:.25:.5];
        ax.YLabel.String        = {'Eccentricity Difference';' [AUROC_P1 - AUROC_P2]'};
    elseif i == 2
        df_subset               = auc.acc1 - auc.acc2;
        sc2                     = scatter(acc_df.solo,df_subset,'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc2.SizeData            = 10;
        
        [x_fit, y_fit, r(i,2), pv(i,2)]	= regr_line(ax,acc_df.solo,df_subset,rcol,.13);
        ax.XLim                 = [-.1 .1];
        ax.YLim                 = [-.2 .2];
        ax.YTick                = [-.2:.1:.2];
        ax.YLabel.String        = {'Accuracy Difference';' [AUROC_P1 - AUROC_P2]'};
    end
    
    ax.YLabel.Interpreter       = 'none';
    ax.Position(1)              = ax.Position(1);
    ax.FontSize                 = lb_fs;
    ax.XLabel.String            = {'Within-dyad difference [P1 - P2]'};
    
end
print(f, [dest_dir '/FIG4_dff_scatter'], '-r500', '-dpng');
print(f, [dest_dir '/FIG4_dff_scatter'], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Solo vs dyadic difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                      	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on

for iPlot = 1:2
    ax                      = subplot(2,2,iPlot); hold on

    if iPlot == 1
        sc1                     = scatter(abs(ecc_df.solo),abs(ecc_df.dyad),'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc1.MarkerFaceAlpha     = .5;
        sc1.SizeData            = 10;
        ax.XLim                 = [0 .4];
        ax.YLim                 = [0 .4];
        subset_str              = 'Eccentricity';
        
        [r_ecc,pv_ecc]          = corrcoef(abs(ecc_df.solo),abs(ecc_df.dyad));
        [p_ecc,h_ecc, stats_ecc] = signrank(abs(ecc_df.solo),abs(ecc_df.dyad));
        
        if pv_ecc(2) < .05
            tx               	= text(0,.3,{['r=' num2str(round(r_ecc(2),2))];['p=' num2str(round(pv_ecc(2),2))]});
            tx.Color         	= [.1 .1 .1];
            tx.FontSize        	= lb_fs;
        end

    elseif iPlot == 2
        sc2                     = scatter(abs(acc_df.solo),abs(acc_df.dyad),'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc2.MarkerFaceAlpha     = .5;
        sc2.SizeData            = 10;
        ax.XLim                 = [0 .12];
        ax.YLim                 = [0 .12];
        subset_str              = 'Accuracy';
                
        [r_acc,pv_acc]          = corrcoef(abs(acc_df.solo),abs(acc_df.dyad));
        [p_acc,h_acc, stats_acc] = signrank(abs(acc_df.solo),abs(acc_df.dyad));

        if pv_acc(2) < .05
            tx                 	= text(0,.3,{['r=' num2str(round(r_acc(2),2))];['p=' num2str(round(pv_acc(2),2))]});
            tx.Color         	= [.1 .1 .1];
            tx.FontSize     	= lb_fs;
        end
    end

    lsl                         = lsline;
    lsl.Color                   = rcol;
    lsl.LineWidth               = lw;

    ln                          = line([0 .4],[0 .4]);
    ln.Color                    = scol;
    ln.LineWidth                = lw;
    ln.LineStyle                = ':';
    
    ax.FontSize                 = lb_fs;
    ax.XLabel.String            = {'Difference [Solo]'; subset_str};
    ax.YLabel.String            = 'Difference [Dyad]';
    ax.XTickLabelRotation       = 0;
    ax.YLabel.Interpreter       = 'none';
    ax.Position(1)              = ax.Position(1);
    ax.FontSize                 = lb_fs;
end

print(f, [dest_dir '/FIG4_dff_solo_dyad'], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG4_dff_solo_dyad'], '-r500', '-dpng', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Solo vs dyadic correlation between players
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f                      	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on

for i = 1:4
    
    if i == 1
        dat1            = raw.acc1;
        dat2            = raw.acc2;
        str             = 'Solo Accuracy';
    elseif i == 2
        dat1            = raw.ecc1;
        dat2            = raw.ecc2;
        str             = 'Solo Eccentricity';
    elseif i == 3
        dat1            = raw.dacc1;
        dat2            = raw.dacc2;
        str             = 'Dyadic Accuracy';
    elseif i == 4
        dat1            = raw.decc1;
        dat2            = raw.decc2;
        str             = 'Dyadic Eccentricity';
    end
    
    ax                  = subplot(2,2,i);
    
    sc                  = scatter(dat1,dat2);
    sc.MarkerFaceColor  = scol;
    sc.MarkerFaceAlpha  = .5;
    sc.MarkerEdgeColor  = 'none';
    sc.SizeData         = 20; 
    
    if i == 2 || i == 4
        ax.XLim     	= [.4 1];
        ax.YLim      	= [.4 1];
        ax.XTick        = [.4:.2:1];
        ax.YTick        = [.4:.2:1];
    else
        ax.XLim        	= [.7 1];
        ax.YLim       	= [.7 1];
        ax.XTick        = [.7:.1:1];
        ax.YTick        = [.7:.1:1];
    end
    
    ax.FontSize         = lb_fs;
    ax.XLabel.String    = 'Player1';
    ax.YLabel.String    = 'Player2';
    
    lsl                 = lsline;
    lsl.Color          	= rcol;
    lsl.LineWidth     	= lw;
    
    [r,pv]              = corrcoef(dat1,dat2)
    title({str, [' r = ' num2str(round(r(2),2)) ' pv = ' num2str(round(pv(2),2))]})
end

print(f, [dest_dir '/FIG4_corr_solo_dyad'], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG4_corr_solo_dyad'], '-r500', '-dpng', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Solo difference vs average social modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f                      	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on

for i = 1:4
    
    if i == 1
        dat1            = abs(ecc_df.solo);
        dat2            = mean([auc.ecc1; auc.ecc2]);
%         dat2            = mean(abs([auc.ecc1; auc.ecc2] - .5));
        str             = 'Eccentricity difference';
    elseif i == 2
        dat1            = abs(acc_df.solo);
        dat2            = mean([auc.ecc1; auc.ecc2]);
%         dat2            = mean(abs([auc.ecc1; auc.ecc2] - .5));
        str             = 'Accuracy difference';
    elseif i == 3
        dat1            = abs(ecc_df.solo);
        dat2            = mean([auc.acc1; auc.acc2]);
%         dat2            = mean(abs([auc.acc1; auc.acc2] - .5));
        str             = 'Eccentricity difference';
    elseif i == 4
        dat1            = abs(acc_df.solo);
        dat2            = mean([auc.acc1; auc.acc2]);
%         dat2            = mean(abs([auc.acc1; auc.acc2] - .5));
        str             = 'Accuracy difference';
    end
    
    ax                  = subplot(2,2,i);
    
    sc                  = scatter(dat1,dat2);
    sc.MarkerFaceColor  = scol;
    sc.MarkerFaceAlpha  = .5;
    sc.MarkerEdgeColor  = 'none';
    sc.SizeData         = 20; 
    
    ax.FontSize         = lb_fs;
    ax.XLabel.String    = str;
    ax.YLabel.String    = 'Avg. AUC';
    
    lsl                 = lsline;
    lsl.Color          	= rcol;
    lsl.LineWidth     	= lw;
    
    ln                     	= line([0 ax.XLim(2)],[.5 .5]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    [r,pv]              = corrcoef(dat1,dat2);
    title({[' r = ' num2str(round(r(2),2)) ' pv = ' num2str(round(pv(2),2))]})
end


print(f, [dest_dir '/FIG4_corr_solo_auc'], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG4_corr_solo_auc'], '-r500', '-dpng', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Solo difference vs dyadic performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f                      	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on

for i = 1:4
    
    if i == 1
        dat1            = abs(ecc_df.solo);
%         dat1            = abs(ecc_df.dyad);
        dat2            = mean(score.dyad,2);
        xstr            = 'Eccentricity difference';
        ystr            = 'Dyadic score';
    elseif i == 2
        dat1            = abs(acc_df.solo);
%         dat1            = abs(acc_df.dyad);
        dat2            = mean(score.dyad,2);
        xstr             = 'Accuracy difference';
        ystr            = 'Dyadic score';
    elseif i == 3
        dat1            = abs(ecc_df.solo);
%         dat1            = abs(ecc_df.dyad);
        dat2            = mean(hir.dyad,2);
        xstr            = 'Eccentricity difference';
        ystr            = 'Dyadic hit rate';
        
    elseif i == 4
        dat1            = abs(acc_df.solo);
%         dat1            = abs(acc_df.dyad);
        dat2            = mean(hir.dyad,2);
        xstr            = 'Accuracy difference';
        ystr            = 'Dyadic hit rate';
    end
    
    ax                  = subplot(2,2,i);
    
    sc                  = scatter(dat1,dat2);
    sc.MarkerFaceColor  = scol;
    sc.MarkerFaceAlpha  = .5;
    sc.MarkerEdgeColor  = 'none';
    sc.SizeData         = 20; 
    
    ax.FontSize         = lb_fs;
    ax.XLabel.String    = xstr;
    ax.YLabel.String    = ystr;
    
    lsl                 = lsline;
    lsl.Color          	= rcol;
    lsl.LineWidth     	= lw;
   
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    [r,pv]              = corrcoef(dat1,dat2);
    title({[' r = ' num2str(round(r(2),2)) ' pv = ' num2str(round(pv(2),2))]})
end

print(f, [dest_dir '/FIG4_corr_solo_score'], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG4_corr_solo_score'], '-r500', '-dpng', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reported stats in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Different from 0.5
median([auc.ecc1,auc.ecc2])
[p,h,stats] = signrank([auc.ecc1,auc.ecc2],.5);

median([auc.acc1,auc.acc2])
[p,h,stats] = signrank([auc.acc1,auc.acc2],.5);

% Regression to mean problem
nRep                        = 5000;
auc_str                     = 'acc';
solo_str                    = 'ecc';
n                           = 1;
[shuffl_coeff, true_coeff]  = regr_ci(nRep,auc,raw,auc_str,solo_str);
p_actual                    = mean(shuffl_coeff(:,1) <= true_coeff(1));

%%% See coefficients here
% figure
% histogram(shuffl_coeff(:,n)) 
% line([true_coeff(n) true_coeff(n)],[0 1000],'Color',[0 0 0], 'LineWidth', 2, 'LineStyle','--')
% ylim([0 500])

% Dyadic improvement
ecc_class       = sum([auc.ecc1' auc.ecc2']>.5,2);
ecc_both_better = sum(ecc_class == 2) / length(ecc_class);
ecc_both_worse  = sum(ecc_class == 0) / length(ecc_class);
ecc_mixed       = sum(ecc_class == 1) / length(ecc_class);

acc_class       = sum([auc.acc1' auc.acc2']>.5,2);
acc_both_better = sum(acc_class == 2) / length(acc_class);
acc_both_worse  = sum(acc_class == 0) / length(acc_class);
acc_mixed       = sum(acc_class == 1) / length(acc_class);

% Conversion
sum(abs(ecc_df.solo) > abs(ecc_df.dyad)) / length(ecc_df.dyad)
sum(abs(acc_df.solo) > abs(acc_df.dyad)) / length(acc_df.dyad)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Covert to table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = 'ecc';
df = ecc_df;
dat_unsorted = [df.solo',df.dyad',...
    raw.([str '1'])', raw.(['d' str '1'])',auc.([str '1'])',...
    raw.([str '2'])',raw.(['d' str '2'])',auc.([str '2'])'];

[~, idx] = sort(df.solo);

dat = dat_unsorted(idx,:);

var_names = {'Solo Difference','Dyadic Difference',...
    'Solo Player1','Dyadic Player1','AUC Player1',...
    'Solo Player2','Dyadic Player2','AUC Player2'};

for i = 1:length(acc_df.solo)
    row_names{i} = ['Dyad' num2str(i)];
end

summary_table = array2table(dat, 'VariableNames',var_names,'RowNames',row_names);
writetable(summary_table, 'summary_table.xlsx')

% T = summary_table;
% % Get the table in string form.
% TString = evalc('disp(T)');
% % Use TeX Markup for bold formatting and underscores.
% TString = strrep(TString,'<strong>','\bf');
% TString = strrep(TString,'</strong>','\rm');
% TString = strrep(TString,'_','\_');
% % Get a fixed-width font.
% FixedWidth = get(0,'FixedWidthFontName');
% % Output the table using the annotation command.
% annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
%     'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [acc_df, ecc_df, auc, score, hir, ply1_flag, ply2_flag,raw] = dyad_effect_size(in_solo, in_dyad, dyad_pooled)

for iSubj = 1:size(in_solo,2)
    solo_id{iSubj}          = in_solo{iSubj}.id;
end

for iDyad = 1:size(in_dyad,1)
    % Find solo data of subjects
    idx_ply1                = cellfun(@(x) strcmp(x,in_dyad{iDyad,1}.id),solo_id);
    idx_ply2                = cellfun(@(x) strcmp(x,in_dyad{iDyad,2}.id),solo_id);
    
    score.solo(iDyad,:)   	= [mean(in_solo{idx_ply1}.trg_all_score) mean(in_solo{idx_ply2}.trg_all_score)];
    score.dyad(iDyad,:)  	= [mean(in_dyad{iDyad,1}.trg_all_score) mean(in_dyad{iDyad,2}.trg_all_score)];
    
    hir.solo(iDyad,:)    	= [sum(in_solo{idx_ply1}.trg_all_outc)/length(in_solo{idx_ply1}.trg_all_outc) sum(in_solo{idx_ply2}.trg_all_outc)/length(in_solo{idx_ply2}.trg_all_outc)];
    hir.dyad(iDyad,:)     	= [sum(in_dyad{iDyad,1}.trg_all_outc)/length(in_dyad{iDyad,1}.trg_all_outc) sum(in_dyad{iDyad,2}.trg_all_outc)/length(in_dyad{iDyad,2}.trg_all_outc)];

    % Label subjects of interest
    % soi                     = ID_subj_negAUC(in_solo, dyad_pooled);
    % ply1_flag(iDyad)      	= logical(sum(cellfun(@(x) strcmp(x,in_dyad{iDyad,1}.id),soi)));
    % ply2_flag(iDyad)      	= logical(sum(cellfun(@(x) strcmp(x,in_dyad{iDyad,2}.id),soi)));
    ply1_flag(iDyad) = false;
    ply2_flag(iDyad) = false;
    
    % Raw performance
    raw.acc1(iDyad)         = mean(in_solo{idx_ply1}.macc_trg); % Mean(Median per coherence)
    raw.acc2(iDyad)         = mean(in_solo{idx_ply2}.macc_trg);
    raw.ecc1(iDyad)         = mean(in_solo{idx_ply1}.mecc_state);
    raw.ecc2(iDyad)         = mean(in_solo{idx_ply2}.mecc_state);
        
    raw.dacc1(iDyad)        = mean(in_dyad{iDyad,1}.macc_trg);
    raw.dacc2(iDyad)        = mean(in_dyad{iDyad,2}.macc_trg);
    raw.decc1(iDyad)        = mean(in_dyad{iDyad,1}.mecc_state);
    raw.decc2(iDyad)        = mean(in_dyad{iDyad,2}.mecc_state);
    
    % Performance difference
    acc_df.solo(iDyad)      = mean(in_solo{idx_ply1}.macc_trg) - mean(in_solo{idx_ply2}.macc_trg);
    ecc_df.solo(iDyad)     	= mean(in_solo{idx_ply1}.mecc_state) - mean(in_solo{idx_ply2}.mecc_state);
    acc_df.dyad(iDyad)      = mean(in_dyad{iDyad,1}.macc_trg) - mean(in_dyad{iDyad,2}.macc_trg);
    ecc_df.dyad(iDyad)     	= mean(in_dyad{iDyad,1}.mecc_state) - mean(in_dyad{iDyad,2}.mecc_state);
       
    % Effect size: Solo vs Dyadic
    [auc.ecc1(iDyad), auc.ecc2(iDyad)] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, 'ecc_state');
    [auc.acc1(iDyad), auc.acc2(iDyad)] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, 'acc_trg');
end
end

function [out_ply1, out_ply2] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, str)

for i = 1:length(in_solo{idx_ply1}.(str))
    out1(i)     = getAUROC(in_solo{idx_ply1}.(str){i},in_dyad{iDyad,1}.(str){i});
    out2(i)     = getAUROC(in_solo{idx_ply2}.(str){i},in_dyad{iDyad,2}.(str){i}); 
end

out_ply1 = mean(out1);
out_ply2 = mean(out2);

end

function [out] = getAUROC(in1, in2)

if size(in1,1) == 1
    in1             = in1';
    in2             = in2';
end

lab                 = [zeros(length(in1),1); ones(length(in2),1)];
[~,~,~,out]         = perfcurve(lab,[in1; in2],1);

% figure
% histogram(in1,20)
% hold on
% histogram(in2,20)
% title([median(in1) median(in2) out])

end

function out = ID_subj_negAUC(solo_perf, dyad_perf)

for j = 1:size(dyad_perf,2)
    if isempty(dyad_perf{j})
        ply_id{j} = nan;
    else
        ply_id{j} = dyad_perf{j}.id;
    end
end

cnt = 0;
for iSubj = 1:size(solo_perf,2)
    
    d.sp                  	= solo_perf{iSubj}; % Performance data
    sIdx                    = cellfun(@(x) strcmp(x,solo_perf{iSubj}.id),ply_id);
    
    if sum(sIdx) == 0
        continue
    end
    
    d.dp                  	= dyad_perf{sIdx};
    
    if isempty(d.dp)
        continue
    end
    
    cnt = cnt+1;
    id_auc{cnt}             = dyad_perf{sIdx}.id;
    
    for iCoh = 1:length(d.sp.carr)
        auc_ecc(cnt,iCoh)	= getAUROC(d.sp.ecc_state{iCoh},d.dp.ecc_state{iCoh});
    end
end

% Identify subjects
out                         = id_auc(mean(auc_ecc,2)<.5);

end

function plot_scatter(df,auc1,auc2,label_neg,ply1_flag,ply2_flag)

if ~label_neg
    ply1_flag = false([1 length(df)]); % Overwrite flags
    ply2_flag = false([1 length(df)]);
end

col                     = [.3 .3 .3];

for i = 1:length(df)
    ln                  = line([df(i) df(i)],[auc1(i) auc2(i)]);
    ln.Color            = col;
    ln.LineWidth        = .5;
    
    if ply1_flag(i)
        sc1 = scatter(df(i),auc1(i),'MarkerFaceColor',[1 .3 .3],'MarkerFaceAlpha',.5, 'MarkerEdgeColor',col,'MarkerEdgeAlpha',1,'Marker','o');
    else
        sc1 = scatter(df(i),auc1(i),'MarkerFaceColor',col,'MarkerFaceAlpha',.5, 'MarkerEdgeColor',col,'MarkerEdgeAlpha',1,'Marker','o');
    end
    
    if ply2_flag(i)
        sc2 = scatter(df(i),auc2(i),'MarkerFaceColor','none','MarkerEdgeColor',[1 .3 .3],'MarkerEdgeAlpha',1,'Marker','o');
    else
        sc2 = scatter(df(i),auc2(i),'MarkerFaceColor','none','MarkerEdgeColor',col,'MarkerEdgeAlpha',1,'Marker','o');
    end
    
    sc1.SizeData                = 10;
    sc2.SizeData                = 10;
end
end

function [x_fit, y_fit, r, pv] = regr_line(ax,x_in,y_in,col,height)

if nargin < 5
    height = 0;
end

axes(ax); hold on

[r,pv]= corrcoef(x_in,y_in);
r = r(2);
pv = pv(2);

linearCoefficients = polyfit(x_in, y_in, 1);
x_fit = linspace(-.5, .5, 50);
y_fit = polyval(linearCoefficients, x_fit);
pl = plot(x_fit, y_fit, '-', 'LineWidth', 1,'Color',col);

if pv < .05
    pl.Color                    = col;
    tx                          = text(0,height,{['r=' num2str(round(r,2))];['p=' num2str(round(pv,2))]});
    tx.Color                    = [.1 .1 .1];
    tx.FontSize                 = 8;
end

end

function [f, pv, r] = scatter_perf_dff(subset_str,ply1_flag,ply2_flag,auc,acc_df,ecc_df)

f                               = figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
marker                          = {'x','+'};
col                             = {[230,97,1]./255;[94,60,153]./255};
lw                              = 1;
lb_fs                           = 8;

if strcmp(subset_str, 'pos')
    idx                   	= ~ply1_flag & ~ply2_flag; % pos
elseif strcmp(subset_str, 'neg')
    idx                     = ply1_flag | ply2_flag; % neg
else
    idx                     = true([length(acc_df) 1]);
end

auc_ecc1_subset             = auc.ecc1(idx);
auc_ecc2_subset             = auc.ecc2(idx);
auc_acc1_subset             = auc.acc1(idx);
auc_acc2_subset             = auc.acc2(idx);
acc_df_subset               = acc_df(idx);
ecc_df_subset               = ecc_df(idx);

for i = 1:2
    ax = subplot(2,2,i); hold on
    
    ln                     	= line([0 0],[-1 1]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    ln                     	= line([-1 1],[0 0]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    if i == 1
        df_subset               = auc_ecc1_subset - auc_ecc2_subset;
        sc1                     = scatter(ecc_df_subset,df_subset,'LineWidth',lw,'MarkerEdgeColor',col{1},'Marker',marker{1});
        sc2                     = scatter(acc_df_subset,df_subset,'LineWidth',lw,'MarkerEdgeColor',col{2},'Marker',marker{2});
        [x_fit, y_fit, r(i,1), pv(i,1)]	= regr_line(ax,ecc_df_subset,df_subset,col{1});
        [x_fit, y_fit, r(i,2), pv(i,2)]	= regr_line(ax,acc_df_subset,df_subset,col{2},.13);
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [-.55 .55];
        ax.YTick                = [-.5:.25:.5];
        ax.YLabel.String        = {'Eccentricity Difference';' [AUROC_P1 - AUROC_P2]'};
    elseif i == 2
        df_subset               = auc_acc1_subset - auc_acc2_subset;
        sc1                     = scatter(ecc_df_subset,df_subset,'LineWidth',lw,'MarkerEdgeColor',col{1},'Marker',marker{1});
        sc2                     = scatter(acc_df_subset,df_subset,'LineWidth',lw,'MarkerEdgeColor',col{2},'Marker',marker{2});
        [x_fit, y_fit, r(i,1), pv(i,1)]	= regr_line(ax,ecc_df_subset,df_subset,col{1});
        [x_fit, y_fit, r(i,2), pv(i,2)]	= regr_line(ax,acc_df_subset,df_subset,col{2},.13);
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [-.25 .25];
        ax.YTick                = [-.2:.1:.2];
        ax.YLabel.String        = {'Accuracy Difference';' [AUROC_P1 - AUROC_P2]'};
    end
    
    ax.YLabel.Interpreter       = 'none';
    ax.Position(1)              = ax.Position(1);
    ax.FontSize                 = lb_fs;
    ax.XLabel.String            = {'Within-dyad difference [P1 - P2]'};
    
    if i == 1
        lg                      = legend([sc1,sc2],{'Eccentricity','Accuracy'});
        lg.FontSize             = lb_fs;
        lg.Location             = 'southoutside';
        lg.Box                  = 'off';
        lg.Position(2)          = lg.Position(2)-.11;
    end
end
end

function [shuffl_coeff, true_coeff] = regr_ci(nRep,auc,raw,auc_str,solo_str)

auc1                    = auc.([auc_str num2str(1)])';
auc2                    = auc.([auc_str num2str(2)])';
solo1                   = raw.([solo_str num2str(1)])';
solo2                   = raw.([solo_str num2str(2)])';

true_coeff              = polyfit(solo1-solo2, auc1-auc2, 1);

for iRep = 1:nRep
    idx                 = randperm(size(auc1,1));
    shuffled_auc2       = auc2(idx);
    shuffled_solo2      = solo2(idx);
    auc_df              = auc1 - shuffled_auc2;
    solo_df            	= solo1 - shuffled_solo2;
    
    excl                = solo_df == 0;
    auc_df              = auc_df(~excl);
    solo_df             = solo_df(~excl);
    shuffl_coeff(iRep,:)= polyfit(solo_df, auc_df, 1);
end
end

function [x_avg, y_avg] = avg_win_value(x,y,win_size)

[x_val,idx]     = sort(x);
y_val           = y(idx);

x_avg           = movmean(x_val,win_size);
y_avg           = movmean(y_val,win_size);

p_avg       	= plot(x_avg,y_avg);
p_avg.Color     = [.1 .1 .1];
p_avg.LineWidth = 1;
p_avg.LineStyle = '-';

% for i = 1:length(y)-win_size
%     x_avg(i) = mean(x_val(i:i+win_size));
%     y_avg(i) = mean(y_val(i:i+win_size));
% end
end

% function [avgY, binEdges] = avg_bin_value(x,y)
% % Define bin edges
% binEdges = linspace(min(x), max(x), 10);  % Adjust numberOfBins as needed
% 
% % Use histcounts to bin the data
% [counts, binEdges] = histcounts(x, binEdges);
% 
% % Initialize arrays to store the sum and count for each bin
% sumY = zeros(1, length(binEdges) - 1);
% countY = zeros(1, length(binEdges) - 1);
% 
% % Loop through each bin and accumulate the sum and count of Y-values
% for i = 1:length(binEdges) - 1
%     indices = x >= binEdges(i) & x < binEdges(i + 1);
%     sumY(i) = sum(y(indices));
%     countY(i) = sum(indices);
% end
% 
% % Calculate the average Y-value for each bin
% avgY = sumY ./ countY;
% 
% % % Plot the results
% % bar(binEdges(1:end-1), averageY)
% % xlabel('X-axis')
% % ylabel('Average Y-value')
% % title('Binned Data with Average Y-values')
% end
