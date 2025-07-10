close all
clear all

addpath /Users/fschneider/Documents/MATLAB/piermorel-gramm-b0fc592
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

% Adjust path
source_pth = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/var_plot/';
load([source_pth '/solo_correlation.mat'])
load([source_pth '/solo_performance.mat'])
load([source_pth '/hh_dyad_pairwise_correlation.mat'])
load([source_pth 'hh_dyad_pairwise_performance.mat'])
load([source_pth '/hc_dyad_pairwise_correlation.mat'])
load([source_pth '/hc_dyad_pairwise_performance.mat'])
load([source_pth '/hh_dyad_performance.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Within-dyad effect size %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alignment_str = 'state';
% alignment_str = 'trg';

[acc_df, ecc_df, auc, score, hir, raw] = dyad_effect_size(solo_perf, dyad_pw_perf, alignment_str);

%%% PLOT PARAMS %%%
dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/FIG_pairwise_modulation/raw/';
rcol                       	= [.6 .1 .1];
scol                        = [.1 .1 .1];
hcol                        = [.1 .1 .1];
lw                         	= 1;
lb_fs                    	= 8;
sc_size                     = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE: Solo vs dyadic joystick difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear g

%%% Eccentricity
p1_better               = raw.ecc1 > raw.ecc2;
ecc_solo_dff            = [raw.ecc1(p1_better) - raw.ecc2(p1_better) raw.ecc2(~p1_better) - raw.ecc1(~p1_better)];
ecc_dyad_dff            = [raw.decc1(p1_better) - raw.decc2(p1_better) raw.decc2(~p1_better) - raw.decc1(~p1_better)];
[r_ecc,pv_ecc]          = corrcoef(abs(ecc_df.solo),abs(ecc_df.dyad));
[p_ecc,h_ecc, stats_ecc]= signrank(abs(ecc_df.solo),abs(ecc_df.dyad));
g(1,1)                  = gramm('x',abs(ecc_solo_dff)','y',abs(ecc_dyad_dff)');

sum(ecc_solo_dff > ecc_dyad_dff)

g(1,1).geom_point();
g(1,1).stat_cornerhist('edges',-.3:0.015:.3,'aspect',.3,'fill','face');
g(1,1).geom_abline();
g(1,1).set_names('x','Solo interplayer distance','y','Dyadic interplayer distance');
g(1,1).set_color_options('map',[.3 .3 .3]);

%%% Accuracy
p1_better               = raw.acc1 > raw.acc2;
acc_solo_dff            = [raw.acc1(p1_better) - raw.acc2(p1_better) raw.acc2(~p1_better) - raw.acc1(~p1_better)];
acc_dyad_dff            = [raw.dacc1(p1_better) - raw.dacc2(p1_better) raw.dacc2(~p1_better) - raw.dacc1(~p1_better)];
[r_acc,pv_acc]          = corrcoef(abs(acc_df.solo),abs(acc_df.dyad));
[p_acc,h_acc, stats_acc]= signrank(abs(acc_df.solo),abs(acc_df.dyad));
g(1,2)                  = gramm('x',abs(acc_solo_dff),'y',abs(acc_dyad_dff));

sum(acc_solo_dff > acc_dyad_dff)

g(1,2).geom_point();
g(1,2).stat_cornerhist('edges',-.1:0.005:.1,'aspect',.3,'fill','face');
g(1,2).geom_abline();
g(1,2).set_names('x','Solo interplayer distance','y','Dyadic interplayer distance');
g(1,2).set_color_options('map',[.3 .3 .3]);

figure('Position',[100 100 600 300]);
g.draw();

axis normal
axis square

print(gcf, [dest_dir '/FIG_js_dff_' alignment_str], '-r500', '-dsvg', '-painters');
print(gcf, [dest_dir '/FIG_js_dff_' alignment_str], '-r500', '-dpng', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE: P_better minus P-worse AUC difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /Users/fschneider/Documents/MATLAB/cbrewer

n                       = 50;
subsets                 = nan(n,1);
subsets(auc.ecc1 > .5 & auc.ecc2 > .5) = 1;
subsets(auc.ecc1 < .5 & auc.ecc2 < .5) = 2;
subsets(auc.ecc1 < .5 & auc.ecc2 > .5) = 3;
subsets(auc.ecc1 > .5 & auc.ecc2 < .5) = 3;
cval={'both better' 'both worse' 'mixed'};


f                       = figure;
ax                      = subplot(2,2,1);hold on
x                       = 0:.05:.35;
% cmap                  = cbrewer('qual', 'Dark2', 3);
cmap                    = [42 182 115;...
                            190 30 45;...
                            247 147 29]./255;

p1_better               = raw.ecc1 > raw.ecc2; % Solo
ecc_better_minus_worse  = [auc.ecc1(p1_better) - auc.ecc2(p1_better) auc.ecc2(~p1_better) - auc.ecc1(~p1_better)];
ecc_solo_dff            = [raw.ecc1(p1_better) - raw.ecc2(p1_better) raw.ecc2(~p1_better) - raw.ecc1(~p1_better)];
ecc_dyad_dff            = [raw.decc1(p1_better) - raw.decc2(p1_better) raw.decc2(~p1_better) - raw.decc1(~p1_better)];
conv                    = ecc_solo_dff > ecc_dyad_dff;
subsets_remap           = [subsets(p1_better); subsets(~p1_better)];
c                       = cval(subsets_remap);

for i=1:3
    sc                  = scatter(abs(ecc_solo_dff(subsets_remap==i)), ecc_better_minus_worse(subsets_remap==i));
    sc.MarkerFaceColor  = cmap(i,:);
    sc.MarkerEdgeColor  = 'none';
    sc.MarkerFaceAlpha  = .75;
    
    text(.1,i/10,num2str(i), 'Color', cmap(i,:))
end
fit_ecc                 = polyfit(abs(ecc_solo_dff),ecc_better_minus_worse,1);
pl                      = plot(x,polyval(fit_ecc, x),'LineWidth',1.5, 'Color',[.2 .2 .2]);
[r_ecc,pv_ecc]          = corrcoef(abs(ecc_solo_dff),ecc_better_minus_worse);
[r_ecc_conv,pv_ecc_conv] = corrcoef(abs(ecc_solo_dff(conv)),ecc_better_minus_worse(conv));
ax.XLim                 = [0 .35];
ax.YLim                 = [-.55 .55];
ax.YTick                = [-.5 -.25 0  .25 .5];
ax.XLabel.String        = {'Solo interplayer distance'; 'Tilt.better - Tilt.worse'};
ax.YLabel.String        = {'Social modulation difference';'AUC.better - AUC.worse'};
ax.FontSize             = lb_fs;

text(.2,0,['All: r=' num2str(round(r_ecc(2),3)) ' p=' num2str(pv_ecc(2))])
text(.2,.15,['Conv: r=' num2str(round(r_ecc_conv(2),3)) ' p=' num2str(pv_ecc_conv(2))])
axis square

ax                      = subplot(2,2,2);hold on
x                       = 0:.05:.15;
p1_better               = raw.acc1 > raw.acc2; % Solo
acc_better_minus_worse  = [auc.acc1(p1_better) - auc.acc2(p1_better) auc.acc2(~p1_better) - auc.acc1(~p1_better)];
acc_solo_dff            = [raw.acc1(p1_better) - raw.acc2(p1_better) raw.acc2(~p1_better) - raw.acc1(~p1_better)];
acc_dyad_dff            = [raw.dacc1(p1_better) - raw.dacc2(p1_better) raw.dacc2(~p1_better) - raw.dacc1(~p1_better)];
conv                    = acc_solo_dff > acc_dyad_dff;

subsets                 = nan(n,1);
subsets(auc.acc1 > .5 & auc.acc2 > .5) = 1;
subsets(auc.acc1 < .5 & auc.acc2 < .5) = 2;
subsets(auc.acc1 < .5 & auc.acc2 > .5) = 3;
subsets(auc.acc1 > .5 & auc.acc2 < .5) = 3;
subsets_remap           = [subsets(p1_better); subsets(~p1_better)];
c                       = cval(subsets_remap);

for i=1:3
    sc                  = scatter(abs(acc_solo_dff(subsets_remap==i)), acc_better_minus_worse(subsets_remap==i));
    sc.MarkerFaceColor  = cmap(i,:);
    sc.MarkerEdgeColor  = 'none';
    sc.MarkerFaceAlpha  = .75;
end
fit_ecc                 = polyfit(abs(acc_solo_dff),acc_better_minus_worse,1);
pl                      = plot(x,polyval(fit_ecc, x),'LineWidth',1.5, 'Color',[.2 .2 .2]);
[r_acc,pv_acc]          = corrcoef(abs(acc_solo_dff),acc_better_minus_worse);
[r_acc_conv,pv_acc_conv] = corrcoef(abs(acc_solo_dff(conv)),acc_better_minus_worse(conv));
ax.XLim                 = [0 .15];
ax.YLim                 = [-.15 .15];
ax.YTick                = [-.15 -.1 -.05 0 .05 .1 .15];
ax.XLabel.String        = {'Solo interplayer distance'; 'Accuracy.better - Accuracy.worse'};
ax.YLabel.String        = {'Social modulation difference';'AUC.better - AUC.worse'};
ax.FontSize             = lb_fs;

text(.05,0,['All: r=' num2str(round(r_acc(2),3)) ' p=' num2str(pv_acc(2))])
text(.05,.05,['Conv: r=' num2str(round(r_acc_conv(2),3)) ' p=' num2str(pv_acc_conv(2))])
axis square

print(f, [dest_dir '/FIG_auc_dff_vs_solo_' alignment_str], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG_auc_dff_vs_solo_' alignment_str], '-r500', '-dpng', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE: Magnitude/step size of social modulation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                       = figure('units','centimeters','position',[0 0 25 10]);hold on
ax1                     = subplot(1,2,1);
[ax1,p1_better_conf]    = modulation_violin(ax1,raw.ecc1,raw.ecc2,auc.ecc1,auc.ecc2);
ax2                     = subplot(1,2,2);
[ax2,p1_better_acc]     = modulation_violin(ax2,raw.acc1,raw.acc2,auc.acc1,auc.acc2);

print(f, [dest_dir '/FIG_better_worse_modulation_' alignment_str], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG_better_worse_modulation_' alignment_str], '-r500', '-dpng', '-painters');

% CONSISTENT RELATIONSHIP?
sum((p1_better_acc - p1_better_conf) == 0) / length((p1_better_acc - p1_better_conf))

% %% Sanity check
% figure
% ax1                     = subplot(1,2,1);
% [ax1,p1_better_conf]    = modulation_violin(ax1,raw.ecc1,raw.ecc2,auc.acc1,auc.acc2);
% ax2                     = subplot(1,2,2);
% [ax2,p1_better_acc]     = modulation_violin(ax2,raw.acc1,raw.acc2,auc.ecc1,auc.ecc2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SFIGURE: Scatter dyadic modulation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                           = figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
ofs                         = .025;
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
        plot_scatter(ecc_df.solo,auc.ecc1,auc.ecc2)
        
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
        plot_scatter(acc_df.solo,auc.acc1,auc.acc2)
        
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
        ax.XLim                 = [-.15 .15];
        ax.YLim                 = [.3 .7];
        ax0v.YLim               = [.3 .7];

    % AUC distance between players        
    elseif iPlot == 3
        x                       = ecc_df.solo';
        y                       = abs(auc.ecc1-auc.ecc2)';
        sc                      = scatter(x,y,'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc.SizeData             = sc_size;
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
        [x_avg, y_avg]          = avg_win_value(x,y,win_size);
    elseif iPlot == 4
        x                       = acc_df.solo';
        y                       = abs(auc.acc1-auc.acc2)';
        sc                      = scatter(x,y,'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc.SizeData             = sc_size;
        [ft1,gof1]              = fit(x,y,'poly1');
        [ft2,gof2]              = fit(x,y,'poly2');
        acc_rsquared            = [gof1.adjrsquare gof2.adjrsquare];
        x_vec                   = -.1:.01:.1;
        p_fit                   = plot(x_vec,ft2(x_vec), 'Color', rcol, 'LineWidth',lw);
        tx                      = text(.001, .175,['adjR^{2}: ' num2str(round(gof2.adjrsquare,3))]);
        tx.FontSize             = lb_fs;
        ax.YLabel.String        = 'AUC difference';
        str                     = 'Accuracy';
        ax.XLim                 = [-.15 .15];
        ax.YLim                 = [0 .2];  
        [x_avg, y_avg]          = avg_win_value(x,y,win_size);
    end
    
    ax.Position(1)              = ax.Position(1) - ofs;
    ax.FontSize                 = lb_fs;
    ax.XLabel.String            = {'Within-dyad difference [P1 - P2]',str};
end

print(f, [dest_dir '/FIG_auc_' alignment_str], '-r500', '-dpng');
print(f, [dest_dir '/FIG_auc_' alignment_str], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SFIGURE: Solo vs dyadic correlation between players %% SUPPLEMENTARY
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
    sc.SizeData         = sc_size; 
    
    if i == 2 || i == 4
        ax.XLim     	= [.4 1];
        ax.YLim      	= [.4 1];
        ax.XTick        = [.4:.2:1];
        ax.YTick        = [.4:.2:1];
    else
        ax.XLim        	= [.8 1];
        ax.YLim       	= [.8 1];
        ax.XTick        = [.8:.1:1];
        ax.YTick        = [.8:.1:1];
    end
    
    ax.FontSize         = lb_fs;
    ax.XLabel.String    = 'Player1';
    ax.YLabel.String    = 'Player2';
    
    lsl                 = lsline;
    lsl.Color          	= rcol;
    lsl.LineWidth     	= lw;
    
    [r,pv]              = corrcoef(dat1,dat2);
    title({str, [' r = ' num2str(round(r(2),2)) ' pv = ' num2str(round(pv(2),2))]})
end

print(f, [dest_dir '/FIG_corr_solo_dyad_' alignment_str], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG_corr_solo_dyad_' alignment_str], '-r500', '-dpng', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SFIGURE: Solo difference vs average social modulation %% SUPPLEMENTARY
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
    sc.SizeData         = sc_size; 
    
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

print(f, [dest_dir '/FIG_corr_solo_auc_' alignment_str], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG_corr_solo_auc_' alignment_str], '-r500', '-dpng', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SFIGURE: Solo difference vs dyadic performance %% SUPPLEMENTARY
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
    sc.SizeData         = sc_size; 
    
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

print(f, [dest_dir '/FIG_corr_solo_score_' alignment_str], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG_corr_solo_score_' alignment_str], '-r500', '-dpng', '-painters');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dyadic Reward
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
s_score             = sum(score.solo,2);
d_score             = sum(score.dyad,2);
[p,h,z] = signrank(s_score, d_score);

vl                  = violinplot([s_score, d_score]);

for iDyad = 1:length(s_score)
    d_greater(iDyad)= s_score(iDyad) < d_score(iDyad);
    sd_dff(iDyad)   = s_score(iDyad) - d_score(iDyad);
end

ax = gca;
ax.FontSize = 24;
ax.YLabel.String = 'Sum of scores [P1 + P2]';
ax.XTick = [1 2];
ax.XTickLabel = {'Solo','Dyadic'};
ax.Title.String = {['Dyad > Solo: ' num2str((sum(d_greater)/length(d_greater))*100) '%; p<0.001']; ['Mean(Solo-Dyad): ' num2str(mean(sd_dff))]};
ax.Title.FontSize = ax.FontSize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reported stats in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Different from 0.5
median([auc.ecc1,auc.ecc2])
[p,h,stats] = signrank([auc.ecc1,auc.ecc2],.5);

median([auc.acc1,auc.acc2])
[p,h,stats] = signrank([auc.acc1,auc.acc2],.5);

% Dyadic improvement
ecc_class       = sum([auc.ecc1' auc.ecc2']>.5,2);
ecc_both_better = sum(ecc_class == 2) / length(ecc_class);
ecc_both_worse  = sum(ecc_class == 0) / length(ecc_class);
ecc_mixed       = sum(ecc_class == 1) / length(ecc_class);

acc_class       = sum([auc.acc1' auc.acc2']>.5,2);
acc_both_better = sum(acc_class == 2) / length(acc_class);
acc_both_worse  = sum(acc_class == 0) / length(acc_class);
acc_mixed       = sum(acc_class == 1) / length(acc_class);

% Convergence
sum(abs(ecc_df.solo) > abs(ecc_df.dyad)) / length(ecc_df.dyad)
sum(abs(acc_df.solo) > abs(acc_df.dyad)) / length(acc_df.dyad)

% Convergence - regression panel B 
sum(ecc_better_minus_worse<0) 
sum(acc_better_minus_worse<0) 
length(ecc_better_minus_worse) 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [acc_df, ecc_df, auc, score, hir, raw] = dyad_effect_size(in_solo, in_dyad, str)

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
    
    % Raw performance
    raw.acc1(iDyad)         = mean(in_solo{idx_ply1}.(['macc_' str])); % Mean(Median per coherence)
    raw.acc2(iDyad)         = mean(in_solo{idx_ply2}.(['macc_' str]));
    raw.ecc1(iDyad)         = mean(in_solo{idx_ply1}.(['mecc_' str]));
    raw.ecc2(iDyad)         = mean(in_solo{idx_ply2}.(['mecc_' str]));
        
    raw.dacc1(iDyad)        = mean(in_dyad{iDyad,1}.(['macc_' str]));
    raw.dacc2(iDyad)        = mean(in_dyad{iDyad,2}.(['macc_' str]));
    raw.decc1(iDyad)        = mean(in_dyad{iDyad,1}.(['mecc_' str]));
    raw.decc2(iDyad)        = mean(in_dyad{iDyad,2}.(['mecc_' str]));
    
    % Performance difference
    acc_df.solo(iDyad)      = mean(in_solo{idx_ply1}.(['macc_' str])) - mean(in_solo{idx_ply2}.(['macc_' str]));
    ecc_df.solo(iDyad)     	= mean(in_solo{idx_ply1}.(['mecc_' str])) - mean(in_solo{idx_ply2}.(['mecc_' str]));
    acc_df.dyad(iDyad)      = mean(in_dyad{iDyad,1}.(['macc_' str])) - mean(in_dyad{iDyad,2}.(['macc_' str]));
    ecc_df.dyad(iDyad)     	= mean(in_dyad{iDyad,1}.(['mecc_' str])) - mean(in_dyad{iDyad,2}.(['mecc_' str]));
       
    % Effect size: Solo vs Dyadic
    [auc.ecc1(iDyad), auc.ecc2(iDyad)] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, ['ecc_' str]);
    [auc.acc1(iDyad), auc.acc2(iDyad)] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, ['acc_' str]);

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

function plot_scatter(df,auc1,auc2)
sc_size              	= 10;
col                     = [.3 .3 .3];

for i = 1:length(df)
    ln                  = line([df(i) df(i)],[auc1(i) auc2(i)]);
    ln.Color            = col;
    ln.LineWidth        = .5;
    

    sc1                 = scatter(df(i),auc1(i),'MarkerFaceColor',col,'MarkerFaceAlpha',.5, 'MarkerEdgeColor',col,'MarkerEdgeAlpha',1,'Marker','o');
    sc2                 = scatter(df(i),auc2(i),'MarkerFaceColor','none','MarkerEdgeColor',col,'MarkerEdgeAlpha',1,'Marker','o');
    
    sc1.SizeData     	= sc_size;
    sc2.SizeData      	= sc_size;
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
end

function [ax, p1_better] = modulation_violin(ax,raw1,raw2,auc1,auc2)

p1_better = raw1 > raw2;
better_solo_player = [auc1(p1_better) auc2(~p1_better)];
worse_solo_player = [auc2(p1_better) auc1(~p1_better)];

disp('n: '); length(better_solo_player)
disp('signrank(better_solo_player,worse_solo_player)')
[p,h,stats]=signrank(better_solo_player,worse_solo_player)

disp('signrank(abs(better_solo_player-.5),abs(worse_solo_player-.5))')
[p,h,stats]=signrank(abs(better_solo_player-.5),abs(worse_solo_player-.5))

line([.5 2.5],[.5 .5], 'LineStyle',':','LineWidth', 1.5, 'Color', [0 0 0])
vl = violinplot([worse_solo_player' better_solo_player']);
vl(1).ViolinColor{1} = [.1 .1 .1];
vl(2).ViolinColor{1} = [.5 .5 .5];

for iDyad = 1:length(auc1)
    ln = line([vl(1).ScatterPlot.XData(iDyad) vl(2).ScatterPlot.XData(iDyad)],[vl(1).ScatterPlot.YData(iDyad) vl(2).ScatterPlot.YData(iDyad)]);
    ln.Color = [.8 .8 .8];
    ln.LineWidth = 1;
    uistack(ln, 'bottom');
end

ax.XLim = [.6 2.4];
ax.YLim = [0 1];
ax.XTick = [1 2];
ax.XTickLabel = {'Worse solo player', 'Better solo player'};
ax.YLabel.String = 'Social modulation [AUC]';
ax.FontSize = 8;

box off

bins                    = [0:.05:1];
h_ofs                   = .35;
ax0v                    = axes('Position', [ax.Position(1)+h_ofs ax.Position(2) ax.Position(3)/5 ax.Position(4)]); hold on
[v,edg]                 = histcounts([auc1,auc2],bins);
cntr                    = edg(2:end) - diff(edg)./2;
% st                      = stairs(v,cntr);
st                      = stairs(v,edg(2:end)); hold on
st.LineWidth            = 2;
st.Color                = [0 0 0];
mln                     = line([0 max(v)],[median([auc1,auc2]) median([auc1,auc2])]);
mln.LineWidth           = 2;
mln.Color                = [.6 0 0];
% ax0v.XAxis.Visible      = 'off';
% ax0v.YAxis.Visible      = 'off';
ax0v.YLim               = [0 1];

disp('n: '); length([auc1,auc2])

disp('signrank([auc1,auc2],.5)')
[p,h,stats]=signrank([auc1,auc2],.5)

disp('median([auc1,auc2])')
median([auc1,auc2])

end