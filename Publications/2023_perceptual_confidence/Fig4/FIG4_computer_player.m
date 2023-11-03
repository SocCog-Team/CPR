% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/comp_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hc_dyad_pairwise_performance.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Within-dyad effect size %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[acc_df, ecc_df, auc_ecc1, auc_ecc2, auc_acc1, auc_acc2, ply1_flag, ply2_flag] = ...
    dyad_effect_size(solo_perf, hc_dyad_pw_perf, comp_perf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT: Scatter dyadic modulation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG4/';
f                           = figure('units','centimeters','position',[0 0 15 15]); hold on
lw                          = 1;
lb_fs                       = 8;
ofs                         = .025;
label_neg                   = false;

for iPlot = 1:4
    ax = subplot(2,2,iPlot); hold on
    ln                     	= line([0 0],[0 1]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    ln                     	= line([-.35 .35],[.5 .5]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    if iPlot == 1
        plot_scatter(acc_df',auc_ecc1',auc_ecc2',label_neg, ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Eccentricity [AUROC: Solo vs Dyadic]';
        str                     = 'Accuracy';
        ax.YLim                 = [.1  .9];
        
    elseif iPlot == 2
        plot_scatter(ecc_df,auc_ecc1,auc_ecc2,label_neg,ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Eccentricity [AUROC: Solo vs Dyadic]';
        str                     = 'Eccentricity';
        ax.YLim                 = [.1 .9];
        
    elseif iPlot == 3
        plot_scatter(acc_df,auc_acc1,auc_acc2,label_neg,ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Accuracy [AUROC: Solo vs Dyadic]';
        str                     = 'Accuracy';
        ax.YLim                 = [.3 .7];
        
    elseif iPlot == 4
        plot_scatter(ecc_df,auc_acc1,auc_acc2,label_neg,ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Accuracy [AUROC: Solo vs Dyadic]';
        str                     = 'Eccentricity';
        ax.YLim                 = [.3 .7];
    end
    
    ax.XLim                     = [-.4 .4];
    ax.Position(1)              = ax.Position(1) - ofs;
    ax.FontSize                 = lb_fs;
    ax.XLabel.String            = {'Within-dyad difference [P1 - P2]',str};
end

print(f, [dest_dir '/FIG4_HC_auc'], '-r500', '-dpng');
print(f, [dest_dir '/FIG4_HC_auc'], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Effect size difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subset_str = 'all';
[f, pv, r] = scatter_perf_dff(subset_str,ply1_flag,ply2_flag,auc_ecc1,auc_ecc2,auc_acc1,auc_acc2,acc_df,ecc_df);

print(f, [dest_dir '/FIG4_HC_dff_scatter_' subset_str], '-r500', '-dpng');
print(f, [dest_dir '/FIG4_HC_dff_scatter_' subset_str], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Performance convergence over time?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reported stats in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all = size(pv,1) * size(pv,2);
pv_corr = pv < (.05 / all);
r
pv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [acc_df, ecc_df, auc_ecc1, auc_ecc2, auc_acc1, auc_acc2, ply1_flag, ply2_flag] = dyad_effect_size(in_solo, in_dyad, in_comp)

for iSubj = 1:size(in_solo,2)
    solo_id{iSubj}          = in_solo{iSubj}.id;
end

for iDyad = 1:size(in_dyad,1)    
    % Find solo data of subjects
    idx                 = cellfun(@(x) strcmp(x,in_dyad{iDyad,2}.id),solo_id);
    
    % Label subjects of interest
    ply1_flag(iDyad)  	= false;
    ply2_flag(iDyad) 	= false;
    
    % Performance difference
    acc_df(iDyad)      	= mean(in_comp{idx}.macc_trg) - mean(in_solo{idx}.macc_trg);
    ecc_df(iDyad)     	= mean(in_comp{idx}.mecc_state) - mean(in_solo{idx}.mecc_state);
    
    % Effect size: Solo vs Dyadic
    [auc_ecc1(iDyad), auc_ecc2(iDyad)] = calcAUROC(in_solo, in_dyad, in_comp, idx, iDyad, 'ecc_state');
    [auc_acc1(iDyad), auc_acc2(iDyad)] = calcAUROC(in_solo, in_dyad, in_comp, idx, iDyad, 'acc_trg');
end
end

function [out1, out2] = calcAUROC(in_solo, in_dyad, in_comp, idx, iDyad, str)

for i = 1:length(in_solo{idx}.(str))
    
    comp_pooled = [];
    for iComp = 1:size(in_comp,2)
        if isempty(in_comp{iComp})
            continue
        end
        
        if contains(str, 'ecc')
            comp_pooled = [comp_pooled; in_comp{iComp}.(str){i}];
        else
            comp_pooled = [comp_pooled, in_comp{iComp}.(str){i}];
        end
    end
    
    out1(i)     = getAUROC(comp_pooled,in_dyad{iDyad,1}.(str){i});
    out2(i)     = getAUROC(in_solo{idx}.(str){i},in_dyad{iDyad,2}.(str){i});
end

out1 = mean(out1);
out2 = mean(out2);

end

function [out] = getAUROC(in1, in2)

if size(in1,1) == 1
    in1             = in1';
    in2             = in2';
end

lab                 = [zeros(length(in1),1); ones(length(in2),1)];
[~,~,~,out]         = perfcurve(lab,[in1; in2],1);

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
    tx                          = text(.25,height,{['r=' num2str(round(r,2))];['p=' num2str(round(pv,2))]});
    tx.Color                    = col;
    tx.FontSize                 = 8;
end

end

function [f, pv, r] = scatter_perf_dff(subset_str,ply1_flag,ply2_flag,auc_ecc1,auc_ecc2,auc_acc1,auc_acc2,acc_df,ecc_df)

f                               = figure('units','centimeters','position',[0 0 15 15]); hold on
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

auc_ecc1_subset             = auc_ecc1(idx);
auc_ecc2_subset             = auc_ecc2(idx);
auc_acc1_subset             = auc_acc1(idx);
auc_acc2_subset             = auc_acc2(idx);
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
        ax.XLim                 = [-.4 .4];
        ax.YLim                 = [-.55 .55];
        ax.YTick                = [-.5:.25:.5];
        ax.YLabel.String        = {'Eccentricity Difference';' [AUROC_P1 - AUROC_P2]'};
    elseif i == 2
        df_subset               = auc_acc1_subset - auc_acc2_subset;
        sc1                     = scatter(ecc_df_subset,df_subset,'LineWidth',lw,'MarkerEdgeColor',col{1},'Marker',marker{1});
        sc2                     = scatter(acc_df_subset,df_subset,'LineWidth',lw,'MarkerEdgeColor',col{2},'Marker',marker{2});
        [x_fit, y_fit, r(i,1), pv(i,1)]	= regr_line(ax,ecc_df_subset,df_subset,col{1});
        [x_fit, y_fit, r(i,2), pv(i,2)]	= regr_line(ax,acc_df_subset,df_subset,col{2},.13);
        ax.XLim                 = [-.4 .4];
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

