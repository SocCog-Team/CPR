% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/comp_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hc_dyad_pairwise_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hc_dyad_pairwise_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_performance.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Within-dyad effect size %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[acc_df, ecc_df, auc, score, ply1_flag, ply2_flag, raw] = ...
    dyad_effect_size(solo_perf, dyad_pw_perf, dyad_perf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Scatter dyadic modulation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG4/SFIG4/';
f                           = figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
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
        plot_scatter(acc_df.solo',auc.ecc1',auc.ecc2',label_neg, ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Eccentricity [AUC: Solo vs Dyadic]';
        str                     = 'Accuracy';
        ax.XLim                 = [-.1 .1];
        ax.YLim                 = [.1 .9];
        
    elseif iPlot == 2
        plot_scatter(ecc_df.solo,auc.ecc1,auc.ecc2,label_neg,ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Eccentricity [AUC: Solo vs Dyadic]';
        str                     = 'Eccentricity';
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [.1 .9];
        
    elseif iPlot == 3
        plot_scatter(acc_df.solo,auc.acc1,auc.acc2,label_neg,ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Accuracy [AUC: Solo vs Dyadic]';
        str                     = 'Accuracy';
        ax.XLim                 = [-.1 .1];
        ax.YLim                 = [.3 .7];
        
    elseif iPlot == 4
        plot_scatter(ecc_df.solo,auc.acc1,auc.acc2,label_neg,ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Accuracy [AUC: Solo vs Dyadic]';
        str                     = 'Eccentricity';
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [.3 .7];
    end
    
    ax.Position(1)              = ax.Position(1) - ofs;
    ax.FontSize                 = lb_fs;
    ax.XLabel.String            = {'Within-dyad difference [P1 - P2]',str};
end

if label_neg
    print(f, [dest_dir '/SFIG4_auc_labelled'], '-r500', '-dpng');
    print(f, [dest_dir '/SFIG4_auc_labelled'], '-r500', '-dsvg', '-painters');
else
    print(f, [dest_dir '/SFIG4_auc'], '-r500', '-dpng');
    print(f, [dest_dir '/SFIG4_auc'], '-r500', '-dsvg', '-painters');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Scatter dyadic modulation corrected for Player1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG4/SFIG4/';
f                           = figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
lw                          = 1;
lb_fs                       = 8;
ofs                         = .025;
label_neg                   = false;
scol                        = [.1 .1 .1];   
rcol                        = [.6 .1 .1];

for iPlot = 1:4
    ax = subplot(2,2,iPlot); hold on
    ln                     	= line([0 0],[-.5 .5]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';

    if iPlot == 1
        x                       = acc_df.solo';
        y                       = abs(auc.ecc1-auc.ecc2)';
        sc                      = scatter(x,y,'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc.SizeData             = 10;
        [ft,gof]                = fit(x,y,'poly2');
        x_vec                   = -.35:.01:.35;
        p_fit                   = plot(x_vec,ft(x_vec), 'Color', rcol, 'LineWidth',lw);
        tx                      = text(.001, .4,['rsme' num2str(round(gof.rmse,3))]);
        tx.FontSize             = lb_fs;
        ax.YLabel.String        = 'AUC difference';
        str                     = 'Accuracy';
        ax.XLim                 = [-.1 .1];
        ax.YLim                 = [0 .5];
     
    elseif iPlot == 2
        x                       = ecc_df.solo';
        y                       = abs(auc.ecc1-auc.ecc2)';
        sc                      = scatter(x,y,'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc.SizeData             = 10;
        [ft,gof]                = fit(x,y,'poly2');
        x_vec                   = -.35:.01:.35;
        p_fit                   = plot(x_vec,ft(x_vec), 'Color', rcol, 'LineWidth',lw);
        tx                      = text(.001, .4,['rsme' num2str(round(gof.rmse,3))]);
        tx.FontSize             = lb_fs;
        ax.YLabel.String        = 'AUC difference';
        str                     = 'Eccentricity';
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [0 .5];
    elseif iPlot == 3
        x                       = acc_df.solo';
        y                       = abs(auc.acc1-auc.acc2)';
        sc                      = scatter(x,y,'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc.SizeData             = 10;
        [ft,gof]                = fit(x,y,'poly2');
        x_vec                   = -.1:.01:.1;
        p_fit                   = plot(x_vec,ft(x_vec), 'Color', rcol, 'LineWidth',lw);
        tx                      = text(.001, .175,['rsme' num2str(round(gof.rmse,3))]);
        tx.FontSize             = lb_fs;
        ax.YLabel.String        = 'AUC difference';
        str                     = 'Accuracy';
        ax.XLim                 = [-.1 .1];
        ax.YLim                 = [0 .2];  
    elseif iPlot == 4
        x                       = ecc_df.solo';
        y                       = abs(auc.acc1-auc.acc2)';
        sc                      = scatter(x,y,'MarkerFaceColor',scol,'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none');
        sc.SizeData             = 10;
        [ft,gof]                = fit(x,y,'poly2');
        x_vec                   = -.35:.01:.35;
        p_fit                   = plot(x_vec,ft(x_vec), 'Color', rcol, 'LineWidth',lw);
        tx                      = text(.001, .175,['rsme' num2str(round(gof.rmse,3))]);
        tx.FontSize             = lb_fs;
        ax.YLabel.String        = 'AUC difference';
        str                     = 'Eccentricity';
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [0 .2];  
    end
   
    ax.Position(1)              = ax.Position(1) - ofs;
    ax.FontSize                 = lb_fs;
    ax.XLabel.String            = {'Within-dyad difference [P1 - P2]',str};
end


if label_neg
    print(f, [dest_dir '/SFIG4_auc_labelled_corrP1'], '-r500', '-dpng');
    print(f, [dest_dir '/SFIG4_auc_labelled_corrP1'], '-r500', '-dsvg', '-painters');
else
    print(f, [dest_dir '/SFIG4_auc_corrP1'], '-r500', '-dpng');
    print(f, [dest_dir '/SFIG4_auc_corrP1'], '-r500', '-dsvg', '-painters');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Dyadic modulation - histogram %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG4/SFIG4/';
f                           = figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
lb_fs                       = 8;
ofs                         = .025;
nbin                        = 10;

for iPlot = 1:4
    ax                      = subplot(2,2,iPlot); hold on
    
    if iPlot == 1
        auc1                    = auc.ecc1 - .5;
        auc2                    = auc.ecc2 - .5;
        h                       = histogram([auc1 auc2],nbin,'FaceColor',[.1 .1 .1], 'EdgeColor','none');
        h.BinLimits             = [-.5 .5];
        ax.XLabel.String        = 'Eccentricity change';
        ax.YLabel.String        = '# Players';
        ax.XLim                 = h.BinLimits;
        ax.YLim                 = [0 max(h.BinCounts)];
    elseif iPlot == 2
        auc1                    = auc.acc1 - .5;
        auc2                    = auc.acc2 - .5;
        h                       = histogram([auc1 auc2],nbin,'FaceColor',[.1 .1 .1], 'EdgeColor','none');
        h.BinLimits             = [-.2 .2];
        ax.XLabel.String        = 'Accuracy change';
        ax.YLabel.String        = '# Players';
        ax.XLim                 = h.BinLimits;
        ax.YLim                 = [0 max(h.BinCounts)];
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
        sc.MarkerFaceColor      = [0 0 0];
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
            pl.Color         	= [.3 .3 .3];
            tx                 	= text(.4,.15,{['r=' num2str(round(r,2))];['p=' num2str(round(pv,2))]});
            tx.Color          	= [0 0 0];
            tx.FontSize      	= 8;
        end
    end
        
    ax.Position(1)              = ax.Position(1) - ofs;
    ax.FontSize                 = lb_fs;
end

if label_neg
    print(f, [dest_dir '/SFIG4_auc_labelled_hist'], '-r500', '-dpng');
    print(f, [dest_dir '/SFIG4_auc_labelled_hist'], '-r500', '-dsvg', '-painters');
else
    print(f, [dest_dir '/SFIG4_auc_hist'], '-r500', '-dpng');
    print(f, [dest_dir '/SFIG4_auc_hist'], '-r500', '-dsvg', '-painters');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Effect size difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subset_str = 'all';
[f, pv, r] = scatter_perf_dff(subset_str,ply1_flag,ply2_flag,auc,acc_df.solo,ecc_df.solo);

print(f, [dest_dir '/SFIG4_dff_scatter_' subset_str], '-r500', '-dpng');
print(f, [dest_dir '/SFIG4_dff_scatter_' subset_str], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Solo vs dyadic difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                      	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
ax                      = subplot(2,2,1); hold on
col                   	= {[230,97,1]./255;[94,60,153]./255};
marker                  = {'x','+'};

sc1                     = scatter(abs(ecc_df.solo),abs(ecc_df.dyad),'LineWidth',lw,'MarkerEdgeColor',col{1},'Marker',marker{1});
sc1.MarkerFaceAlpha     = .5;
sc2                     = scatter(abs(acc_df.solo),abs(acc_df.dyad),'LineWidth',lw,'MarkerEdgeColor',col{2},'Marker',marker{2});
sc2.MarkerFaceAlpha     = .5;

ax.XLim                 = [0 .4];
ax.YLim                 = [0 .4];

lsl                     = lsline;
lsl(1).Color            = col{2};
lsl(2).Color            = col{1};

ax.FontSize             = lb_fs;
ax.XLabel.String        = 'Difference [Solo]';
ax.YLabel.String        = 'Difference [Dyad]';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reported stats in paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regression to mean problem
nRep                        = 5000;
auc_str                     = 'ecc';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [acc_df, ecc_df, auc, score, ply1_flag, ply2_flag,raw] = dyad_effect_size(in_solo, in_dyad, dyad_pooled)

for iSubj = 1:size(in_solo,2)
    solo_id{iSubj}          = in_solo{iSubj}.id;
end

for iDyad = 1:size(in_dyad,1)
    score(iDyad,:)          = [mean(in_dyad{iDyad,1}.trg_mscore) mean(in_dyad{iDyad,2}.trg_mscore)];
    
    % Find solo data of subjects
    idx_ply1                = cellfun(@(x) strcmp(x,in_dyad{iDyad,1}.id),solo_id);
    idx_ply2                = cellfun(@(x) strcmp(x,in_dyad{iDyad,2}.id),solo_id);
    
    if strcmp(in_dyad{iDyad,1}.id,'agnt')
        idx_ply1 = idx_ply2;
    end
    
    % Label subjects of interest
    soi                     = ID_subj_negAUC(in_solo, dyad_pooled);
    ply1_flag(iDyad)      	= logical(sum(cellfun(@(x) strcmp(x,in_dyad{iDyad,1}.id),soi)));
    ply2_flag(iDyad)      	= logical(sum(cellfun(@(x) strcmp(x,in_dyad{iDyad,2}.id),soi)));
    
    % Raw performance
    raw.acc1(iDyad)         = mean(in_solo{idx_ply1}.macc_trg);
    raw.acc2(iDyad)         = mean(in_solo{idx_ply2}.macc_trg);
    raw.ecc1(iDyad)         = mean(in_solo{idx_ply1}.mecc_state);
    raw.ecc2(iDyad)         = mean(in_solo{idx_ply2}.mecc_state);
    
    % Performance difference
    acc_df.solo(iDyad)      = mean(in_solo{idx_ply1}.macc_trg) - mean(in_solo{idx_ply2}.macc_trg);
    ecc_df.solo(iDyad)     	= mean(in_solo{idx_ply1}.mecc_state) - mean(in_solo{idx_ply2}.mecc_state);
    acc_df.dyad(iDyad)      = mean(in_dyad{idx_ply1}.macc_trg) - mean(in_dyad{idx_ply2}.macc_trg);
    ecc_df.dyad(iDyad)     	= mean(in_dyad{idx_ply1}.mecc_state) - mean(in_dyad{idx_ply2}.mecc_state);
       
    % Effect size: Solo vs Dyadic
    [auc.ecc1(iDyad), auc.ecc2(iDyad)] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, 'ecc_state');
    [auc.acc1(iDyad), auc.acc2(iDyad)] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, 'acc_trg');
end
end

function [out1, out2] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, str)

for i = 1:length(in_solo{idx_ply1}.(str))
    out1(i)     = getAUROC(in_solo{idx_ply1}.(str){i},in_dyad{iDyad,1}.(str){i});
    out2(i)     = getAUROC(in_solo{idx_ply2}.(str){i},in_dyad{iDyad,2}.(str){i});
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
    tx                          = text(.25,height,{['r=' num2str(round(r,2))];['p=' num2str(round(pv,2))]});
    tx.Color                    = col;
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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% REMOVED: Performance convergence over time?
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for iSubj = 1:size(solo_perf,2)
%     solo_id{iSubj}          = solo_perf{iSubj}.id;
% end
% addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
% 
% for iDyad = 1:length(dyad_pw_perf) 
% %     % Find solo data of subjects
% %     idx_ply1                = cellfun(@(x) strcmp(x,dyad_pw_perf{iDyad,1}.id),solo_id);
% %     idx_ply2                = cellfun(@(x) strcmp(x,dyad_pw_perf{iDyad,2}.id),solo_id);
% % 
% %     % hit_idx ??? always confounded by evidence accumulation
% %     solo1 =  solo_perf{idx_ply1}.trg_all_ecc - solo_perf{idx_ply1}.trg_all_js_ecc100;
% %     solo2 =  solo_perf{idx_ply2}.trg_all_ecc - solo_perf{idx_ply2}.trg_all_js_ecc100;
% %     dyad1 =  dyad_pw_perf{iDyad,1}.trg_all_ecc - dyad_pw_perf{iDyad,1}.trg_all_js_ecc100;
% %     dyad2 =  dyad_pw_perf{iDyad,2}.trg_all_ecc - dyad_pw_perf{iDyad,2}.trg_all_js_ecc100;
% %     auc(iDyad,1) = getAUROC(solo1, dyad1) % is there a difference?
% %     auc(iDyad,2) = getAUROC(solo2, dyad2) % is there a difference?
% %     figure;plot(dyad1);lsline % Is there a slope?
% 
%     %High/Low confidence, high convergence of js
%     js_dev =  rad2deg(circ_dist(deg2rad(dyad_pw_perf{iDyad,1}.trg_all_js_dir),deg2rad(dyad_pw_perf{iDyad,2}.trg_all_js_dir)))
%     js_dev100 =  rad2deg(circ_dist(deg2rad(dyad_pw_perf{iDyad,1}.trg_all_js_dir100),deg2rad(dyad_pw_perf{iDyad,2}.trg_all_js_dir100)))
%     
%     eidx = dyad_pw_perf{iDyad,1}.trg_all_ecc > median(dyad_pw_perf{iDyad,1}.trg_all_ecc)
% %         figure; hold on
% %     scatter(abs(js_dev(eidx)),abs(js_dev100(eidx)))
%     
%     hit_vec = dyad_pw_perf{iDyad,1}.trg_all_outc & dyad_pw_perf{iDyad,2}.trg_all_outc;
% figure;
% subplot(2,2,1)
% scatter(abs(js_dev(hit_vec&eidx)), abs(dyad_pw_perf{iDyad,1}.trg_all_ecc(hit_vec&eidx) - dyad_pw_perf{iDyad,2}.trg_all_ecc(hit_vec&eidx)))
% title('Dual hit & High confidence')
% ylabel('abs JS ecc dff')
% lsline
% set(gca,'fontsize',16);
% 
% subplot(2,2,2)
% scatter(abs(js_dev(hit_vec&~eidx)), abs(dyad_pw_perf{iDyad,1}.trg_all_ecc(hit_vec&~eidx) - dyad_pw_perf{iDyad,2}.trg_all_ecc(hit_vec&~eidx)))
% title('Dual hit & Low confidence')
% lsline
% set(gca,'fontsize',16);
% 
% subplot(2,2,3)
% scatter(abs(js_dev(~hit_vec&eidx)), abs(dyad_pw_perf{iDyad,1}.trg_all_ecc(~hit_vec&eidx) - dyad_pw_perf{iDyad,2}.trg_all_ecc(~hit_vec&eidx)))
% title('Dual miss & High confidence')
% ylabel('abs JS ecc dff')
% xlabel('abs JS deviation at Target')
% lsline
% set(gca,'fontsize',16);
% 
% subplot(2,2,4)
% scatter(abs(js_dev(~hit_vec&~eidx)), abs(dyad_pw_perf{iDyad,1}.trg_all_ecc(~hit_vec&~eidx) - dyad_pw_perf{iDyad,2}.trg_all_ecc(~hit_vec&~eidx)))
% title('Dual miss & Low confidence')
% xlabel('abs JS deviation at Target')
% lsline
% set(gca,'fontsize',16);
% 
% end

