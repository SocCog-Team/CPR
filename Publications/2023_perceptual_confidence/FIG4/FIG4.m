% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hc_dyad_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hc_dyad_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/comp_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/comp_performance.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Within-dyad effect size %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSubj = 1:size(solo_perf,2)
    solo_id{iSubj}          = solo_perf{iSubj}.id;
end

for iDyad = 1:size(dyad_pw_perf,1)
    score(iDyad,:)          = [mean(dyad_pw_perf{iDyad,1}.trg_mscore) mean(dyad_pw_perf{iDyad,2}.trg_mscore)];
    
    
    % Find solo data of subjects
    idx_ply1                = cellfun(@(x) strcmp(x,dyad_pw_perf{iDyad,1}.id),solo_id);
    idx_ply2                = cellfun(@(x) strcmp(x,dyad_pw_perf{iDyad,2}.id),solo_id);
    
    % Label subjects of interest
    soi                     = ID_subj_negAUC(solo_perf, dyad_perf);
    ply1_flag(iDyad)      	= logical(sum(cellfun(@(x) strcmp(x,dyad_pw_perf{iDyad,1}.id),soi)));
    ply2_flag(iDyad)      	= logical(sum(cellfun(@(x) strcmp(x,dyad_pw_perf{iDyad,2}.id),soi)));
    
    % Performance difference
    %     acc_df(iDyad)           = calcDiff(dyad_pw_perf,iDyad,'acc_trg');
    %     ecc_df(iDyad)           = calcDiff(dyad_pw_perf,iDyad,'ecc');
    acc_df(iDyad)           = mean(solo_perf{idx_ply1}.macc_trg) - mean(solo_perf{idx_ply2}.macc_trg);
    ecc_df(iDyad)           = mean(solo_perf{idx_ply1}.mecc) - mean(solo_perf{idx_ply2}.mecc);
    
    % AUROC: solo vs dyadic
    [auc_ecc1(iDyad), auc_ecc2(iDyad)] = calcAUROC(solo_perf, dyad_pw_perf, idx_ply1, idx_ply2, iDyad, 'ecc');
    [auc_acc1(iDyad), auc_acc2(iDyad)] = calcAUROC(solo_perf, dyad_pw_perf, idx_ply1, idx_ply2, iDyad, 'acc_trg');
end

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG4/';
lw                          = 1;
lb_fs                       = 8;
f                           = figure('units','centimeters','position',[0 0 15 15]); hold on
ofs                         = .025;

for iResponse = 1:4
    ax = subplot(2,2,iResponse); hold on
    ln                     	= line([0 0],[0 1]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    ln                     	= line([-.35 .35],[.5 .5]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    if iResponse == 1
        % Contrast eccentricity solo vs dyadic for each player
        plot_scatter(acc_df',auc_ecc1',auc_ecc2',score,ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Eccentricity [AUROC: Solo vs Dyadic]';
        str                     = 'Accuracy';
        ax.XLim                 = [-.1 .1];
        ax.YLim                 = [.1 .9];
        
    elseif iResponse == 2
        plot_scatter(ecc_df,auc_ecc1,auc_ecc2,score,ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Eccentricity [AUROC: Solo vs Dyadic]';
        str                     = 'Eccentricity';
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [.1 .9];
        
    elseif iResponse == 3
        % Contrast accuracy solo vs dyadic for each player
        plot_scatter(acc_df,auc_acc1,auc_acc2,score,ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Accuracy [AUROC: Solo vs Dyadic]';
        str                     = 'Accuracy';
        ax.XLim                 = [-.1 .1];
        ax.YLim                 = [.3 .7];
        
    elseif iResponse == 4
        plot_scatter(ecc_df,auc_acc1,auc_acc2,score,ply1_flag,ply2_flag)
        ax.YLabel.String        = 'Accuracy [AUROC: Solo vs Dyadic]';
        str                     = 'Eccentricity';
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [.3 .7];
    end
    
    ax.Position(1)              = ax.Position(1) - ofs;
    ax.FontSize                 = lb_fs;
    ax.XLabel.String            = {'Within-dyad difference [P1 - P2]',str};
end

% cb                              = colorbar;
% cb.Position(1:2)                = [.93 .75];
% cb.Position(4)                  = .1;
% cb.Title.String                 = {'Hit score [%]'};
% cb.FontSize                     = lb_fs;
% cb.TickLabels                   = {'25', '50', '75'};
% n                               = 100;
% col                             = [linspace(202/255,5/255,n)',linspace(0/255,113/255,n)',linspace(32/255,176/255,n)'];
% colormap(col)

print(f, [dest_dir '/FIG4_auc'], '-r500', '-dpng');
print(f, [dest_dir '/FIG4_auc'], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Effect size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f                           = figure('units','centimeters','position',[0 0 15 15]); hold on
marker                      = {'x','+'};
col                         = {[230,97,1]./255;[94,60,153]./255};

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
        df                      = auc_ecc1 - auc_ecc2;
        sc1                     = scatter(ecc_df,df,'LineWidth',lw,'MarkerEdgeColor',col{1},'Marker',marker{1});
        sc2                     = scatter(acc_df,df,'LineWidth',lw,'MarkerEdgeColor',col{2},'Marker',marker{2});
        [x_fit, y_fit, r(i,1), pv(i,1)]	= regr_line(ax,ecc_df,df,col{1});
        [x_fit, y_fit, r(i,2), pv(i,2)]	= regr_line(ax,acc_df,df,col{2},.13);
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [-.55 .55];
        ax.YTick                = [-.5:.25:.5];
        ax.YLabel.String        = {'Eccentricity Difference';' [AUROC_P1 - AUROC_P2]'};
    elseif i == 2
        df                      = auc_acc1 - auc_acc2;
        sc1                     = scatter(ecc_df,df,'LineWidth',lw,'MarkerEdgeColor',col{1},'Marker',marker{1});
        sc2                     = scatter(acc_df,df,'LineWidth',lw,'MarkerEdgeColor',col{2},'Marker',marker{2});
        [x_fit, y_fit, r(i,1), pv(i,1)]	= regr_line(ax,ecc_df,df,col{1});
        [x_fit, y_fit, r(i,2), pv(i,2)]	= regr_line(ax,acc_df,df,col{2},.13);
        ax.XLim                 = [-.35 .35];
        ax.YLim                 = [-.25 .25];
        ax.YTick                = [-.2:.1:.2];
        ax.YLabel.String        = {'Accuracy Difference';' [AUROC_P1 - AUROC_P2]'};
    end
    
    ax.YLabel.Interpreter       = 'none';
    ax.Position(1)              = ax.Position(1) - ofs;
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

print(f, [dest_dir '/FIG4_dff_scatter'], '-r500', '-dpng');
print(f, [dest_dir '/FIG4_dff_scatter'], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exclude subjects with negative eccentricity AUC from FIG2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
marker                          = {'x','+'};
col                             = {[230,97,1]./255;[94,60,153]./255};

for iSubset = 1:2
    f                       	= figure('units','centimeters','position',[0 0 15 15]); hold on

    if iSubset == 1
        idx                   	= ~ply1_flag & ~ply2_flag; % pos
    else
        idx                     = ply1_flag | ply2_flag; % neg
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
            [x_fit, y_fit, r(i,1,iSubset), pv(i,1,iSubset)]	= regr_line(ax,ecc_df_subset,df_subset,col{1});
            [x_fit, y_fit, r(i,2,iSubset), pv(i,2,iSubset)]	= regr_line(ax,acc_df_subset,df_subset,col{2},.13);
            ax.XLim                 = [-.35 .35];
            ax.YLim                 = [-.55 .55];
            ax.YTick                = [-.5:.25:.5];
            ax.YLabel.String        = {'Eccentricity Difference';' [AUROC_P1 - AUROC_P2]'};
        elseif i == 2
            df_subset               = auc_acc1_subset - auc_acc2_subset;
            sc1                     = scatter(ecc_df_subset,df_subset,'LineWidth',lw,'MarkerEdgeColor',col{1},'Marker',marker{1});
            sc2                     = scatter(acc_df_subset,df_subset,'LineWidth',lw,'MarkerEdgeColor',col{2},'Marker',marker{2});
            [x_fit, y_fit, r(i,1,iSubset), pv(i,1,iSubset)]	= regr_line(ax,ecc_df_subset,df_subset,col{1});
            [x_fit, y_fit, r(i,2,iSubset), pv(i,2,iSubset)]	= regr_line(ax,acc_df_subset,df_subset,col{2},.13);
            ax.XLim                 = [-.35 .35];
            ax.YLim                 = [-.25 .25];
            ax.YTick                = [-.2:.1:.2];
            ax.YLabel.String        = {'Accuracy Difference';' [AUROC_P1 - AUROC_P2]'};
        end
        
        ax.YLabel.Interpreter       = 'none';
        ax.Position(1)              = ax.Position(1) - ofs;
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
    
    if iSubset == 1
        print(f, [dest_dir '/FIG4_dff_scatter_subset_pos'], '-r500', '-dpng');
        print(f, [dest_dir '/FIG4_dff_scatter_subset_pos'], '-r500', '-dsvg', '-painters');
    else
        print(f, [dest_dir '/FIG4_dff_scatter_subset_neg'], '-r500', '-dpng');
        print(f, [dest_dir '/FIG4_dff_scatter_subset_neg'], '-r500', '-dsvg', '-painters');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Performance convergence over time?
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
function out = calcDiff(in_dyad,iDyad,str)

% Extract within-dyad accuracy difference
if strcmp(str, 'ecc')
    dat1 = mean(cell2mat(in_dyad{iDyad,1}.(str)'));
    dat2 = mean(cell2mat(in_dyad{iDyad,2}.(str)'));
elseif strcmp(str, 'hir_pool')
    dat1 = in_dyad{iDyad,1}.(str);
    dat2 = in_dyad{iDyad,2}.(str);
else
    dat1 = mean(cell2mat(in_dyad{iDyad,1}.(str)));
    dat2 = mean(cell2mat(in_dyad{iDyad,2}.(str)));
end

out = dat1-dat2;

end

function [out1, out2] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, str)

% if strcmp(str, 'ecc')
%     out1            = getAUROC(cell2mat(in_solo{idx_ply1}.(str)'),cell2mat(in_dyad{iDyad,1}.(str)'));
%     out2            = getAUROC(cell2mat(in_solo{idx_ply2}.(str)'),cell2mat(in_dyad{iDyad,2}.(str)'));
% else
%     out1            = getAUROC(cell2mat(in_solo{idx_ply1}.(str)),cell2mat(in_dyad{iDyad,1}.(str)));
%     out2            = getAUROC(cell2mat(in_solo{idx_ply2}.(str)),cell2mat(in_dyad{iDyad,2}.(str)));
% end

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

function plot_scatter(df,auc1,auc2,score,ply1_flag,ply2_flag)

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
    
    
    %     ln                  = line([df(i) df(i)],[auc1(i) auc2(i)]);
    %     ln.Color            = [.5 .5 .5];
    %     ln.LineWidth        = .5;
    %
    %     n = 50;
    %
    %     col                 = [linspace(202/255,5/255,n)',linspace(0/255,113/255,n)',linspace(32/255,176/255,n)'];
    %     col_x               = [repmat(col(1,:),n/2,1);col;repmat(col(end,:),n/2,1)];
    %     col_score1          = int64(round(score(i,1),2)*100);
    %     col_score2          = int64(round(score(i,2),2)*100);
    %
    %     sc1 = scatter(df(i),auc1(i),'MarkerFaceColor',col_x(col_score1,:),'MarkerFaceAlpha',.85,'MarkerEdgeColor','none','Marker','^');
    %     sc2 = scatter(df(i),auc2(i),'MarkerFaceColor',col_x(col_score2,:),'MarkerFaceAlpha',.85,'MarkerEdgeColor','none','Marker','o');
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

function out = ID_subj_negAUC(solo_perf, dyad_perf)

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
        auc_ecc(cnt,iCoh)	= getAUROC(d.sp.ecc{iCoh},d.dp.ecc{iCoh});
    end
end

% Identify subjects
out                         = id_auc(mean(auc_ecc,2)<.5);

end