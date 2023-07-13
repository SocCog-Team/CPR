% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/WIP/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

% Import subject summary table
pth                         = '/Volumes/T7_Shield/CPR_psychophysics/';
x                           = readtable([pth 'Subjects_summary.xlsx']);
sbj_lst                     = x.Abbreviation;
sbj_lst(cellfun(@isempty,sbj_lst)) = [];
c                           = 0;
cc                          = 0;

sbj_lst(cellfun(@(x) strcmp(x, 'RoH'),sbj_lst)) = [];

% Extract file names of solo experiments
for iSubj = 1:length(sbj_lst)
    data_pth                = [pth sbj_lst{iSubj} '/summary/'];
    
    if isdir(data_pth)
        cd(data_pth)
        mat_files        	= dir('*.mat');
        
        for iFile = 1:length(mat_files)
            if contains(mat_files(iFile).name,'CPRsolo')
                c               = c+1;
                fname_solo{c} 	= mat_files(iFile).name;
            end
            if contains(mat_files(iFile).name,'CPRagent')
                cc              = cc+1;
                fname_agnt{cc} 	= mat_files(iFile).name;
            end
        end
    end
end

%% Extract subject data

nSample                                 = 30;
scnt                                    = 0;
nLag                                    = 150;

for iSub = 1:length(sbj_lst)
    
    for iExp = 1:2
        if iExp == 1
            fname                       = fname_solo;
        else
            fname                       = fname_agnt;
        end  
      
        fname_id                        = cellfun(@(x) x(1:12), fname, 'UniformOutput', false);
        fidx                            = find(cellfun(@(x) contains(x,lower(sbj_lst{iSub})),fname_id));
        tc                              = 0;
        t                               = [];
        all                             = [];
        tmp                             = {};
        
        if isempty(fidx)
            lag(iSub,iExp)              = nan;
            id{iSub,iExp}               = nan;
            continue
        end
        
        % Load experimental blocks
        for iBlock = 1:length(fidx)
            tbl                         = [];
            tbl                         = load([pth sbj_lst{iSub} '/summary/' fname{fidx(iBlock)}]);
            
            if sum(unique(tbl.t.rdp_coh) < .1) > 2
                continue
            end
            
            t                           = [t; tbl.t];                       % Concatenate data tables
            [tmp, tc]                   = extractTrials(tbl, tmp, tc);     	% Extract trials for correlation analysis
        end
        
        clear cr ps
        snr                         = unique(t.rdp_coh);
        id{iSub,iExp}            	= sbj_lst{iSub};
        nLag                        = 150;
        [cr,ps]                     = CPR_correlation_analysis_WIP(tmp, nLag, false);
        lag(iSub,iExp)           	= median(cr.lag);
        
        perf{iSub,iExp}           	= response_readout(t, nSample);
        
    end
    
    %%% Calculate difference between conditions
    if isempty(perf{iSub,2})
        continue
    end
    scnt                            = scnt+1;
    hir_df(scnt,:)                  = perf{iSub,2}.hir - perf{iSub,1}.hir;
    
    for iCoh = 1:length(snr)
        auc_acc(scnt,iCoh)          = f_auroc(perf{iSub,1}.acc_trg{iCoh},perf{iSub,2}.acc_trg{iCoh});
        auc_ecc(scnt,iCoh)          = f_auroc(perf{iSub,1}.ecc{iCoh},perf{iSub,2}.ecc{iCoh});
        auc_score(scnt,iCoh)        = f_auroc(perf{iSub,1}.trg_score{iCoh},perf{iSub,2}.trg_score{iCoh});
    end
    
    %%% Load all dyadic sessions of subjects
    dyad_tbl                        = getDyadicSessions(pth,sbj_lst{iSub});
    if ~isempty(dyad_tbl)
        perf{iSub,3}             	= response_readout(dyad_tbl, nSample);
    end
end

% Calculate area under receiver operating characteristics
% Value > .5 --> 2nd input located to the right of first
% input, i.e. more positive in this case
% a = normrnd(0,1,[1 1000]);
% b = normrnd(1,1,[1 1000]);
% figure; hold on
% histogram(a, 50)
% histogram(b, 50)
% f_auroc(a, b)

%% Load all agent sessions

t_agnt                              = [];
sxc                                 = [];
coh                                 = [];

for iSub = 1:length(sbj_lst)
    
    cd([pth sbj_lst{iSub} '/summary/'])
    mat_files                           = dir('*.mat');
    tc                                  = 0;
    tmp                                 = {};
    
    for iFile = 1:length(mat_files)
        if contains(mat_files(iFile).name,'agnt')
            tbl                         = load(mat_files(iFile).name);
            t_agnt                      = [t_agnt; tbl.t];
            
            if sum(unique(tbl.t.rdp_coh) < .1) > 2
                continue
            end
            
            [tmp, tc]                   = extractTrials(tbl, tmp, tc);     	% Extract trials for correlation analysis  
            id{iSub}                    = sbj_lst{iSub};
            nLag                        = 150;
            [cr,ps]                     = CPR_correlation_analysis_WIP(tmp, nLag, false);
            lag_agnt(iSub)            	= median(cr.lag);
        end
    end
  
    agnt(iSub)                       	= response_readout(t_agnt, nSample);
    sxc                                 = [sxc; cr.sxc];
    coh                                 = [coh; cr.coh'];
end

for iSess = 1:length(agnt)
    macc_agnt(iSess,:)  	= agnt(iSess).macc_trg;
    mecc_agnt(iSess,:)    	= agnt(iSess).mecc;
    hir_agnt(iSess,:)      	= agnt(iSess).hir;
    trg_score_agnt(iSess,:)	= agnt(iSess).trg_mscore;
end

%% PLOT

f                           = figure('units','normalized','position',[0 0 .5 1]);
height                    	= fliplr(linspace(.03,.8,4));
clmns                      	= linspace(.11,.75,3);
lb_fs                       = 14;
lg_fs                       = 10;
lw                          = 3;
frme_ms                     = 1000/120;
alp                         = .35;
avg_mult                    = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Example agent joystick response %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [0.2 0.2]*1.5;
row                         = .65;
clm                         = .14; 
hofs                        = .1;
cmap                        = [0 0 0; gray(256)];
steps                       = .05;
bins                        = 0:steps:1;
c                           = 0;

ax0v                       	= axes('Position', [clm row-hofs dim(1) dim(2)/5]); hold on
ax0h                        = axes('Position', [clm-hofs row dim(1)/5 dim(2)]); hold on
ax0                       	= axes('Position', [clm row dim]); hold on

dte                         = unique(t_agnt.date);
t_agnt_plot                 = t_agnt(t_agnt.date == dte(6),:);

% Extracte experimental data
clear trg_acc trg_conf trg_coh trg_hit
for iState = 1:size(t_agnt_plot,1)
    if t_agnt_plot.trg_shown(iState) == false
        continue
    end
    
    for iTarget = 1:length(t_agnt_plot.trg_ts{iState})
        c                   = c+1;
        trg_acc(c)          = t_agnt_plot.trg_acc{iState}(iTarget);
        trg_conf(c)         = t_agnt_plot.trg_ecc{iState}(iTarget);
        trg_coh(c)          = t_agnt_plot.rdp_coh(iState);
        trg_hit(c)          = t_agnt_plot.trg_hit{iState}(iTarget);
        
        if trg_acc(c) < .5 && trg_hit(c) == 1
            disp([num2str(iState) ' ' num2str(iTarget) ' ' num2str(c) ])
        end
    end
end

% Calculate reward matrix
acc                         = 0:.001:1;
conf                        = 0:.001:1;
rew                         = acc' .* conf;

% Determine arc width for each confidence level
for j = 1:length(conf)
    arc(j)                  = 180 - (180 * conf(j));
end

% Cap arc width at target width (2dva == 12.7587deg at chosen position)
aidx                        = arc < 12.7587;
arc(aidx)                   = 12.7587;

% For each confidence level, calculate minimum accuracy required to hit
% the target at given arc width - normalised values
hit_width_acc               = 1 - ((arc/2) / 180);
hit_width_acc(aidx)         = 1 - (12.7587/2)/180; % arc width fixed

% Remove position that cannot yield reward from reward matrix
for iAcc = 1:length(acc)
    indx                    = conf < hit_width_acc(iAcc);
    rew(iAcc,indx)          = nan;
end

% Plot reward matrix
hold on
im                          = imagesc(acc,conf,rew);
ax0.XLabel.String           = 'Accuracy [norm]';
ax0.YLabel.String           = 'Eccentricity [norm]';
ax0.FontSize                = lb_fs;
ax0.XLim                    = [0 1];
ax0.YLim                    = [0 1];
cb                          = colorbar;
cb.Label.String             = '% Reward';
cb.Location                 = 'eastoutside';
ax0.XLim                    = [0 1];
ax0.XTick                   = [0:.2:1];
ax0.YTick                   = [0:.2:1];
ax0.XTickLabelRotation      = 0;
ax0.YLabel.Position(1)      = -.35;
ax0.XLabel.Position(2)      = -.22;

colormap(cmap)

cmap_coh                    = cool(size(snr,1));

for iCoh = 1:length(snr)
    cidx                    = trg_coh == snr(iCoh);
    sc                      = scatter(trg_acc(cidx), trg_conf(cidx), 'filled');
    sc.CData            	= cmap_coh(iCoh,:);
    sc.SizeData             = 20;
    sc.MarkerFaceAlpha      = .9;
end

ax0.Position                = [clm row dim];

nBin                        = 40;
axes(ax0h)
[h, edg]                    = histcounts(trg_conf,nBin);
cntr                        = edg(1:end-1) + diff(edg) ./ 2;
st                          = stairs(-h,cntr);
st.LineWidth                = lw/1.5;
st.Color                    = [0 0 0];
ax0h.YLim                   = [0 1];
ax0h.XAxis.Visible          = 'off';
ax0h.YAxis.Visible          = 'off';

axes(ax0v)
[v, edg]                    = histcounts(trg_acc,nBin);
cntr                        = edg(1:end-1) + diff(edg) ./ 2;
st                          = stairs(cntr,-v);
st.LineWidth                = lw/1.5;
st.Color                    = [0 0 0];
ax0v.XLim                   = [0 1];
ax0v.XAxis.Visible          = 'off';
ax0v.YAxis.Visible          = 'off';
uistack(ax0v,'bottom')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Performance AGNT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yof                         = .075;
xof                         = .075;
ax4                         = axes('Position', [clm-xof row-yof dim(1)/1.25 dim(2)/2]); hold on
cmap                        = cbrewer('div', 'PRGn', 4, 'PCHIP');

[pl1, ci1]                 	= plotAvgCI(macc_agnt, cmap(1,:), alp, lw);
[pl2, ci2]               	= plotAvgCI(hir_agnt, cmap(2,:), alp, lw);
[pl3, ci3]                	= plotAvgCI(mecc_agnt, cmap(3,:), alp, lw);
[pl4, ci4]                	= plotAvgCI(trg_score_agnt, cmap(4,:), alp, lw);

ax4.YLim                    = [40 100];
ax4.XLim                    = [1 size(hir_agnt,2)];
ax4.XLabel.String           = 'Coherence [%]';
ax4.YLabel.String           = 'Avg performance [%]';
ax4.XTick                   = 1:length(snr);
ax4.XTickLabel              = round(snr,2)*100;
ax4.FontSize                = lb_fs;
ax4.XTickLabelRotation      = 0;
ax4.Position                = [clmns(3)-xof row-yof dim(1)/1.25 dim(2)/2];
% ax4.XAxis.Visible           = 'off';

box off

lg                          = legend([pl1 pl2 pl3 pl4],{'Accuracy', 'Hit rate','Eccentricity' 'Score'}, 'NumColumns', 2)';
lg.Box                      = 'off';
lg.FontSize                 = lg_fs;
lg.Location                 = 'northeast';
lg.Position(1)              = .675;
lg.Position(2)              = .67;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Crosscorrelation AGNT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yof                         = .16;
ax3                         = axes('Position', [clm-xof row+yof dim(1)/1.25 dim(2)/2]); hold on
cmap                        = cool(size(snr,1));

for iCoh = 1:size(snr,1)
    avg_sxc{iCoh} = mean(sxc(coh == snr(iCoh),:));
end

clear pl
for iCoh = 1:size(avg_sxc,2)
    mmat(iCoh, :)           = avg_sxc{iCoh};
    pl(iCoh)               	= plot(mmat(iCoh, :),'LineWidth', lw, 'Color', [cmap(iCoh,:) .75]);
    lg_str{iCoh}            = num2str(round(snr(iCoh)*100));
end

avg_peak_pos                = find(mean(mmat) == max(mean(mmat)));
px                          = plot(avg_peak_pos, max(mean(mmat)), 'kx');
px.MarkerSize               = 8;
px.LineWidth                = lw/2;

ln                          = line([avg_peak_pos avg_peak_pos],[0 max(mean(mmat))]);
ln.LineStyle                = ':';
ln.LineWidth                = lw/2;
ln.Color                    = [.5 .5 .5];

ln                          = line([150 150],[0 .2]);
ln.LineStyle                = ':';
ln.LineWidth                = lw;
ln.Color                    = [0 0 0];

ax3.YLabel.String           = 'XC Coef';
ax3.XLabel.String           = 'Lag [ms]';
ax3.XLim                    = [0 301];
ax3.YLim                    = [0 .16];
ax3.XTick                   = [0 150 avg_peak_pos 300];
ax3.FontSize              	= lb_fs;
ax3.XTickLabel              = round((cellfun(@str2num, ax3.XTickLabel)-nLag) * frme_ms);
ax3.Position                = [clmns(3)-xof row+yof dim(1)/1.25 dim(2)/2];
ax3.XTickLabelRotation      = 0;

lg                          = legend(pl,lg_str,'Location','northwest','NumColumns', 2);
lg.Box                      = 'off';
lg.TextColor                = [.99 .99 .99];
lg.FontSize                 = lg_fs;
lg.Position(1)              = .17;
lg.Position(2)              = .67;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBPLOT: Lag difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax6                         = axes('Position', [clmns(1) height(3) dim]); hold on
% 
% for iL = 1:length(lag)
%     pl(iL)               	= plot([1 2],[lag(iL,1) lag(iL,2)]);
%     pl(iL).Color          	= [.5 .5 .5 alp];
%     pl(iL).LineWidth      	= lw;
% end
% 
% bx                          = boxplot(lag, 'Colors', 'k');
% set(bx,'MarkerEdgeColor','k')
% set(bx, {'linew'},{lw})
% 
% ax6.YLabel.String           = 'Time [ms]';
% ax6.XTick                   = [1 2];
% ax6.XLim                    = [.5 2.5];
% ax6.XTickLabel              = {'Solo','Agnt'};
% ax6.FontSize                = lb_fs;
% ax6.XTickLabelRotation      = 0;
% ax6.Position                = [clmns(1) height(3) dim];
% ax6.Box                     = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Score difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [.2 .2];

ax6                         = axes('Position', [clmns(1) height(3) dim]); hold on
ln                          = line([0 1], [0 1]);
ln.LineStyle                = '--';
ln.LineWidth                = lw/2;
ln.Color                    = [0 0 0];


for iSubj = 1:size(perf,1)
    if ~isempty(perf{iSubj,2})
        sc_agnt            	= scatter(mean(cell2mat(perf{iSubj,2}.trg_score)),mean(cell2mat(perf{iSubj,1}.trg_score)), 'MarkerFaceColor', [.5 .5 .5],'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
    end
    
    if ~isempty(perf{iSubj,3})
        sc_dyad            	= scatter(mean(cell2mat(perf{iSubj,3}.trg_score)),mean(cell2mat(perf{iSubj,1}.trg_score)), 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
    end
end

ax6.XLabel.String           = {'Avg reward score', '[Dyad]'};
ax6.YLabel.String           = {'Avg reward score', '[Solo]'};
ax6.FontSize                = lb_fs;

lg                          = legend([sc_agnt(1) sc_dyad(1)], 'AGNT','HUMAN' );
lg.Location                 = 'southeast';
lg.FontSize                 = lg_fs;
lg.Box                      = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Hit rate comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSub = 1:size(perf,1)
    hir.solo(iSub)   	= perf{iSub,1}.hir_pool;
    if ~isempty(perf{iSub,2})
        hir.agnt(iSub) 	= perf{iSub,2}.hir_pool;
    else
        hir.agnt(iSub) 	= nan;
    end
    if ~isempty(perf{iSub,3})
        hir.dyad(iSub) 	= perf{iSub,3}.hir_pool;
    else
        hir.dyad(iSub) 	= nan;
  	end
end

nRep = 500;
e.solo                      = bootci(nRep, {@mean, hir.solo},'alpha', .001);
e.agnt                      = bootci(nRep, {@nanmean, hir.agnt},'alpha', .001);
e.dyad                      = bootci(nRep, {@nanmean, hir.dyad},'alpha', .001);

ax1                      	= axes('Position', [clmns(1) height(4) dim]); hold on
bp                          = bar([1:3],[mean(hir.solo) nanmean(hir.agnt) nanmean(hir.dyad)]);
bp.FaceColor                = [.5 .5 .5];
bp.EdgeColor                = 'none';

er                          = errorbar([1:3],[mean(hir.solo) nanmean(hir.agnt) nanmean(hir.dyad)],[nanmean(hir.solo)-e.solo(2) nanmean(hir.agnt)-e.agnt(2) nanmean(hir.dyad)-e.dyad(2)],[e.solo(1)-mean(hir.solo) e.agnt(1)-nanmean(hir.agnt) e.dyad(1)-nanmean(hir.dyad)]);    
er.Color                    = [0 0 0];                            
er.LineStyle                = 'none';
er.LineWidth                = lw/1.5;

ax1.YLabel.String           = 'Avg hit rate';
ax1.XLim                    = [.5 3.5];
ax1.YLim                    = [.3 .6];
ax1.YTick                 	= .3:.1:.6;
ax1.XTickLabel              = {'Solo','Agnt','Dyad'};
ax1.FontSize                = lb_fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Hit rate difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax10                      	= axes('Position', [clmns(2) height(3) dim]); hold on

ln                          = line([0 8],[0 0], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax10                        = plotData(ax10, hir_df, 'Hit rate diff', snr, alp, lw, lb_fs, avg_mult,true);
ax10.YLim                 	= [-.1 .4];
ax10.YLim                 	= [-.2 .4];
ax10.YTick                  = [-.2 0 .2 .4];
ax10.Position             	= [clmns(2) height(3) dim];
ax10.XAxis.Visible          = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Score difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax11                    	= axes('Position', [clmns(3) height(3) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax11                        = plotData(ax11, auc_score, {'Score diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax11.Position           	= [clmns(3) height(3) dim];
ax11.XAxis.Visible          = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Eccentricity difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax12                        = axes('Position', [clmns(2) height(4) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax12                        = plotData(ax12, auc_ecc, {'Eccentricity diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax12.Position           	= [clmns(2) height(4) dim];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Accuracy difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax13                     	= axes('Position', [clmns(3) height(4) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax13                        = plotData(ax13, auc_acc, {'Accuracy diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax13.Position             	= [clmns(3) height(4) dim];

text(3.75,0.825, 'Agnt > Solo', 'Parent', ax13, 'FontSize', lb_fs, 'Color', [0 0 0])
text(3.75,0.175, 'Agnt < Solo', 'Parent', ax13, 'FontSize', lb_fs, 'Color', [0 0 0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax0                         = axes('Position',[0 0 1 1],'Visible','off');
lofs                        = .225;
text(clmns(1)-.1,height(1)+lofs, 'A', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(clmns(3)-.175,height(1)+lofs, 'B', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(clmns(3)-.175,height(2)+lofs, 'C', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(clmns(1)-.1,height(3)+lofs, 'D', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(clmns(2)-.1,height(3)+lofs, 'E', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')

text(0.2,.75, 'Coherence', 'Parent', ax0, 'FontSize', lg_fs, 'Color', [.99 .99 .99])

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG3/FIG3', '-r400', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG3/FIG3', '-r400', '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tmp, tc] = extractTrials(tbl, tmp, tc)
for iTrl = 1:tbl.t.cyc_no(end)
    tidx                    = tbl.t.cyc_no == iTrl;             % Trial index
    tt                      = tbl.t(tidx,:);                    % Relevant stimulus cycles
    
    % Initiate cell
    tc                      = tc+1;                             % Trial counter
    tmp.frme_ts{tc}         = [];
    tmp.rdp_dir{tc}         = [];
    tmp.rdp_coh{tc}     	= [];
    tmp.js_dir{tc}      	= [];
    tmp.js_ecc{tc}      	= [];
    tmp.refresh{tc}         = [];
    
    % Extract and add data
    for iState = 1:size(tt,1)
        tmp.frme_ts{tc}  	= [tmp.frme_ts{tc} tt.frme_ts{iState}];
        tmp.rdp_dir{tc}  	= [tmp.rdp_dir{tc} repmat(tt.rdp_dir(iState),1,length(tt.frme_ts{iState}))];
        tmp.rdp_coh{tc}  	= [tmp.rdp_coh{tc} repmat(tt.rdp_coh(iState),1,length(tt.frme_ts{iState}))];
        tmp.js_dir{tc}  	= [tmp.js_dir{tc} tt.js_dir{iState}];
        tmp.js_ecc{tc}  	= [tmp.js_ecc{tc} tt.js_ecc{iState}];
        tmp.refresh{tc}     = [tmp.refresh{tc} median(diff(tmp.frme_ts{tc}))];
    end
end
end

function dyad_tbl = getDyadicSessions(pth,ID)
    dyad_cnt                    = 0;
    dyad_tbl                    = [];
    
    for iDyad = 19:61
        cd([pth ['Dyad' num2str(iDyad)] '/summary/'])
        mat_files              	= dir('*.mat');
        
        if isempty(mat_files)
            continue
        end
        
        dyad_cnt                = dyad_cnt+1;
        tmp1                    = split(mat_files(1).name,'_');
        tmp3                    = split(mat_files(3).name,'_');
        id_dyad(dyad_cnt,:)    	= [{tmp1{2}}, {tmp3{2}}];
        n_dyad(dyad_cnt,:)      = iDyad;
        
        % Check if subject part of this dyad
        for iFile = 1:size(mat_files,1)
            is_player(iFile) = contains(mat_files(iFile).name,lower(ID));
            
            if is_player(iFile)
                tmp_tbl      	= load(mat_files(iFile).name);
                dyad_tbl      	= [dyad_tbl; tmp_tbl.t];                % Organise dyadic data
            end
        end
    end
end

function out = response_readout(in, nSample)

c                           = 0;
snr                         = unique(in.rdp_coh);

% Target score
clear tscore score_cum score_hi score_coh score_dte score_exp trg_states
for iState = 1:size(in.trg_ts,1)
    for iTrg = 1:length(in.trg_ts{iState})
        try
            c                   = c +1;
            score_hi(c)         = in.trg_hit{iState}(iTrg);
            score(c)            = double(in.trg_score{iState}(iTrg));
            score_coh(c)        = in.rdp_coh(iState);
            score_dte(c)        = in.date(iState);
            score_exp(c)        = in.exp(iState);
        catch
            warning(['Skipped state/target: ' num2str(iState) '/' num2str(iTrg)])
        end
    end
end

trg_states                 	= in.trg_hit(logical(in.trg_shown));
out.hir_pool                = sum(cellfun(@sum,trg_states)) / sum(cellfun(@numel,trg_states));

dte = unique(score_dte);
dte(ismissing(dte)) = [];

for iDate = 1:length(dte)
    clear dte_idx score_cum dte_score
    dte_idx                 = score_dte == dte(iDate);
    out.score_final_exp(iDate) = unique(score_exp(dte_idx));
    dte_score               = score(dte_idx);
    score_cum               = cumsum(dte_score(~isnan(dte_score)));
    out.score_norm(iDate)   = score_cum(end) ./ length(dte_score(~isnan(dte_score)));
end

for iCoh = 1:length(snr)
    clear cIdx tIdx nhi ntrg
    
    cIdx = in.rdp_coh == snr(iCoh);
    tIdx = logical(in.trg_shown);
    
    tIdx(cellfun(@length,in.js_ecc) < 100) = false;
    cIdx(cellfun(@length,in.js_ecc) < 100) = false;
    
    % Hit rate
    nhi                     = sum(cellfun(@sum,in.trg_hit(cIdx & tIdx)));
    ntrg                    = sum(cellfun(@numel,in.trg_hit(cIdx & tIdx)));
    out.hir(iCoh)           = nhi / ntrg;
    
    % Target score [hits only]
    out.trg_mscore(iCoh)	= nanmean(score(score_coh  == snr(iCoh) & score_hi == true));
    out.trg_score{iCoh}  	= score(score_coh  == snr(iCoh) & score_hi == true);
    
    % Joystick displacement
    out.mecc(iCoh)         	= nanmedian(cellfun(@(x) nanmedian(x(end-nSample:end)), in.js_ecc(cIdx)));
    out.ecc{iCoh}           = cellfun(@(x) nanmedian(x(end-nSample:end)), in.js_ecc(cIdx));
    
    % Joystick accuracy [state-wise]
    for iState = 1:length(in.rdp_dir)
        % At least 100 frames
        if length(in.js_dir{iState}) < 100
            continue
        end
        js_dev              = rad2deg(circ_dist(deg2rad(in.js_dir{iState}(end-nSample:end)),deg2rad(in.rdp_dir(iState))));  % Minimum RDP-Joystick difference
        js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));         	% Joystick accuracy
    end
    out.macc_state(iCoh)  	= nanmedian(js_acc);
    out.acc_state{iCoh}    	= js_acc;
    
    % Joystick accuracy [before first target]
    t1_ts                   = cellfun(@(x) x(1), in.trg_ts);
    f1_ts                   = cellfun(@(x) x(1), in.frme_ts);
    trgIdx                  = (t1_ts-f1_ts) >= 1e6;
    rdp_dir                 = in.rdp_dir(cIdx & in.trg_shown & trgIdx);
    js_dir                  = in.js_dir(cIdx & in.trg_shown & trgIdx);
    frmes                   = in.frme_ts(cIdx & in.trg_shown & trgIdx);
    trg1_ts                 = t1_ts(cIdx & in.trg_shown & trgIdx);
    
    clear js_acc
    for iState = 1:length(rdp_dir)
        clear js_dev
        smpl_idx            = find(frmes{iState} < trg1_ts(iState),1,'last')-nSample : find(frmes{iState} < trg1_ts(iState),1,'last');
        js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iState}(smpl_idx)),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
        js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));         	% Joystick accuracy
    end
    
    out.macc_trg(iCoh)      = nanmedian(js_acc);
    out.acc_trg{iCoh}     	= js_acc;
    out.carr(iCoh)         	= snr(iCoh);
end
end

function [pl, ci] = plotAvgCI(dat, col, alp, lw)

ci                	= bootci(500,{@median, dat.*100}, 'alpha', .01);
pt                  = patch([1:size(dat,2) fliplr(1:size(dat,2))], [ci(1,:) fliplr(ci(2,:))], col, 'FaceAlpha',alp, 'EdgeColor','none');

pl                 	= plot(median(dat).*100);
pl.Color           	= col;
pl.Marker       	= '.';
pl.MarkerSize       = 15;
pl.LineWidth       	= lw/2;
end

function [ax] = plotData(ax, dat, ylab, snr, alp, lw, lb_fs, avg_mult, flag)

if nargin < 9
    flag = false;
end

% if flag == false
%     pt                      = patch([0 8 8 0], [.5 .5 1 1], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
%     pt                      = patch([0 8 8 0], [0 0 .5 .5], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
% else
%     pt                      = patch([0 8 8 0], [0 0 .5 .5], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
%     pt                      = patch([0 8 8 0], [-.5 -.5 0 0], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
% end

for iL = 1:size(dat,1)
    pl(iL)               	= plot(dat(iL,:));
    pl(iL).Color          	= [.5 .5 .5 alp];
    pl(iL).LineWidth      	= lw/2;
end

% Boostrap confidence intervals
nRep                        = 1000;
[CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.05);

% Prepare filled area
vec                         = 1:length(CI);
x_spacing                   = [vec fliplr(vec)];
ci                          = [CI(1,:) fliplr(CI(2,:))];

% Overlay confidence intervals
fl                          = fill(x_spacing,ci,[.3 0 0],'EdgeColor','none', 'FaceAlpha', alp);

% Plot mean curve
pm                          = plot(nanmedian(dat),'LineWidth', lw*avg_mult, 'Color', [0 0 0]);

ax.XLim                     = [1 size(dat,2)];
ax.XLabel.String            = 'Coherence [%]';
ax.YLabel.String            = ylab;
ax.XTick                    = [1 4 length(snr)];
ax.XTickLabel               = round([snr(1) snr(4) snr(end)],2)*100;
ax.FontSize                 = lb_fs;
ax.YLim                     = [.1 .9];
ax.YTick                    = [.25 .5 .75];
ax.XTickLabelRotation       = 0;

box off
end

