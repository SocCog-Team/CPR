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
s_cnt                                   = 0;
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
        
        % Target score
        clear tscore score_cum score_hi score_coh
        c                           = 0;
        for iState = 1:size(t.trg_ts,1)
            for iTrg = 1:length(t.trg_ts{iState})
                c                   = c +1;
                score_cum(c)        = t.trg_score{iState}(iTrg);
                score_coh(c)        = t.rdp_coh(iState);
                score_hi(c)         = t.trg_hit{iState}(iTrg);
            end
        end
        
        tmp_score                   = score_cum(~isnan(score_cum));
        tscore(~isnan(score_cum))   = [0 diff(tmp_score)];
        tscore(tscore < 0)          = nan;
        
        
        trg_states                 	= t.trg_hit(logical(t.trg_shown));
        hir_pool(iSub,iExp)       	= sum(cellfun(@sum,trg_states)) / sum(cellfun(@numel,trg_states));

        for iCoh = 1:length(snr)
            clear cIdx tIdx nhi ntrg
            
            cIdx = t.rdp_coh == snr(iCoh);
            tIdx = logical(t.trg_shown);
            
            tIdx(cellfun(@length,t.js_str) < 100) = false;
            cIdx(cellfun(@length,t.js_str) < 100) = false;
            
            % Hit rate
            nhi                     = sum(cellfun(@sum,t.trg_hit(cIdx & tIdx)));
            ntrg                    = sum(cellfun(@numel,t.trg_hit(cIdx & tIdx)));
            hir(iSub,iExp,iCoh)  	= nhi / ntrg;
            
            % Target score
            trg_score{iSub,iExp,iCoh} 	= tscore(score_coh  == snr(iCoh) & score_hi == true);
            
            % Joystick displacement
            mstr{iSub,iExp,iCoh}        = cellfun(@(x) nanmedian(x(end-nSample:end)), t.js_str(cIdx));
            
            % Joystick accuracy before first target
            t1_ts                   = cellfun(@(x) x(1), t.trg_ts);
            f1_ts                   = cellfun(@(x) x(1), t.frme_ts);
            trgIdx                  = (t1_ts-f1_ts) > 1e6;
            rdp_dir                 = t.rdp_dir(cIdx & t.trg_shown & trgIdx);
            js_dir                  = t.js_dir(cIdx & t.trg_shown & trgIdx);
            frmes                   = t.frme_ts(cIdx & t.trg_shown & trgIdx);
            trg1_ts                 = t1_ts(cIdx & t.trg_shown & trgIdx);
            
            clear js_acc js dev
            for iState = 1:length(rdp_dir)
                clear js_dev
                smpl_idx            = find(frmes{iState} < trg1_ts(iState),1,'last')-nSample : find(frmes{iState} < trg1_ts(iState),1,'last');
                js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iState}(smpl_idx)),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
                js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));         	% Joystick accuracy
            end
            
            carr(iSub,iExp,iCoh)  	= snr(iCoh);                                % Coherence
            macc{iSub,iExp,iCoh}  	= js_acc;                                   % Avg accuracy
            sxc{iSub,iExp,iCoh}   	= cr.sxc(cr.coh == snr(iCoh),:);            % Smoothed cross-correlation
            mR{iSub,iExp,iCoh}   	= cr.maxR(cr.coh == snr(iCoh));             % Avg crosscorrelation coefficient
            mcc{iSub,iExp,iCoh}   	= cr.cc(cr.coh == snr(iCoh));               % Avg circular correlation coefficient
            mpk{iSub,iExp,iCoh}   	= cr.posPk(cr.coh == snr(iCoh));            % Avg crosscorrelation peak
            spk(iSub,iExp,iCoh)   	= std(cr.posPk(cr.coh == snr(iCoh)));       % Avg standard deviation of cross-correlation peak
        end
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

for iSub = 1:length(sbj_lst)
    for iCoh = 1:length(snr)
        try
            auc_acc(iSub,iCoh)      = f_auroc(macc{iSub,1,iCoh},macc{iSub,2,iCoh});
            auc_ecc(iSub,iCoh)      = f_auroc(mstr{iSub,1,iCoh},mstr{iSub,2,iCoh});
            auc_score(iSub,iCoh)    = f_auroc(trg_score{iSub,1,iCoh},trg_score{iSub,2,iCoh});
            auc_cc(iSub,iCoh)       = f_auroc(mcc{iSub,1,iCoh},mcc{iSub,2,iCoh});
            auc_xcp(iSub,iCoh)      = f_auroc(mpk{iSub,1,iCoh},mpk{iSub,2,iCoh});
            auc_xc(iSub,iCoh)       = f_auroc(mR{iSub,1,iCoh},mR{iSub,2,iCoh});
        catch
            auc_acc(iSub,iCoh)      = nan;
            auc_ecc(iSub,iCoh)      = nan;
            auc_score(iSub,iCoh)    = nan;
            auc_cc(iSub,iCoh)       = nan;
            auc_xcp(iSub,iCoh)   	= nan;
            auc_xc(iSub,iCoh)       = nan;
        end
    end
end

%% Load agent data

t_agnt                              = [];

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
    
    % Target score
    clear score score_hi score_coh
    c                           = 0;
    for iState = 1:size(t_agnt.trg_ts,1)
        for iTrg = 1:length(t_agnt.trg_ts{iState})
            c                   = c +1;
            score_cum(c)        = t_agnt.trg_score{iState}(iTrg);
            score_coh(c)        = t_agnt.rdp_coh(iState);
            score_hi(c)         = t_agnt.trg_hit{iState}(iTrg);
        end
    end
    
    tmp_score                   = score_cum(~isnan(score_cum));
    tscore(~isnan(score_cum))   = [0 diff(tmp_score)];
    tscore(tscore < 0)          = nan;
    
    for iCoh = 1:length(snr)
        clear cIdx rdp_dir t.js_dir
        
        cIdx = t_agnt.rdp_coh == snr(iCoh);
        cIdx(cellfun(@length,t_agnt.js_str) < 100) = false;
        
        % Hit rate
        nhi                         = sum(cellfun(@sum,t_agnt.trg_hit(cIdx)));
        ntrg                        = sum(cellfun(@numel,t_agnt.trg_hit(cIdx)));
        hir_agnt(iSub,iCoh)         = nhi / ntrg;
        
        % Target score
        trg_score_agnt(iSub,iCoh)	= nanmean(tscore(score_coh  == snr(iCoh) & score_hi == true));
        
        % Joystick displacement
        mstr_agnt(iSub,iCoh)        = nanmedian(cellfun(@(x) nanmedian(x(end-nSample:end)), t_agnt.js_str(cIdx)));
        
        % Joystick accuracy before first target
        t1_ts                       = cellfun(@(x) x(1), t_agnt.trg_ts);
        f1_ts                       = cellfun(@(x) x(1), t_agnt.frme_ts);
        trgIdx                      = (t1_ts-f1_ts) > 1e6;
        rdp_dir                     = t_agnt.rdp_dir(cIdx & t_agnt.trg_shown & trgIdx);
        js_dir                      = t_agnt.js_dir(cIdx & t_agnt.trg_shown & trgIdx);
        frmes                       = t_agnt.frme_ts(cIdx & t_agnt.trg_shown & trgIdx);
        trg1_ts                     = t1_ts(cIdx & t_agnt.trg_shown & trgIdx);
        
        for iState = 1:length(rdp_dir)
            clear js_dev
            smpl_idx                = find(frmes{iState} < trg1_ts(iState),1,'last')-nSample : find(frmes{iState} < trg1_ts(iState),1,'last');
            js_dev                  = rad2deg(circ_dist(deg2rad(js_dir{iState}(smpl_idx)),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
            js_acc(iState)          = nanmean(abs(1 - abs(js_dev) / 180));         	% Joystick accuracy
        end
        
        carr_agnt(iSub,iCoh)         = snr(iCoh);                                   % Coherence
        macc_agnt(iSub,iCoh)         = nanmedian(js_acc);                           % Avg accuracy
        sxc_agnt{iSub,iCoh}          = cr.sxc(cr.coh == snr(iCoh),:);               % Smoothed cross-correlation
        mR_agnt(iSub,iCoh)           = nanmean(cr.maxR(cr.coh == snr(iCoh)));       % Avg crosscorrelation coefficient
        mcc_agnt(iSub,iCoh)          = nanmean(cr.cc(cr.coh == snr(iCoh)));         % Avg circular correlation coefficient
        mpk_agnt(iSub,iCoh)          = nanmedian(cr.posPk(cr.coh == snr(iCoh)));	% Avg crosscorrelation peak
        spk_agnt(iSub,iCoh)          = std(cr.posPk(cr.coh == snr(iCoh)));          % Avg standard deviation of cross-correlation peak
    end
end

%% PLOT

f                           = figure('units','normalized','position',[0 0 .5 1]);
% height                      = [.81 .6 .35 .1];
height                    	= fliplr(linspace(.1,.75,4));
clmns                      	= linspace(.09,.83,4);
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
clm                         = .16;
row                         = .65;
hofs                        = .13;
ax0                       	= axes('Position', [clm row dim]); hold on
cmap                        = [0 0 0; gray(256)];
steps                       = .05;
bins                        = 0:steps:1;
c                           = 0;

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

colormap(cmap)

cmap_coh                    = jet(size(snr,1));

for iCoh = 1:length(snr)
    cidx                    = trg_coh == snr(iCoh);
    sc                      = scatter(trg_acc(cidx), trg_conf(cidx), 'filled');
    sc.CData            	= cmap_coh(iCoh,:);
    sc.SizeData             = 20;
    sc.MarkerFaceAlpha      = .9;
end


ax0.Position                = [clm row dim];

nBin                        = 40;
ax0h                        = axes('Position', [clm-hofs row dim(1)/5 dim(2)]); hold on
[h, edg]                    = histcounts(trg_conf,nBin);
cntr                        = edg(1:end-1) + diff(edg) ./ 2;
st                          = stairs(-h,cntr);
st.LineWidth                = lw/1.5;
st.Color                    = [0 0 0];
ax0h.YLim                   = [0 1];
ax0h.XAxis.Visible          = 'off';
ax0h.YAxis.Visible          = 'off';

ax0v                       	= axes('Position', [clm row-hofs dim(1) dim(2)/5]); hold on
[v, edg]                    = histcounts(trg_acc,nBin);
cntr                        = edg(1:end-1) + diff(edg) ./ 2;
st                          = stairs(cntr,-v);
st.LineWidth                = lw/1.5;
st.Color                    = [0 0 0];
ax0v.XLim                   = [0 1];
ax0v.XAxis.Visible          = 'off';
ax0v.YAxis.Visible          = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Performance AGNT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yof                         = .075;
xof                         = .55;
ax4                         = axes('Position', [clm+xof row-yof dim(1)/1.25 dim(2)/2]); hold on
cmap                        = cbrewer('div', 'PRGn', 4, 'PCHIP');

[pl1, ci1]                 	= plotAvgCI(macc_agnt, cmap(1,:), alp, lw);
[pl2, ci2]               	= plotAvgCI(hir_agnt, cmap(2,:), alp, lw);
[pl3, ci3]                	= plotAvgCI(mstr_agnt, cmap(3,:), alp, lw);
[pl4, ci4]                	= plotAvgCI(trg_score_agnt, cmap(4,:), alp, lw);

ax4.YLim                    = [40 100];
ax4.XLim                    = [1 size(hir_agnt,2)];
ax4.XLabel.String           = 'Coherence [%]';
ax4.YLabel.String           = 'Avg performance [%]';
ax4.XTick                   = 1:length(snr);
ax4.XTickLabel              = round(snr,2)*100;
ax4.FontSize                = lb_fs;
ax4.XTickLabelRotation      = 0;
ax4.Position                = [clm++xof row-yof dim(1)/1.25 dim(2)/2];
% ax4.XAxis.Visible           = 'off';

box off

lg                          = legend([pl1 pl2 pl3 pl4],{'Accuracy', 'Hit rate','Eccentricity' 'Score'}, 'NumColumns', 2)';
lg.Box                      = 'off';
lg.FontSize                 = lg_fs;
lg.Location                 = 'northeast';
lg.Position(1)              = .71;
lg.Position(2)              = .67;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Crosscorrelation AGNT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yof                         = .16;
ax3                         = axes('Position', [clm+xof row+yof dim(1)/1.25 dim(2)/2]); hold on
cmap                        = jet(size(snr,1));

msxc                       = cellfun(@mean,sxc_agnt,'UniformOutput',false);
for iCoh = 1:size(msxc,2)
    avg_sxc{iCoh}        	= cell2mat(msxc(:,iCoh));
end

clear pl
for iCoh = 1:size(msxc,2)
    mmat(iCoh, :)           = mean(avg_sxc{iCoh});
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
ax3.Position                = [clm+xof row+yof dim(1)/1.25 dim(2)/2];
ax3.XTickLabelRotation      = 0;

lg                          = legend(pl,lg_str,'Location','northwest','NumColumns', 2);
lg.Box                      = 'off';
lg.TextColor                = [.99 .99 .99];
lg.FontSize                 = lg_fs;
lg.Position(1)              = .17;
lg.Position(2)              = .67;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Lag difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [.15 .15];
ax6                         = axes('Position', [clmns(1) height(3) dim]); hold on

for iL = 1:length(lag)
    pl(iL)               	= plot([1 2],[lag(iL,1) lag(iL,2)]);
    pl(iL).Color          	= [.5 .5 .5 alp];
    pl(iL).LineWidth      	= lw;
end

bx                          = boxplot(lag, 'Colors', 'k');
set(bx,'MarkerEdgeColor','k')
set(bx, {'linew'},{lw})

ax6.YLabel.String           = 'Time [ms]';
ax6.XTick                   = [1 2];
ax6.XLim                    = [.5 2.5];
ax6.XTickLabel              = {'Solo','Agnt'};
ax6.FontSize                = lb_fs;
ax6.XTickLabelRotation      = 0;
ax6.Position                = [clmns(1) height(3) dim];
ax6.Box                     = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: XC coefficient difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax7                         = axes('Position', [clmns(2) height(3) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax7                         = plotData(ax7, auc_xc, {'XC coef diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax7.Position                = [clmns(2) height(3) dim];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Peak position difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax8                         = axes('Position', [clmns(1) height(4) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax8                         = plotData(ax8, auc_xcp, {'XC peak diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax8.Position                = [clmns(1) height(4) dim];

text(1.25,.825, 'Agnt > Solo', 'FontSize', lb_fs, 'Color', 'k')
text(1.25,.175, 'Agnt < Solo', 'FontSize', lb_fs, 'Color', 'k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Circular correlation difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax9                         = axes('Position', [clmns(2) height(4) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax9                         = plotData(ax9, auc_cc, {'CC coef diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax9.Position                = [clmns(2) height(4) dim];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Hit rate difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax10                      	= axes('Position', [clmns(3) height(3) dim]); hold on

ln                          = line([0 8],[0 0], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

excl                        = sum(squeeze(hir(:,1,:)),2) == 0 | sum(squeeze(hir(:,2,:)),2) == 0;
dat                         = squeeze(hir(~excl,2,:) - hir(~excl,1,:));

ax10                        = plotData(ax10, dat, 'Hit rate diff', snr, alp, lw, lb_fs, avg_mult,true);
ax10.YLim                 	= [-.15 .22];
ax10.Position             	= [clmns(3) height(3) dim];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Score difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax11                    	= axes('Position', [clmns(4) height(3) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax11                        = plotData(ax11, auc_score, {'Score diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax11.Position           	= [clmns(4) height(3) dim];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Eccentricity difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax12                        = axes('Position', [clmns(3) height(4) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax12                        = plotData(ax12, auc_ecc, {'Eccentricity diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax12.Position           	= [clmns(3) height(4) dim];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Accuracy difference AGNT - SOLO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax13                     	= axes('Position', [clmns(4) height(4) dim]); hold on

ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

ax13                        = plotData(ax13, auc_acc, {'Accuracy diff','[AUROC]'}, snr, alp, lw, lb_fs, avg_mult);
ax13.Position             	= [clmns(4) height(4) dim];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax0                         = axes('Position',[0 0 1 1],'Visible','off');
lofs                        = .2;
text(0.03,height(1)+lofs, 'A', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(0.6,height(1)+lofs, 'B', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(0.6,height(2)+lofs, 'C', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')

lofs                        = .175;
text(0.002,height(3)+lofs, 'D', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(0.49,height(3)+lofs, 'E', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')

text(0.2,.75, 'Coherence', 'Parent', ax0, 'FontSize', lg_fs, 'Color', [.99 .99 .99])

print(f, '/Users/fschneider/ownCloud/Documents/Publications/CPR_psychophysics/Figures/FIG3/FIG3', '-r400', '-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tmp, tc] = extractTrials(tbl, tmp, tc)
for iTrl = 1:tbl.t.trl_no(end)
    tidx                    = tbl.t.trl_no == iTrl;             % Trial index
    tt                      = tbl.t(tidx,:);                    % Relevant stimulus cycles
    
    % Initiate cell
    tc                      = tc+1;                             % Trial counter
    tmp.frme_ts{tc}         = [];
    tmp.rdp_dir{tc}         = [];
    tmp.rdp_coh{tc}     	= [];
    tmp.js_dir{tc}      	= [];
    tmp.js_str{tc}      	= [];
    tmp.refresh{tc}         = [];
    
    % Extract and add data
    for iState = 1:size(tt,1)
        tmp.frme_ts{tc}  	= [tmp.frme_ts{tc} tt.frme_ts{iState}];
        tmp.rdp_dir{tc}  	= [tmp.rdp_dir{tc} repmat(tt.rdp_dir(iState),1,length(tt.frme_ts{iState}))];
        tmp.rdp_coh{tc}  	= [tmp.rdp_coh{tc} repmat(tt.rdp_coh(iState),1,length(tt.frme_ts{iState}))];
        tmp.js_dir{tc}  	= [tmp.js_dir{tc} tt.js_dir{iState}];
        tmp.js_str{tc}  	= [tmp.js_str{tc} tt.js_str{iState}];
        tmp.refresh{tc}     = [tmp.refresh{tc} median(diff(tmp.frme_ts{tc}))];
    end
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

if flag == false
    pt                      = patch([0 8 8 0], [.5 .5 1 1], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
    pt                      = patch([0 8 8 0], [0 0 .5 .5], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
else
    pt                      = patch([0 8 8 0], [0 0 .5 .5], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
    pt                      = patch([0 8 8 0], [-.5 -.5 0 0], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
end

for iL = 1:size(dat,1)
    pl(iL)               	= plot(dat(iL,:));
    pl(iL).Color          	= [.5 .5 .5 alp];
    pl(iL).LineWidth      	= lw/2;
end

pm                          = plot(nanmedian(dat),'LineWidth', lw*avg_mult, 'Color', [0 0 0]);

ax.XLim                     = [1 size(dat,2)];
ax.XLabel.String            = 'Coherence [%]';
ax.YLabel.String            = ylab;
ax.XTick                    = [1 4 length(snr)];
ax.XTickLabel               = round([snr(1) snr(4) snr(end)],2)*100;
ax.FontSize                 = lb_fs;
ax.YLim                     = [.1 .9];
ax.XTickLabelRotation       = 0;

box off
end
