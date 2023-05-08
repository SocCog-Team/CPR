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

sbj_lst(cellfun(@(x) strcmp(x, 'IrS'),sbj_lst)) = [];

% Extract file names of solo experiments
for iSubj = 1:30%length(sbj_lst)
    load([pth sbj_lst{iSubj} '/summary/' sbj_lst{iSubj} '_RT.mat'])
    mrt(iSubj)              = median(rt.dat);
    
    load([pth sbj_lst{iSubj} '/summary/' sbj_lst{iSubj} '_SNRfit.mat'])
    psy_df                     = abs(psy_func.model - .5);
    tresh_smpl(iSubj)       = find(psy_df == min(psy_df));
    thresh_snr(iSubj)       = psy_func.model_snr(tresh_smpl(iSubj));
    
    data_pth                = [pth sbj_lst{iSubj} '/summary/'];
    
    if isdir(data_pth)
        cd(data_pth)
        mat_files        	= dir('*.mat');
        
        for iFile = 1:length(mat_files)
            if contains(mat_files(iFile).name,'CPRsolo')
                c        	= c+1;
                fname{c}   	= mat_files(iFile).name;
            end
        end
    end
end

%% Extract data

nSample                     = 30;
s_cnt                       = 0;
nLag                        = 150;

for iSub = 1:30%length(sbj_lst)
    
    fname_crop                      = cellfun(@(x) x(1:12), fname, 'UniformOutput', false);
    fidx                            = find(cellfun(@(x) contains(x,lower(sbj_lst{iSub})),fname_crop));
    cc                              = 0;
    t                               = [];
    all                             = [];
    tmp                             = {};
    
    if isempty(fidx)
        lag(iSub) = nan;
        continue
    end
    
    s_cnt                           = s_cnt + 1;
    
    for iBlock = 1:length(fidx)
        clear tbl
        tbl = load([pth sbj_lst{iSub} '/summary/' fname{fidx(iBlock)}]);
        
        if sum(unique(tbl.t.rdp_coh) < .1) > 2
            continue
        end
        
        t                       = [t; tbl.t];
        [tmp, tc]               = extractTrials(tbl, tmp, cc);
      
    end
    
    snr                         = unique(t.rdp_coh);
    carr(s_cnt,:)               = nan(1,13);
    id{s_cnt}                   = sbj_lst{iSub};
    [cr,ps]                     = CPR_correlation_analysis_WIP(tmp, nLag, false);
    lag(iSub)                   = median(cr.lag);
    
    % Target score
    clear score score_hi score_coh
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
    
    for iCoh = 1:length(snr)
        clear cIdx rdp_dir t.js_dir
        
        cIdx = t.rdp_coh == snr(iCoh);
        cIdx(cellfun(@length,t.js_str) < 100) = false;
        
        % Hit rate
        nhi                     = sum(cellfun(@sum,t.trg_hit(cIdx)));
        ntrg                    = sum(cellfun(@numel,t.trg_hit(cIdx)));
        hir(s_cnt,iCoh)         = nhi / ntrg;
        
        % Target score
        trg_score(s_cnt,iCoh) 	= nanmean(tscore(score_coh  == snr(iCoh) & score_hi == true));
        
        % Joystick displacement
        mstr(s_cnt, iCoh)       = nanmedian(cellfun(@(x) nanmedian(x(end-nSample:end)), t.js_str(cIdx)));
        
        % Joystick accuracy before first target
        t1_ts                   = cellfun(@(x) x(1), t.trg_ts);
        f1_ts                   = cellfun(@(x) x(1), t.frme_ts);
        trgIdx                  = (t1_ts-f1_ts) > 1e6;
        rdp_dir                 = t.rdp_dir(cIdx & t.trg_shown & trgIdx);
        js_dir                  = t.js_dir(cIdx & t.trg_shown & trgIdx);
        frmes                   = t.frme_ts(cIdx & t.trg_shown & trgIdx);
        trg1_ts                 = t1_ts(cIdx & t.trg_shown & trgIdx);
        
        for iState = 1:length(rdp_dir)
            clear js_dev
            smpl_idx            = find(frmes{iState} < trg1_ts(iState),1,'last')-nSample : find(frmes{iState} < trg1_ts(iState),1,'last');
            js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iState}(smpl_idx)),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
            js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));                           % Joystick accuracy
        end
        
        %         rdp_dir                 = t.rdp_dir(cIdx & t.trg_shown);
        %         js_dir                  = t.js_dir(cIdx & t.trg_shown);
        %         for iState = 1:length(rdp_dir)
        %             clear js_dev
        %             js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iState}(end-nSample:end)),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
        %             js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));                           % Joystick accuracy
        %         end
        
        carr(s_cnt,snr == snr(iCoh))   	= snr(iCoh);
        macc(s_cnt,snr == snr(iCoh))    = nanmedian(js_acc);
        sxc{s_cnt,snr == snr(iCoh)}   	= cr.sxc(cr.coh == snr(iCoh),:);
        mR(s_cnt,snr == snr(iCoh))     	= nanmean(cr.maxR(cr.coh == snr(iCoh)));
        mcc(s_cnt,snr == snr(iCoh))   	= nanmean(cr.cc(cr.coh == snr(iCoh)));
        mpk(s_cnt,snr == snr(iCoh))   	= nanmedian(cr.posPk(cr.coh == snr(iCoh)));
        spk(s_cnt,snr == snr(iCoh))    	= std(cr.posPk(cr.coh == snr(iCoh)));
    end
end

%% PLOT

f                           = figure('units','normalized','position',[0 0 .5 1]);
height                      = [.78 .59 .325 .135 .06];
lb_fs                       = 14;
lg_fs                       = 10;
lw                          = 3;
nLag                        = 150;
frme_ms                     = 1000/120;
alp                         = .35;
avg_mult                    = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT A: Example subject joystick response %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [0.2 0.2]*1.33;
clm                         = .12;
row                         = .69;
hofs                        = .115;
ax0                       	= axes('Position', [clm row dim]); hold on
cmap                        = [0 0 0; gray(256)];
steps                       = .05;
bins                        = 0:steps:1;
c                           = 0;
tbl                         = load([pth 'AnM/summary/20220629_anm_CPRsolo_block2_tbl.mat']);
exmpl_id                    = find(strcmp(sbj_lst, 'AnM'));
exmpl_col                   = [.5 0 .5 .5];
t                           = tbl.t;

% Extracte experimental data
clear trg_acc trg_conf trg_coh trg_hit
for iState = 1:size(t,1)
    if t.trg_shown(iState) == false
        continue
    end
    
    for iTarget = 1:length(t.trg_ts{iState})
        c                   = c+1;
        trg_acc(c)          = t.trg_acc{iState}(iTarget);
        trg_conf(c)         = t.trg_ecc{iState}(iTarget);
        trg_coh(c)          = t.rdp_coh(iState);
        trg_hit(c)          = t.trg_hit{iState}(iTarget);
        
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
ax0.XLabel.String           = 'Accuracy';
ax0.YLabel.String           = 'Eccentricity';
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
    sc.SizeData             = 15;
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
%%% SUBPLOT: HIT RATE & Score Raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [0.175 0.175];
ax4                         = axes('Position', [.575 height(1) dim]); hold on
[ax4,pl]                    = plotData(ax4,hir,false,lw,alp,exmpl_id,exmpl_col,avg_mult);
ax4.YLim                    = [35 73];
ax4.XLim                    = [1 size(hir,2)];
ax4.XLabel.String           = 'Coherence [%]';
ax4.YLabel.String           = 'Hit rate [%]';
ax4.XTick                   = 1:length(snr);
ax4.XTickLabel              = round(snr,2)*100;
ax4.FontSize                = lb_fs;
ax4.XTickLabelRotation      = 0;
ax4.XAxis.Visible           = 'off';


ax5                         = axes('Position', [.575 height(2) dim]); hold on
[ax5,pl]                    = plotData(ax5,trg_score,false,lw,alp,exmpl_id,exmpl_col,avg_mult);
ax5.YLim                    = [20 92];
ax5.XLim                    = [1 size(trg_score,2)];
ax5.XLabel.String           = 'Coherence [%]';
ax5.YLabel.String           = 'Score [%]';
ax5.XTick                   = 1:length(snr);
ax5.YTick                   = [30 60 90];
ax5.XTickLabel              = round(snr,2)*100;
ax5.FontSize                = lb_fs;
ax5.XTickLabelRotation      = 0;
% ax5.XAxis.Visible           = 'off';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: HIT RATE & Score Normalised %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax44                        = axes('Position', [.815 height(1) dim]); hold on
ln44                    	= line([0 8],[1 1],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
[ax44,pl]                    = plotData(ax44,hir,true,lw,alp,exmpl_id,exmpl_col,avg_mult);
ax44.XLim                    = [1 size(hir,2)];
ax44.XLabel.String           = 'Coherence [%]';
ax44.YLabel.String           = 'Hit rate [norm]';
ax44.XTick                   = 1:length(snr);
ax44.XTickLabel              = round(snr,2)*100;
ax44.YTick                   = [.8 1 1.2];
ax44.FontSize                = lb_fs;
ax44.XTickLabelRotation      = 0;
ax44.XAxis.Visible           = 'off';

ax55                        = axes('Position', [.815 height(2) dim]); hold on
ln55                     	= line([0 8],[1 1],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
[ax55,pl]                   = plotData(ax55,trg_score,true,lw,alp,exmpl_id,exmpl_col,avg_mult);
ax55.XLim                    = [1 size(trg_score,2)];
ax55.XLabel.String           = 'Coherence [%]';
ax55.YLabel.String           = 'Score [norm]';
ax55.XTick                   = 1:length(snr);
ax55.YTick                   = [.6 1 1.4];
ax55.XTickLabel              = round(snr,2)*100;
ax55.FontSize                = lb_fs;
ax55.XTickLabelRotation      = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Crosscorrelation Pooled %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [.175 .175];
yof                         = .06;
ax3                         = axes('Position', [.075 height(3)+yof dim]); hold on
cmap                        = jet(size(snr,1));

msxc                        = cellfun(@mean,sxc,'UniformOutput',false);
for iCoh = 1:size(msxc,2)
    avg_sxc{iCoh}        	= cell2mat(msxc(:,iCoh));
end

for iCoh = 1:size(msxc,2)
    pl(iCoh)               	= plot(mean(avg_sxc{iCoh}),'LineWidth', lw, 'Color', cmap(iCoh,:));
    lg_str{iCoh}            = num2str(round(snr(iCoh)*100));
end

ax3.YLabel.String           = 'XCorr Coeff';
ax3.XLabel.String           = 'Lag [ms]';
ax3.XLim                    = [0 301];
ax3.YLim                    = [0 max(mean(avg_sxc{iCoh}))];
ax3.XTick                   = [0 150 300];
ax3.FontSize              	= lb_fs;
ax3.XTickLabel              = round((cellfun(@str2num, ax3.XTickLabel)-nLag) * frme_ms);
ax3.Position                = [.075 height(3)+yof dim(1)*1.33 dim(2)];

ln                          = line([150 150],[0 .2]);
ln.LineStyle                = ':';
ln.LineWidth                = lw;
ln.Color                    = [0 0 0];

% lg                          = legend(pl,lg_str,'Location','northwest','NumColumns',2);
lg                          = legend(pl,lg_str,'Location','northwest');
lg.Box                      = 'off';
lg.FontSize                 = lg_fs;
lg.Position(1)              = .09;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Decision time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [.1 .175];
xof                         = .38;
ax6                         = axes('Position', [xof height(4) dim]); hold on

for iS = 1:length(lag)
    pl(iS)               	= plot([1 2],[mrt(iS) lag(iS)]);
    pl(iS).Color          	= [.5 .5 .5 alp];
    pl(iS).LineWidth      	= lw/2;
    if iS == exmpl_id
        %         pl(iS).LineStyle 	= ':';
        pl(iS).LineWidth 	= lw;
        pl(iS).Color    	= exmpl_col;
    end
end

bx                          = boxplot([mrt' lag'], 'Colors', 'k');
% set(bx(end,:),'Visible','off')
set(bx,'MarkerEdgeColor','k')
set(bx, {'linew'},{lw/2})

ax6.YLabel.String           = 'Time [ms]';
ax6.XTick                   = [1 2];
% ax6.YTick                   = [.6 .8 1 1.2];
ax6.XLim                    = [.5 2.5];
ax6.XTickLabel              = {'RT','CPR'};
ax6.FontSize                = lb_fs;
ax6.XTickLabelRotation      = 0;
ax6.Position                = [xof height(4) dim];
ax6.Box                     = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Crosscorrelation Peak Times %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim                         = [.175 .175];
dim                         = [dim(1)*1.33 dim(2)*.5];
xof                         = .075;
yof                         = .09;

%%% Median
ax33                        = axes('Position', [xof height(4) dim]); hold on
bx                          = boxplot(abs(mpk*frme_ms), 'Color' ,'k');

ax33.YLabel.String          = 'Lag [ms]';
% ax33.YLabel.Position(1:2)   =  [-1.25 28 ];
ax33.XLabel.String          = 'Coherence [%]';
ax33.XTick                  = 1:length(snr);
ax33.YLim                   = [300 1100];
ax33.YTick                  = [400 600 800 1000];
ax33.XTickLabel             = round(snr,2)*100;
ax33.FontSize              	= lb_fs;
ax33.XTickLabelRotation     = 0;
ax33.Box                    = 'off';
ax33.Position               = [xof height(4)+yof dim];
ax33.XAxis.Visible           = 'off';
set(bx, {'linew'},{lw/2})
set(bx,'MarkerEdgeColor','k')

%%% Standard deviation
ax34                        = axes('Position', [xof height(4) dim]); hold on
bx                          = boxplot(spk, 'Color' ,'k');

ax34.YLabel.String          = 'Std';
ax34.XLabel.String          = 'Coherence [%]';
ax34.XTick                  = 1:length(snr);
ax34.XTickLabel             = round(snr,2)*100;
ax34.FontSize              	= lb_fs;
ax34.XTickLabelRotation     = 0;
ax34.Box                    = 'off';
ax34.Position               = [xof height(4) dim];
set(bx, {'linew'},{lw/2})
set(bx,'MarkerEdgeColor','k')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Peak-to-average ratio %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [0.1 0.175];
xof                         = .38;
yof                         = .06;
ax7                         = axes('Position', [xof height(3)+yof dim]); hold on
msxc                     	= cellfun(@mean,sxc,'UniformOutput',false);

for iSub = 1:size(sxc,1)
    for iCoh = 1:size(sxc,2)
        x                	= msxc{iSub,iCoh};
        par(iSub,iCoh)   	= max(x) / mean(x);
    end
end

[ax7,pl]                   = plotData(ax7,par,false,lw,alp,exmpl_id,exmpl_col,avg_mult);
ax7.XLim                    = [1 7];
ax7.XTick                   = [1 4 7];
ax7.XTickLabel              = round([snr(1) snr(4) snr(7)],2)*100;
ax7.XLabel.String           = {'Coherence'; '[%]'};
ax7.YLabel.String           = {'Peak-Average ratio'};
ax7.FontSize                = lb_fs;
ax7.XTickLabelRotation      = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Accuracy & Eccentricity Raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim                         = [.175 .175];
ax1                         = axes('Position', [.575 height(4) dim]); hold on
[ax1,pl]                   = plotData(ax1,mstr,false,lw,alp,exmpl_id,exmpl_col,avg_mult);
% ax1.YLim                    = [.5 1.5];
ax1.XLim                    = [1 size(mstr,2)];
ax1.XLabel.String           = 'Coherence [%]';
ax1.YLabel.String           = 'Ecc [%]';
ax1.XTick                   = 1:length(snr);
ax1.XTickLabel              = round(snr,2)*100;
ax1.FontSize                = lb_fs;
ax1.XTickLabelRotation      = 0;
% ax1.XAxis.Visible           = 'off';

ax2                         = axes('Position', [.575 height(3) dim]); hold on
[ax2,pl]                   = plotData(ax2,macc,false,lw,alp,exmpl_id,exmpl_col,avg_mult);
% ax2.YLim                    = [.5 1.5];
ax2.XLim                    = [1 size(macc,2)];
ax2.XLabel.String           = 'Coherence [%]';
ax2.YLabel.String           = 'Acc [%]';
ax2.XTick                   = 1:length(snr);
% ax2.YTick                   = [30 60 90];
ax2.XTickLabel              = round(snr,2)*100;
ax2.FontSize                = lb_fs;
ax2.XTickLabelRotation      = 0;
ax2.XAxis.Visible           = 'off';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Accuracy & Eccentricity Normalised %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax11                        = axes('Position', [.815 height(4) dim]); hold on
ln11                        = line([0 8],[1 1],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
[ax11,pl]                   = plotData(ax11,mstr,true,lw,alp,exmpl_id,exmpl_col,avg_mult);
ax11.XLim                    = [1 size(mstr,2)];
ax11.XLabel.String           = 'Coherence [%]';
ax11.YLabel.String           = 'Ecc [norm]';
ax11.XTick                   = 1:length(snr);
ax11.XTickLabel              = round(snr,2)*100;
ax11.YTick                   = [.6 1 1.4];
ax11.YLim                    = [.5 1.42];
ax11.FontSize                = lb_fs;
ax11.XTickLabelRotation      = 0;
% ax11.XAxis.Visible           = 'off';

ax22                        = axes('Position', [.815 height(3) dim]); hold on
ln22                     	= line([0 8],[1 1],'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
[ax22,pl]                   = plotData(ax22,macc,true,lw,alp,exmpl_id,exmpl_col,avg_mult);
ax22.XLim                    = [1 size(macc,2)];
ax22.XLabel.String           = 'Coherence [%]';
ax22.YLabel.String           = 'Acc [norm]';
ax22.XTick                   = 1:length(snr);
ax22.YTick                   = [.6 .8 1 1.2];
ax22.YLim                    = [.6 1.22];
ax22.XTickLabel              = round(snr,2)*100;
ax22.FontSize                = lb_fs;
ax22.XTickLabelRotation      = 0;
ax22.XAxis.Visible           = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lofs                        = .185;
ax0                         = axes('Position',[0 0 1 1],'Visible','off');
text(0.002,height(1)+lofs, 'A', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(0.49,height(1)+lofs, 'B', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(0.49,height(2)+lofs, 'C', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(0.002,height(3)+lofs+.06, 'D', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(0.002,height(4)+lofs, 'E', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(0.3,height(4)+lofs, 'F', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(0.49,height(3)+lofs, 'G', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')
text(0.49,height(4)+lofs, 'H', 'Parent', ax0, 'FontSize', 30, 'Color', 'k')

text(0.09,.57, 'Coherence', 'Parent', ax0, 'FontSize', lg_fs, 'Color', 'k')

print(f, '/Users/fschneider/ownCloud/Documents/Publications/CPR_psychophysics/Figures/FIG2/FIG2', '-r400', '-dpng');

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

function [ax,pl] = plotData(ax,dat,norm_flag,lw,alp,exmpl_id,exmpl_col,avg_mult)
for iL = 1:size(dat,1)
    if norm_flag
        pl(iL)                  = plot((dat(iL,:)./mean(dat(iL,:))), 'LineWidth', lw/2, 'Color', [.5 .5 .5 alp]);
    else
        pl(iL)                  = plot(dat(iL,:)*100, 'LineWidth', lw/2, 'Color', [.5 .5 .5 alp]);    
    end
    
    if iL == exmpl_id
        pl(iL).LineWidth 	= lw/1.25;
        pl(iL).Color      	= exmpl_col;
    end
end

if norm_flag
    plot(mean(dat./mean(dat,2)),'LineWidth', lw*avg_mult, 'Color', [0 0 0])
else
    plot(mean(dat)*100,'LineWidth', lw*avg_mult, 'Color', [0 0 0])
end
end