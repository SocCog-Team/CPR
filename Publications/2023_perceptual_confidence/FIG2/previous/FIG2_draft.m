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
pth                         = '/Users/fschneider/Documents/CPR_psychophysics/';
x                           = readtable([pth 'Subjects_summary.xlsx']);
sbj_lst                     = x.Abbreviation;
sbj_lst(cellfun(@isempty,sbj_lst)) = [];
c                           = 0;

sbj_lst(cellfun(@(x) strcmp(x, 'IrS'),sbj_lst)) = [];
    
% Extract file names of solo experiments 
for iSubj = 1:20%length(sbj_lst)
    load([pth sbj_lst{iSubj} '/summary/' sbj_lst{iSubj} '_RT.mat'])
    mrt(iSubj)              = mean(rt.dat);
    
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

nSample                     = 30;
s_cnt                       = 0;
nLag                        = 150;
                          
for iSub = 1:20%length(sbj_lst)
    
    fname_crop                      = cellfun(@(x) x(1:12), fname, 'UniformOutput', false);
    fidx                            = find(cellfun(@(x) contains(x,lower(sbj_lst{iSub})),fname_crop));
    cc                              = 0;
    t                               = [];
    all                             = [];
    
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
        
        t = [t; tbl.t];
        
        % Extract trials for correlation analysis
        for iTrl = 1:tbl.t.trl_no(end)
            tidx                    = tbl.t.trl_no == iTrl;
            tt                      = tbl.t(tidx,:);
            
            cc                      = cc+1;
            tmp.frme_ts{cc}         = [];
            tmp.rdp_dir{cc}         = [];
            tmp.rdp_coh{cc}     	= [];
            tmp.js_dir{cc}      	= [];
            tmp.js_str{cc}      	= [];
            tmp.refresh{cc}         = [];
            
            for iState = 1:size(tt,1)
                tmp.frme_ts{cc}  	= [tmp.frme_ts{cc} tt.frme_ts{iState}];
                tmp.rdp_dir{cc}  	= [tmp.rdp_dir{cc} repmat(tt.rdp_dir(iState),1,length(tt.frme_ts{iState}))];
                tmp.rdp_coh{cc}  	= [tmp.rdp_coh{cc} repmat(tt.rdp_coh(iState),1,length(tt.frme_ts{iState}))];
                tmp.js_dir{cc}  	= [tmp.js_dir{cc} tt.js_dir{iState}];
                tmp.js_str{cc}  	= [tmp.js_str{cc} tt.js_str{iState}];
                tmp.refresh{cc}     = [tmp.refresh{cc} median(diff(tmp.frme_ts{cc}))];
            end
        end
        
    end
    
    snr                         = unique(t.rdp_coh);
    
    carr(s_cnt,:)               = nan(1,13); 
    id{s_cnt}                   = sbj_lst{iSub};
    nLag                        = 150;
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
        mpk(s_cnt,snr == snr(iCoh))   	= nanmean(cr.posPk(cr.coh == snr(iCoh)));
        spk(s_cnt,snr == snr(iCoh))    	= std(cr.posPk(cr.coh == snr(iCoh)));
    end
end

%% Plot

f                   = figure('units','normalized','outerposition',[0 0 1 1]);
dim                 = [0.24 0.24];
clm                 = linspace(.06,.7,3);
height              = [.75 .4 .097];
fs                  = 10;
lb_fs               = 18;
lg_fs               = 14;
lw                  = 3;
frme_ms             = 1000/120;

%%% Crosscorrelation - Raw
ax1                 = axes('Position', [clm(1) height(1) dim]);
ID                  = 3;
x                   = carr(ID,~isnan(carr(ID,:)));
% cmap                = flipud(cbrewer('div', 'RdYlBu', size(x,2), 'PCHIP'));
cmap                = jet(size(x,2));

hold on
tmpp                = cellfun(@mean,sxc(ID,:),'UniformOutput',false);
for iCoh = 1:size(x,2)
    pl(iCoh)          	= plot(tmpp{iCoh});
    pl(iCoh).LineWidth	= lw;
    pl(iCoh).Color    	= [cmap(iCoh,:)];
    lg_str{iCoh}        = num2str(round(snr(iCoh),2));
end

ax1.YLabel.String   = 'XCorr Coeff';
ax1.XLabel.String   = 'Lag [ms]';
ax1.XLim            = [0 301];
ax1.YLim            = [0 max(tmpp{end})];
ax1.FontSize        = lb_fs;
ax1.XTickLabel      = round((cellfun(@str2num, ax1.XTickLabel)-nLag) * frme_ms);

ln = line([150 150],[0 .2]);
ln.LineStyle = ':';
ln.LineWidth = 2;
ln.Color = [0 0 0];

lg                	= legend(pl,lg_str,'Location','northwest','NumColumns',2);
lg.Position(1:2)    = [clm(1)+.01 height(1)+.115];
lg.Box              = 'off';
lg.FontSize         = lg_fs;


%%% Peak-to-average ratio
ax3                 = axes('Position', [clm(3) height(1) dim]);
msxc                = cellfun(@mean,sxc,'UniformOutput',false);

for iSub = 1:size(sxc,1)
    for iCoh = 1:size(sxc,2)
        x               = msxc{iSub,iCoh};
        par(iSub,iCoh)  = max(x) / mean(x);
    end
end

[~,sidx]            = sort(mean(par,2));

im                  = imagesc(par(sidx,:));
ax3.XTick         	= 1:length(snr);
ax3.XTickLabel    	= round(snr,2)*100;
ax3.XLabel.String   = 'Coherence [%]';
ax3.YLabel.String   = 'Subject';
ax3.FontSize        = lb_fs;
% ax3.XTickLabelRotation = 15;
cb                  = colorbar;
cb.Label.String     = {'Peak-to-Mean Ratio';'"Response reliability"'};
cb.FontSize         = lg_fs;
colormap(gray(256))
set(gca,'YDir','normal')

yyaxis right
% lp                          = plot(mean(par(sidx,:) ./ mean(par(sidx,:),2)));
lp                          = plot(mean(par(sidx,:)));
lp.LineWidth                = 2;
lp.Color                    = [1 0 0];
ax3.YAxis(2).Color          = [1 0 0];
ax3.YAxis(2).Limits         = [1 6];

%%% Average lag
ax2                         = axes('Position', [clm(2) height(1) dim]); hold on
cl                          = cbrewer('qual', 'Paired', size(sbj_lst,1), 'PCHIP');

vl                          = violinplot([mrt lag],[ones(1,length(mrt)) ones(1,length(lag))+1]);
for v = 1:2
vl(v).ShowData = 0;
vl(v).WhiskerPlot.Color     = 'none';
vl(v).BoxPlot.FaceColor     = 'none';
vl(v).BoxPlot.EdgeColor     = 'none';
vl(v).MedianPlot.Marker     = 'none';
vl(v).ViolinColor           = [.5 .5 .5];
end

cl(isnan(lag'),:)          	= [];
mrt(isnan(lag))             = [];
lag(isnan(lag))             = [];

for iSubj = 1:size(lag,2)    
    plot([1 2],[mrt(iSubj) lag(iSubj)], 'Color', [cl(iSubj,:) .5], 'LineWidth',2)
    scatter(1, mrt(iSubj),'filled', 'CData', cl(iSubj,:), 'SizeData', 100, 'MarkerFaceAlpha', .75);
    scatter(2, lag(iSubj),'filled', 'CData', cl(iSubj,:), 'SizeData', 100, 'MarkerFaceAlpha', .75)
end

ax2.XTick                   = [1 2];
ax2.XTickLabel              = {'Reaction Time', 'CPR_solo'};
ax2.TickLabelInterpreter    = 'none';
ax2.XLim                    = [.5 2.5];
ax2.YLim                    = [200 800];
ax2.YLabel.String           = 'Time [ms]';
ax2.FontSize                = lb_fs;

%%% Hit rate
cl                          = cbrewer('qual', 'Paired', size(id,2), 'PCHIP');
ax4                         = axes('Position', [clm(2) height(3) dim]); hold on
im                          = imagesc(hir(sidx,:) ./ mean(hir(sidx,:),2));
% im                          = imagesc(hir(sidx,:));

ax4.YLim                    = [.5 size(hir,1)+.5];
ax4.XLim                    = [.5 size(hir,2)+.5];
ax4.XLabel.String           = 'Coherence [%]';
ax4.YLabel.String           = 'Subject';
% ax4.YTick                   = 1:length(id);
% ax4.YTickLabel              = id(sidx);
% ax4.YTick                   = [];
ax4.XTick                   = 1:length(snr);
ax4.XTickLabel              = round(snr,2)*100;
ax4.FontSize                = fs;
ax4.FontSize                = lb_fs;
% ax4.XTickLabelRotation      = 15;
% ax4.XLabel.FontSize         = lb_fs;
% ax4.YLabel.FontSize         = lb_fs;
cb                          = colorbar;
cb.Label.String             = 'Hit rate [norm]';
cb.FontSize                 = lg_fs;
colormap(gray(256))

caxis([.8 1.2])

yyaxis right
lp                          = plot(mean(hir(sidx,:) ./ mean(hir(sidx,:),2)));
% lp                          = plot(mean(hir(sidx,:)));
lp.LineWidth                = 2;
lp.Color                    = [1 0 0];
ax4.YAxis(2).Color          = [1 0 0];

%%% Target score
ax5                         = axes('Position', [clm(3) height(3) dim]); hold on
im                          = imagesc(trg_score(sidx,:) ./ mean(trg_score(sidx,:),2));
% im                          = imagesc(trg_score(sidx,:));

ax5.YLim                    = [.5 size(trg_score,1)+.5];
ax5.XLim                    = [.5 size(trg_score,2)+.5];
ax5.XLabel.String           = 'Coherence [%]';
ax5.YLabel.String           = 'Subject';
% ax5.YTick                   = 1:length(id);
% ax5.YTickLabel              = id(sidx);
% ax5.YTick                   = [];
ax5.XTick                   = 1:length(snr);
ax5.XTickLabel              = round(snr,2)*100;
ax5.FontSize                = lb_fs;
% ax5.XTickLabelRotation      = 15;
% ax5.XLabel.FontSize         = lg_fs;
% ax5.YLabel.FontSize         = lg_fs;
cb                          = colorbar;
cb.Label.String             = 'Target Score [norm]';
cb.FontSize                 = lg_fs;
colormap(gray(256))

caxis([.8 1.2])

yyaxis right
lp                          = plot(mean(trg_score(sidx,:) ./ mean(trg_score(sidx,:),2)));
% lp                          = plot(mean(trg_score(sidx,:)));
lp.LineWidth                = 2;
lp.Color                    = [1 0 0];
ax5.YAxis(2).Color          = [1 0 0];

%%% Accuracy
ax6                         = axes('Position', [clm(2) height(2) dim]); hold on
im                          = imagesc(macc(sidx,:) ./ mean(macc(sidx,:),2));
% im                          = imagesc(macc(sidx,:));

ax6.YLim                    = [.5 size(trg_score,1)+.5];
ax6.XLim                    = [.5 size(trg_score,2)+.5];
% ax6.XLabel.String           = 'Coherence [%]';
ax6.YLabel.String           = 'Subject';
% ax6.YTick                   = 1:length(id);
% ax6.YTickLabel              = id(sidx);
ax6.XTick                   = 1:length(snr);
ax6.XTickLabel              = round(snr,2)*100;
ax6.FontSize                = lb_fs;
% ax6.XTickLabelRotation      = 15;
% ax6.XLabel.FontSize     	= lb_fs;
% ax6.YLabel.FontSize         = lb_fs;
cb                          = colorbar;
cb.Label.String             = 'Response Accuracy [norm]';
cb.FontSize                 = lg_fs;
colormap(gray(256))

caxis([.9 1.1])

yyaxis right
lp                          = plot(mean(macc(sidx,:) ./ mean(macc(sidx,:),2)));
% lp                          = plot(mean(macc(sidx,:)));
lp.LineWidth                = 2;
lp.Color                    = [1 0 0];
ax6.YAxis(2).Color          = [1 0 0];

%%% Displacement
ax7                         = axes('Position', [clm(3) height(2) dim]); hold on
im                          = imagesc(mstr(sidx,:) ./ mean(mstr(sidx,:),2));
% im                          = imagesc(mstr(sidx,:));

ax7.YLim                    = [.5 size(trg_score,1)+.5];
ax7.XLim                    = [.5 size(trg_score,2)+.5];
% ax7.XLabel.String           = 'Coherence [%]';
ax7.YLabel.String           = 'Subject';
% ax7.YTick                   = 1:length(id);
% ax7.YTickLabel              = id(sidx);
ax7.XTick                   = 1:length(snr);
ax7.XTickLabel              = round(snr,2)*100;
ax7.FontSize                = lb_fs;
% ax7.XTickLabelRotation     	= 15;
% ax7.XLabel.FontSize         = lb_fs;
% ax7.YLabel.FontSize         = lb_fs;

cb                          = colorbar;
cb.Label.String             = 'Response Eccentricity [norm]';
cb.FontSize                 = lg_fs;
colormap(gray(256))

caxis([.8 1.2])

yyaxis right
lp                          = plot(mean(mstr(sidx,:) ./ mean(mstr(sidx,:),2)));
% lp                          = plot(mean(mstr(sidx,:)));
lp.LineWidth                = 2;
lp.Color                    = [1 0 0];
ax7.YAxis(2).Color          = [1 0 0];

% Reward map
ax8                       	= axes('Position', [clm(1) height(2) dim]); hold on
cmap                        = [0 0 0; gray(256)];
steps                       = .05;
bins                        = 0:steps:1;
c                           = 0;
tbl                         = load([pth 'AnM/summary/20220629_anm_CPRsolo_block2_tbl.mat']);
% tbl                         = load([pth 'HeL/summary/20220902_hel_CPRsolo_block2_tbl.mat']);
t                           = tbl.t;

% Extracte experimental data
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
acc              	= 0:.001:1;
conf                = 0:.001:1;
rew                 = acc' .* conf;

% Determine arc width for each confidence level
for j = 1:length(conf)
    arc(j)          = 180 - (180 * conf(j));
end

% Cap arc width at target width (2dva == 12.7587deg at chosen position)
aidx                = arc < 12.7587;
arc(aidx)           = 12.7587;

% For each confidence level, calculate minimum accuracy required to hit
% the target at given arc width - normalised values
hit_width_acc       = 1 - ((arc/2) / 180);
hit_width_acc(aidx) = 1 - (12.7587/2)/180; % arc width fixed

% Remove position that cannot yield reward from reward matrix
for iAcc = 1:length(acc)
    indx            = conf < hit_width_acc(iAcc);
    rew(iAcc,indx)  = nan;
end

% Plot reward matrix
hold on
ax                  = gca;
im                  = imagesc(acc,conf,rew);
ax.XLabel.String    = 'Accuracy';
ax.YLabel.String    = 'Eccentricity';
ax.FontSize         = lb_fs;
ax.XLim             = [0 1];
ax.YLim             = [0 1];
cb                  = colorbar;
cb.Label.String     = '% Reward';

ax.XLim             = [0 1];
ax.XTick            = [0:.2:1];
ax.YTick            = [0:.2:1];
colormap(cmap)

% cmap                = flipud(cbrewer('div', 'RdYlBu', size(snr,1), 'PCHIP'));
cmap                = jet(size(snr,1));

for iCoh = 1:length(snr)
    cidx            = trg_coh == snr(iCoh);
    sc              = scatter(trg_acc(cidx), trg_conf(cidx), 'filled');
    sc.CData        = cmap(iCoh,:);
    sc.SizeData     = 20;
end

% Circ Corr
ax9                         = axes('Position', [clm(1) height(3) dim]); hold on
im                          = imagesc(mcc(sidx,:));

ax9.YLim                    = [.5 size(trg_score,1)+.5];
ax9.XLim                    = [.5 size(trg_score,2)+.5];
ax9.XLabel.String           = 'Coherence [%]';
ax9.YLabel.String           = 'Subject';
% ax7.YTick                   = 1:length(id);
% ax7.YTickLabel              = id(sidx);
ax9.XTick                   = 1:length(snr);
ax9.XTickLabel              = round(snr,2)*100;
ax9.FontSize                = lb_fs;
% ax9.XTickLabelRotation      = 15;

cb                          = colorbar;
cb.Label.String             = 'Response similarity [raw]';
cb.FontSize                 = lg_fs;
colormap(gray(256))


yyaxis right
lp                          = plot(mean(mstr(sidx,:)));
lp.LineWidth                = 2;
lp.Color                    = [1 0 0];
ax9.YAxis(2).Color          = [1 0 0];
ax9.YLim                    = [.5 .8];
set(ax9,'YDir','normal')

print(f, '/Users/fschneider/Desktop/summ', '-r300', '-dpng');