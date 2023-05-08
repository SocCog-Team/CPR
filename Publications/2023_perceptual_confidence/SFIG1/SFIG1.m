%% FIG 1 %%
close all
clear all

% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/WIP/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

%% Import and save RT data

pth                     = '/Volumes/T7_Shield/CPR_psychophysics/';
fname                   = 'Subjects_summary.xlsx';
x                       = readtable([pth fname]);
sbj_lst                 = x.Abbreviation;
dte_format              = 'yyyymmdd';
sbj_lst(cellfun(@isempty,sbj_lst)) = [];

f                    	= figure('units','normalized','position',[0 0 .5 1]);
clm                     = linspace(.1,.79,3);
height                  = linspace(.77,.06,4);
dim                     = [.2 .2];
lb_fs                   = 14;
n_step                  = 30;
vec                     = 0:n_step:360;
alph                    = .4;
col                     = [.5 .5 .5];
lw                      = 3;

ax0                  	= polaraxes('Position', [clm(1) height(1) dim]); hold on
for iSubj = 1:length(sbj_lst)
    load([pth sbj_lst{iSubj} '/summary/' sbj_lst{iSubj} '_RT.mat'])
    
    for iDir = 1:size(vec,2)-1
        didx            = rt.trg_dir > vec(iDir) & rt.trg_dir <= vec(iDir+1);
        mrt(iDir)       = median(rt.dat(didx));
        tdir(iDir)      = mean([vec(iDir) vec(iDir+1)]);
    end
    
    pl                  = polarplot(deg2rad([tdir tdir(1)]), [mrt mrt(1)]);
    pl.LineWidth        = lw;
    pl.Color            = [col alph];
    
    hold on
end

ax0                     = gca;
ax0.RLim                = [0 500];
ax0.FontSize            = lb_fs;

%% SNR %%%

ax1                  	= axes('Position', [clm(2) height(1) dim(1)*2 dim(2)]); hold on
lw                      = 3;

for iSubj = 1:length(sbj_lst)
    try
        load([pth sbj_lst{iSubj} '/summary/' sbj_lst{iSubj} '_SNRfit.mat'])
    catch
        continue
    end
    
    psy_func.model_snr(1) = .0001;
    psy_func.snr(1)       = .01;

    pl                      = plot((psy_func.model_snr),psy_func.model,'-','color',[col alph],'linewidth',lw);
    sc                      = scatter(psy_func.snr,psy_func.hir);
    sc.Marker               = '.';
    sc.MarkerEdgeColor      = col;
    sc.MarkerEdgeAlpha      = alph;
    sc.MarkerFaceAlpha      = alph;
    sc.SizeData             = 100;
end

ax1.XScale                  = 'log';
ax1.YLim                    = [0 1];
ax1.XLim                    = [.01 1];
ax1.FontSize                = lb_fs;
ax1.XTick                   = psy_func.snr;
ax1.YTick                   = [0 .25 .5 .75 1];
ax1.XLabel.String           = 'Coherence [%]';
ax1.YLabel.String           = '% correct';
ln                          = line(ax1.XLim,[.25 .25],'Color', 'k', 'LineStyle', '--');

for i = 1:length(ax1.XTickLabel)
    if i == 1
        ax1.XTickLabel{i}   = 0;
    else
        ax1.XTickLabel{i} 	= round(psy_func.snr(i),2);
    end
end

grid on
grid minor 
grid minor

%% Score over time - different subjects + session types labelled

ax1                  	= axes('Position', [clm(1) height(2) dim]); hold on
cnt                     = 0; 

for iSubj = 1:length(sbj_lst)
    
    cd([pth sbj_lst{iSubj} '/summary/'])
    mat_files                       = dir('*.mat');
    
    for iFile = 1:length(mat_files)
        fstr                        = strsplit(mat_files(iFile).name,'_');
        
        if size(fstr,2) < 3 || length(fstr{2}) ~= 3
            continue
        end
    
        if contains(fstr{3},'CPR')
            load(mat_files(iFile).name)
            [cum_score, trg_score]  = get_trg_score(t);
            
            cnt                     = cnt + 1;
            id{cnt}                 = lower(fstr{2});
            last_score(cnt)         = cum_score(end);
            dte(cnt)                = str2num(fstr{1});
            
            if strcmp(fstr{3},'CPRsolo')
                flg(cnt)          	= 1;
            elseif strcmp(fstr{3},'CPRagent')
                flg(cnt)         	= 2;
            else
                flg(cnt)         	= 999;
            end
        end
    end
end

id(last_score < 50) = [];
dte(last_score < 50) = [];
last_score(last_score < 50) = [];

cl = jet(length(sbj_lst));
for iSubj = 1:length(sbj_lst)
    sidx                         	= cellfun(@(x) strcmp(x, lower(sbj_lst{iSubj})), id);
    exp_dte                       	= dte(sidx);
    sc                          	= scatter(1:length(last_score(sidx)), last_score(sidx));
    sc.MarkerFaceColor              = cl(iSubj,:);
    sc.MarkerEdgeColor            	= 'none';
    
    P                               = polyfit(1:length(last_score(sidx)),last_score(sidx),1);
    slpe(iSubj)                     = P(1);
    yfit                            = P(1)*(1:length(last_score(sidx)))+P(2);  % P(1) is the slope and P(2) is the intercept
    pf                              = plot(1:length(last_score(sidx)),yfit,'k');
end

ax2                 	= axes('Position', [clm(2) height(2) dim]); hold on
bx                      = boxplot(slpe);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% print(f, '/Users/fschneider/ownCloud/Documents/Publications/CPR_psychophysics/Figures/FIG2/FIG2', '-r400', '-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tmp_score, tscore] = get_trg_score(t)
c                           = 0;
for iState = 1:size(t.trg_ts,1)
    if length(t.trg_ts{iState}) ~= length(t.trg_hit{iState})
        continue
    end
    
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
end