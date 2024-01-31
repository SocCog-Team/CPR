% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_pairwise_performance.mat')

% Import subject summary spreadsheet
pth                         = '/Volumes/T7_Shield/CPR_psychophysics/';      % Local hard drive
x                           = readtable([pth 'Subjects_summary.xlsx']);     % Spreadsheet
sbj_lst                     = x.Abbreviation;                               % Subject ID list
sbj_lst(cellfun(@isempty,sbj_lst)) = [];

nSample                     = 30;                                           % Time window size [samples]
nLag                        = 150;                                          % Cross-correlation lag                               

%% Convert relevant data to matrix for simplicity

cnt                         = 0;
for iSubj = 1:length(sbj_lst)
    if isempty(solo_perf{iSubj}.hir)
        continue
    end
    
    cnt                   	= cnt + 1;
    solo_hir(cnt,:)        	= solo_perf{iSubj}.hir;             % Hit rate
    solo_trg_score(cnt,:)  	= solo_perf{iSubj}.trg_mscore;      % Avg target scores
    solo_macc(cnt,:)       	= solo_perf{iSubj}.macc_trg;        % Avg accuracy
    solo_mecc(cnt,:)       	= solo_perf{iSubj}.mecc_state;  	% Avg eccentricity
end

%% PLOT

f                           = figure('units','centimeters','position',[0 0 20 4]);
height                      = [linspace(.555, .056,4) .1];
colmn                       = linspace(.075, .85,4);
lb_fs                       = 8;
lg_fs                       = 8;
lw                          = 3;
frme_ms                     = 1000/120;
alp                         = .4;
dim                         = [.18 .15];
col_dat                     = [0 0 0];
col_ci                      = [.3 0 0];
snr                         = solo_perf{end}.carr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Hit rate raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1                         = subplot(1,4,1); hold on
[ax1,pl]                    = plotData(ax1,solo_hir,true,lw,alp,col_dat,col_ci);
ax1.YLim                    = [10 70];
ax1.XLim                    = [1 size(solo_hir,2)];
ax1.XLabel.String           = 'Coherence [%]';
ax1.YLabel.String           = 'Hit rate [%]';
ax1.YLabel.Position(1)      = -.4;
ax1.XTick                   = 1:length(snr);
ax1.XTickLabel              = round(snr,2)*100;
ax1.FontSize                = lb_fs;
ax1.XTickLabelRotation      = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Lag %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSubj = 1:size(solo_cr,2)
    solo_mlag_id{iSubj} = solo_cr{iSubj}.id;

    for iCoh = 1:length(snr)
        clear cIdx
        cIdx             	= solo_cr{iSubj}.coh == snr(iCoh);
        solo_mlag(iSubj,iCoh) = mean(solo_cr{iSubj}.lag(cIdx));
    end
end

ax2                         = subplot(1,4,4); hold on
[ax2,pl]                    = plotData(ax2,solo_mlag ./ 1e3,false,lw,alp,col_dat,col_ci);
ax2.XLim                    = [1 size(solo_mlag,2)];
ax2.YLim                    = [.3 1];
ax2.XLabel.String           = 'Coherence [%]';
ax2.YLabel.String           = 'Lag [s]';
ax2.YLabel.Position(1)      = -.4;
ax2.XTick                   = 1:length(snr);
ax2.XTickLabel              = round(snr,2)*100;
ax2.FontSize                = lb_fs;
ax2.XTickLabelRotation      = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Avg accuracy raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3                         = subplot(1,4,2); hold on
[ax3,pl]                    = plotData(ax3,solo_macc,true,lw,alp,col_dat,col_ci);
ax3.XLim                    = [1 size(solo_macc,2)];
ax3.YLim                    = [30 100];
ax3.XLabel.String           = 'Coherence [%]';
ax3.YLabel.String           = 'Accuracy [%]';
ax3.YLabel.Position(1)      = -.4;
ax3.XTick                   = 1:length(snr);
ax3.XTickLabel              = round(snr,2)*100;
ax3.FontSize                = lb_fs;
ax3.XTickLabelRotation      = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Avg eccentricity raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax4                         = subplot(1,4,3); hold on
[ax4,pl]                    = plotData(ax4,solo_mecc,true,lw,alp,col_dat,col_ci);
ax4.XLim                    = [1 size(solo_mecc,2)];
ax4.YLim                    = [20 100];
ax4.XLabel.String           = 'Coherence [%]';
ax4.YLabel.String           = 'Eccentricity [%]';
ax4.YLabel.Position(1)      = -.4;
ax4.XTick                   = 1:length(snr);
ax4.XTickLabel              = round(snr,2)*100;
ax4.FontSize                = lb_fs;
ax4.XTickLabelRotation      = 0;

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/raw/';
print(f, [dest_dir '/FIG1_solo'], '-r500', '-dpng');
print(f, [dest_dir '/FIG1_solo'], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ax,pl] = plotData(ax,dat,scale_flag,lw,alp,col_dat,col_ci)

% Remove NaNs and zero rows
dat(sum(isnan(dat),2)>0,:)  = [];
dat(sum(dat==0,2)>0,:)      = [];

if scale_flag
    dat                     = dat.*100;
end

% Boostrap confidence intervals
nRep                        = 1000;
[CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.05);

% Prepare filled area
vec                         = 1:length(CI);
x_spacing                   = [vec fliplr(vec)];
ci                          = [CI(1,:) fliplr(CI(2,:))];

% Plot subject-wise data
hold on
for iL = 1:size(dat,1)
    pl(iL)                  = plot(dat(iL,:), 'LineWidth', lw/lw, 'Color', [.5 .5 .5 alp]);
end

% Overlay confidence intervals
fl                          = fill(x_spacing,ci,col_ci,'EdgeColor','none', 'FaceAlpha', alp);

% Plot mean curve
pl                          = plot(mean(dat), 'LineWidth', lw/1.5, 'Color', col_dat);
end
