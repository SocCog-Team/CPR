close all 
clear all

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/dyad_human_human_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/dyad_human_human_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')

addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

lb_fs = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example subject: eccentricity distributiuon solo vs dyadic %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSubj = 1:size(dyad_perf,2)
    if ~isempty(dyad_perf{iSubj})
        id_dyad{iSubj}	= dyad_perf{iSubj}.id;
    else
        id_dyad{iSubj}	= 'empty';
    end
end

c = 0;
for iExample = [13 20 37]

    
    iSub                	= iExample;
    snr                   	= solo_perf{iSub}.carr;
    coh_col                 = cool(length(snr));
    idx_dy               	= cellfun(@(x) strcmp(x,solo_perf{iSub}.id),id_dyad);
    
    if sum(idx_dy) == 0
        continue
    end
    
    f                      	= figure('units','centimeters','position',[0 0 7 5]);
    c = c+1;
    for iCoh = 1:length(snr)
        lab{iCoh} = num2str(round(snr(iCoh),2)*100);
        auc(c,iCoh) = getAUROC(solo_perf{iSub}.ecc{iCoh}, dyad_perf{idx_dy}.ecc{iCoh});

        for i = 1:2
            if i == 1
                str = 'left';
                dat = solo_perf{iSub}.ecc{iCoh};
            else
                str = 'right';
                dat = dyad_perf{idx_dy}.ecc{iCoh};
            end

            vp = Violin({dat},iCoh,'HalfViolin',str);
            vp.ViolinColor{1} = coh_col(iCoh,:);
            vp.ViolinPlot.LineStyle = 'none';
            vp.ScatterPlot.SizeData = 8;
            vp.ShowBox = false;
            vp.ShowMean = false;
            vp.ShowMedian = false;
            vp.ShowWhiskers = false;
            vp.ShowNotches = false;
            
            if i == 2
                vp.ViolinColor{1} = coh_col(iCoh,:)/2;
            end
        end
    end
    
    ax2                     = gca;
    ax2.XTick = 1:length(snr);
    ax2.XTickLabel = lab;
    ax2.XLabel.String = 'Coherence';
    ax2.YLabel.String = 'Eccentricity';
    ax2.FontSize = lb_fs;
    ax2.YLim = [0 1];
    ax2.Title.String = ['Subject: ' num2str(iExample)];   
    
    print(f, ['/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG2/SFIG2/subject_distribution/SFIG2a' num2str(iExample)], '-r500', '-dpng');
    print(f, ['/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG2/SFIG2/subject_distribution/SFIG2a' num2str(iExample)], '-r500', '-dsvg', '-painters');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Example AUROC %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f                      	= figure('units','centimeters','position',[0 0 7 5]); 
ax = gca; hold on
for iL = 1:size(auc,1)
    pl(iL)               	= plot(auc(iL,:),'Marker', 'x');
    pl(iL).Color          	= [.5 .5 .5];
    pl(iL).LineWidth      	= 1.5;
end

ax.XLim                  	= [1 length(snr)];
ax.YLim                  	= [.1 .9];
ax.FontSize               	= lb_fs;
ax.YLabel.String           	= 'AUROC';
ax.XLabel.String           	= 'Coherence [%]';
ax.XTick                 	= 1:7;
ax.YTick                  	= [.25 .5 .75];
ax.XTickLabel            	= round(snr,2).*100;
ax.XTickLabelRotation      	= 0;
lm                          = line([1 7],[.5 .5], 'Color', 'k', 'LineStyle', ':', 'LineWidth',2);
    
print(f, ['/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG2/SFIG2/SFIG2_auroc'], '-r500', '-dpng');
print(f, ['/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG2/SFIG2/SFIG2_auroc'], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert relevant data to matrix for simplicity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb_fs                       = 8;
lg_fs                       = 8;
lw                          = 3;
frme_ms                     = 1000/120;
alp                         = .4;
dim                         = [.18 .15];
col_dat                     = [0 0 0];
col_ci                      = [.3 0 0];
snr                         = solo_perf{end}.carr;
cnt                         = 0;

for iSubj = 1:length(solo_perf)
    if isempty(solo_perf{iSubj}.hir)
        continue
    end
    
    cnt = cnt+1;
    try
        dyad_hir(cnt,:)       	= dyad_perf{iSubj}.hir;          % Hit rate
        dyad_trg_score(cnt,:)	= dyad_perf{iSubj}.trg_mscore;   % Avg target scores
        dyad_macc(cnt,:)     	= dyad_perf{iSubj}.macc_trg;     % Avg accuracy
        dyad_mecc(cnt,:)    	= dyad_perf{iSubj}.mecc;         % Avg eccentricity
    catch
        dyad_hir(cnt,:)       	= nan;
        dyad_trg_score(cnt,:)	= nan;
        dyad_macc(cnt,:)     	= nan;
        dyad_mecc(cnt,:)    	= nan;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBPLOT: DYADIC - Hit rate raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f                           = figure('units','centimeters','position',[0 0 22.5 5]); hold on
ax1                         = subplot(1,4,1); hold on
[ax1,pl]                    = plotData(ax1,dyad_hir,true,lw,alp,col_dat,col_ci);
ax1.YLim                    = [10 70];
ax1.XLim                    = [1 size(dyad_hir,2)];
ax1.XLabel.String           = 'Coherence [%]';
ax1.YLabel.String           = 'Hit rate [%]';
ax1.XTick                   = 1:length(snr);
ax1.XTickLabel              = round(snr,2)*100;
ax1.FontSize                = lb_fs;
ax1.XTickLabelRotation      = 0;
ax1.XAxis.Visible           = 'on';
% ax1.XGrid                   = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Lag %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 0;
for iSubj = 1:size(dyad_cr,2)
    
    if isempty(dyad_cr{iSubj})
        continue
    end
    
    cnt = cnt+1;
     dyad_mlag_id{cnt} = dyad_cr{iSubj}.id;

    for iCoh = 1:length(snr)
        clear cIdx
        cIdx             	= dyad_cr{iSubj}.coh == snr(iCoh);
        dyad_mlag(cnt,iCoh) = mean(dyad_cr{iSubj}.lag(cIdx));
    end
end

ax2                         = subplot(1,4,2); hold on
[ax2,pl]                    = plotData(ax2,dyad_mlag ./ 1e3,false,lw,alp,col_dat,col_ci);
ax2.XLim                    = [1 size(dyad_mlag,2)];
ax2.YLim                    = [.3 1];
ax2.XLabel.String           = 'Coherence [%]';
ax2.YLabel.String           = 'Lag [s]';
ax2.XTick                   = 1:length(snr);
ax2.XTickLabel              = round(snr,2)*100;
ax2.FontSize                = lb_fs;
ax2.XTickLabelRotation      = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Avg accuracy raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3                         = subplot(1,4,3); hold on
[ax3,pl]                    = plotData(ax3,dyad_macc,true,lw,alp,col_dat,col_ci);
ax3.XLim                    = [1 size(dyad_macc,2)];
ax3.YLim                    = [30 100];
ax3.XLabel.String           = 'Coherence [%]';
ax3.YLabel.String           = 'Accuracy [%]';
ax3.XTick                   = 1:length(snr);
ax3.XTickLabel              = round(snr,2)*100;
ax3.FontSize                = lb_fs;
ax3.XTickLabelRotation      = 0;
ax3.XAxis.Visible           = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Avg eccentricity raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax4                         = subplot(1,4,4); hold on
[ax4,pl]                    = plotData(ax4,dyad_mecc,true,lw,alp,col_dat,col_ci);
ax4.XLim                    = [1 size(dyad_mecc,2)];
ax4.YLim                    = [20 100];
ax4.XLabel.String           = 'Coherence [%]';
ax4.YLabel.String           = 'Eccentricity [%]';
ax4.XTick                   = 1:length(snr);
ax4.XTickLabel              = round(snr,2)*100;
ax4.FontSize                = lb_fs;
ax4.XTickLabelRotation      = 0;
ax4.XAxis.Visible           = 'on';

print(f, ['/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG2/SFIG2/SFIG2_dyad'], '-r500', '-dpng');
print(f, ['/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG2/SFIG2/SFIG2_dyad'], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = getAUROC(in1, in2)

if size(in1,1) == 1
    in1         = in1';
    in2         = in2';
end

lab          	= [zeros(length(in1),1); ones(length(in2),1)];
[~,~,~,out]     = perfcurve(lab,[in1; in2],1);

end

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
