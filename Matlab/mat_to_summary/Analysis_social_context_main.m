% Mains cript
close all
clear all

addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/mat_to_summary/


dest_dir                    = '/Users/fschneider/Desktop/';
file_pth                    = '/Users/fschneider/Desktop/social_context_study/';
num_dyads                   = 7;

for iDyad = 1:num_dyads
    % Extract all files
    h5Files = dir(fullfile([file_pth 'Dyad' num2str(iDyad)], '*.h5'));
    
    for iFile = 1:length(h5Files)
        % Load file
        d                        	= MW_readH5([file_pth '/Dyad' num2str(iDyad) '/' h5Files(iFile).name]); % ...load .h5 file
        exp_info                    = strsplit(h5Files(iFile).name,'_');
        
        % Index variables of interest
        idx.cOn                     = d.event == 'TRIAL_start';
        idx.cEnd                    = d.event == 'TRIAL_end';
        idx.frame                   = d.event == 'STIM_displayUpdate';
        idx.RDP_onset               = d.event == 'STIM_RDP_onset';
        idx.RDP_dir                 = d.event == 'STIM_RDP_direction';
        idx.RDP_coh                 = d.event == 'STIM_RDP_coherence';
        idx.RDP_dot                	= d.event == 'STIM_RDP_dotPositions';
        idx.trg_on                  = d.event == 'STIM_target_onset';
        idx.JS_dir                  = d.event == 'IO_joystickDirection';
        idx.JS_str                  = d.event == 'IO_joystickStrength';
        idx.JS2_dir                 = d.event == 'IO_joystickDirection2';
        idx.JS2_str                 = d.event == 'IO_joystickStrength2';
        idx.fixation              	= d.event == 'IO_fixation_flag';
        idx.outcome                 = d.event == 'TRIAL_outcome';
        idx.outcome2                = d.event == 'TRIAL_outcome2';
        idx.trg                     = d.event == 'TRIAL_reactionTrigger';
        idx.cum_score            	= d.event == 'INFO_Score';
        idx.cum_score2              = d.event == 'INFO_Score2';
        idx.trg_score              	= d.event == 'INFO_TargetScore';
        idx.trg_score2           	= d.event == 'INFO_TargetScore2';
        idx.bonus1                  = d.event == 'INFO_bonus_ply1_cents';
        idx.bonus2                  = d.event == 'INFO_bonus_ply2_cents';
        idx.performance             = d.event == 'INFO_performance_percent';
        idx.performance2            = d.event == 'INFO_performance_percent2';
        
        % Extract subject behavior
        [out_p1{iFile}, out_p2{iFile}] = CPR_extract_response_profile(d,idx,exp_info);
    end
    
    dyad_summary{iDyad} = [out_p1 out_p2];
end

%% Show average data for each dyadic player and condition 

close all
f = figure('units','centimeters','position',[0 0 40 40]);
cmap = jet(num_dyads)./1.3;
cnt = 0;
lw = 1;
sz = 75;

for iDyad = 1:num_dyads
    
    ax1 = subplot(num_dyads,4,(cnt*4)+1); hold on;
    ax2 = subplot(num_dyads,4,(cnt*4)+2); hold on;
    ax3 = subplot(num_dyads,4,(cnt*4)+3); hold on;
    ax4 = subplot(num_dyads,4,(cnt*4)+4); hold on;
    cnt = cnt+1;
    dyad{iDyad} = dyad_summary{iDyad}{1}.dyad;
    
    for iPlayer = 1:2
        for iCondition = 1:3  
            % Data from each block of a given exp. condition (i.e. neutral_psy3 and neutral_psy4)
            out = get_summary(dyad_summary{iDyad}, iPlayer, iCondition);
            bonus_sum(iPlayer,iCondition,iDyad) = mean([out{1}.score out{2}.score]);
            hir(iPlayer,iCondition,iDyad) = mean([out{1}.hir out{2}.hir]);
            acc_mean(iPlayer,iCondition,iDyad) = mean([out{1}.acc_mean out{2}.acc_mean]);
            ecc_mean(iPlayer,iCondition,iDyad) = mean([out{1}.ecc_mean out{2}.ecc_mean]);
        end
        
        axes(ax1)
        plot([1 2 3], bonus_sum(iPlayer,:,iDyad),'Color',cmap(iDyad,:),'LineWidth', lw)
        scatter([1 2 3], bonus_sum(iPlayer,:,iDyad),'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
        
        axes(ax2)
        plot([1 2 3],hir(iPlayer,:,iDyad),'Color',cmap(iDyad,:),'LineWidth', lw)
        scatter([1 2 3], hir(iPlayer,:,iDyad),'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
        
        axes(ax3)
        plot([1 2 3],acc_mean(iPlayer,:,iDyad),'Color',cmap(iDyad,:),'LineWidth', lw)
        scatter([1 2 3], acc_mean(iPlayer,:,iDyad),'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
        
        axes(ax4)
        plot([1 2 3],ecc_mean(iPlayer,:,iDyad),'Color',cmap(iDyad,:),'LineWidth', lw)
        scatter([1 2 3], ecc_mean(iPlayer,:,iDyad),'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)

    end
    change_axes(ax1,'reward score', iDyad, dyad{iDyad}, true)
    change_axes(ax2,'hit rate', iDyad, dyad{iDyad})
    change_axes(ax3,'accuracy',iDyad, dyad{iDyad})
    change_axes(ax4,'eccentricity',iDyad, dyad{iDyad})
end

print(f, [dest_dir '/raw_perf'], '-r500', '-dpng');

%% Show normalised data for each dyadic player and condition

close all
f = figure('units','centimeters','position',[0 0 40 15]);

for iPlot = 1:4
ax(iPlot) = subplot(1,4,iPlot); hold on; line([1 2],[0 0], 'LineStyle',':', 'Color','k', 'LineWidth',2);
end

cmap = jet(num_dyads)./1.3;
lw = 1;
sz = 75;
xvec = 1:2;

for iDyad = 1:num_dyads    
    for iPlayer = 1:2
        axes(ax(1))
        plot(xvec, bonus_sum(iPlayer,2:3,iDyad)-bonus_sum(iPlayer,1,iDyad),'Color',cmap(iDyad,:),'LineWidth', lw)
        scatter(xvec, bonus_sum(iPlayer,2:3,iDyad)-bonus_sum(iPlayer,1,iDyad),'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
        
        axes(ax(2))
        plot(xvec, hir(iPlayer,2:3,iDyad)-hir(iPlayer,1,iDyad),'Color',cmap(iDyad,:),'LineWidth', lw)
        scatter(xvec, hir(iPlayer,2:3,iDyad)-hir(iPlayer,1,iDyad),'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
        
        axes(ax(3))
        plot(xvec, acc_mean(iPlayer,2:3,iDyad)-acc_mean(iPlayer,1,iDyad),'Color',cmap(iDyad,:),'LineWidth', lw)
        scatter(xvec, acc_mean(iPlayer,2:3,iDyad)-acc_mean(iPlayer,1,iDyad),'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
        
        axes(ax(4))
        plot(xvec, ecc_mean(iPlayer,2:3,iDyad)-ecc_mean(iPlayer,1,iDyad),'Color',cmap(iDyad,:),'LineWidth', lw)
        scatter(xvec, ecc_mean(iPlayer,2:3,iDyad)-ecc_mean(iPlayer,1,iDyad),'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
    end
end

for iPlot = 1:4
    if iPlot == 1
        ystr = 'reward score';
        ylim = [-22 22];
    elseif iPlot == 2
        ystr = 'hit rate';
        ylim = [-.14 .14];
    elseif iPlot == 3
        ystr = 'accuracy';
        ylim = [-.04 .04];
    elseif iPlot == 4
        ystr = 'eccentricity';
        ylim = [-.1 .1];
    end
        
    change_axes_dff_plot(ax(iPlot),ystr,ylim,[1 2],{'Coop','Comp'})
end

print(f, [dest_dir '/normalised_perf'], '-r500', '-dpng');

%% Show difference between players

close all
f = figure('units','centimeters','position',[0 0 40 15]);

for iPlot = 1:4
    ax(iPlot) = subplot(1,4,iPlot); hold on;
end

cmap = jet(num_dyads)./1.3;
cnt = 0;
lw = 1;
sz = 75;
xvec = 1:3;

for iDyad = 1:num_dyads
    clear ydat
    
    axes(ax(1))
    ydat = abs(bonus_sum(1,:,iDyad)-bonus_sum(2,:,iDyad));
    p = plot(xvec, ydat,'Color',cmap(iDyad,:),'LineWidth', lw);
    scatter(xvec, ydat,'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
    
    axes(ax(2))
    ydat = abs(hir(1,:,iDyad)-hir(2,:,iDyad));
    plot(xvec, ydat,'Color',cmap(iDyad,:),'LineWidth', lw)
    scatter(xvec, ydat,'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
    
    axes(ax(3))
    ydat = abs(acc_mean(1,:,iDyad)-acc_mean(2,:,iDyad));
    plot(xvec, ydat,'Color',cmap(iDyad,:),'LineWidth', lw)
    scatter(xvec, ydat,'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
    
    axes(ax(4))
    ydat = abs(ecc_mean(1,:,iDyad)-ecc_mean(2,:,iDyad));
    plot(xvec, ydat,'Color',cmap(iDyad,:),'LineWidth', lw)
    scatter(xvec, ydat,'MarkerEdgeColor',cmap(iDyad,:),'LineWidth', lw*3, 'Marker', 'x', 'SizeData', sz)
end

for iPlot = 1:4
    if iPlot == 1
        ystr = 'reward score';
        ylim = [0 22];
    elseif iPlot == 2
        ystr = 'hit rate';
        ylim = [0 .14];
    elseif iPlot == 3
        ystr = 'accuracy';
        ylim = [0 .05];
    elseif iPlot == 4
        ystr = 'eccentricity';
        ylim = [0 .15];
    end

    change_axes_dff_plot(ax(iPlot),ystr,ylim,[1 2 3],{'Neutral','Coop','Comp'})
end
    
print(f, [dest_dir '/abs_partner_dff'], '-r500', '-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = get_summary(dyad_summary, iPlayer, iCondition)
dyad_id = dyad_summary{1}.dyad;

if iPlayer == 1
    player_idx = cellfun(@(x) strcmp(x.player_id,dyad_id(1:3)),dyad_summary);
elseif iPlayer == 2
    player_idx = cellfun(@(x) strcmp(x.player_id,dyad_id(4:end)),dyad_summary);
end

if iCondition == 1
    condition_idx = cellfun(@(x) strcmp(x.condition,'CPRneutral'),dyad_summary);
elseif iCondition == 2
    condition_idx = cellfun(@(x) strcmp(x.condition,'CPRcooperation'),dyad_summary);
elseif iCondition == 3
    condition_idx = cellfun(@(x) strcmp(x.condition,'CPRcompetition'),dyad_summary);
end

out = dyad_summary(condition_idx & player_idx);

end

function change_axes(ax,ystr,iDyad,dyad_id,yflag)
if nargin < 5
    yflag = false;
end

if nargin < 4
    dyad_id = [];
end

ax.FontSize = 14;
ax.XTickLabel = {'Neutr','Coop','Comp'};
ax.XTickLabelRotation = 30;

if yflag == true
    ax.YLabel.String = {['Dyad ' num2str(iDyad)]; dyad_id};
end
if iDyad == 1
    ax.Title.String = ystr;
end

if iDyad < 7
    ax.XAxis.Visible = 'off';
end
end

function change_axes_dff_plot(ax,ystr,ylim,xtick, xticklab)
ax.YLim = ylim;
ax.XTick = xtick;
ax.XTickLabel = xticklab;
ax.Title.String = ystr;
ax.FontSize = 14;
ax.XTickLabelRotation = 30;
end



