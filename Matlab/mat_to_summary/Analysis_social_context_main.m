% Mains cript
close all
clear all
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/mat_to_summary/

file_pth                    = '/Users/fschneider/Desktop/social_context_study/';
num_dyads                   = 5;

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

%%
close all
f = figure('units','centimeters','position',[0 0 40 40]);
cmap = jet(5)./1.3;
cnt = 0;
for iDyad = 1:num_dyads
    ax1 = subplot(5,4,(cnt*4)+1); hold on;
    ax2 = subplot(5,4,(cnt*4)+2); hold on;
    ax3 = subplot(5,4,(cnt*4)+3); hold on;
    ax4 = subplot(5,4,(cnt*4)+4); hold on;
    cnt = cnt+1;
    dyad{iDyad} = dyad_summary{iDyad}{1}.dyad;
    for iPlayer = 1:2
        
        for iCondition = 1:3
            out = get_summary(dyad_summary{iDyad}, iPlayer, iCondition);
            
            bonus_sum(iPlayer,iCondition,iDyad) = out{iPlayer}.score;
            hir(iPlayer,iCondition,iDyad) = out{iPlayer}.hir;
            acc_mean(iPlayer,iCondition,iDyad) = mean(out{iPlayer}.acc_mean);
            ecc_mean(iPlayer,iCondition,iDyad) = mean(out{iPlayer}.ecc_mean);
        end
        
        axes(ax1)
        plot([1 2 3], bonus_sum(iPlayer,:,iDyad),'Color',cmap(iDyad,:))
        
        axes(ax2)
        plot([1 2 3],hir(iPlayer,:,iDyad),'Color',cmap(iDyad,:))
        
        axes(ax3)
        plot([1 2 3],acc_mean(iPlayer,:,iDyad),'Color',cmap(iDyad,:))
        
        axes(ax4)
        plot([1 2 3],ecc_mean(iPlayer,:,iDyad),'Color',cmap(iDyad,:))
    end
    change_axes(ax1,'reward score', iDyad, true)
    change_axes(ax2,'hit rate', iDyad)
    change_axes(ax3,'accuracy',iDyad)
    change_axes(ax4,'eccentricity',iDyad)
end

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

function change_axes(ax,ystr,iDyad,yflag)
if nargin < 4
    yflag = false;
end

ax.FontSize = 14;
ax.XTickLabel = {'Neutr','Coop','Comp'};
ax.XTickLabelRotation = 30;

if yflag == true
    ax.YLabel.String = ['Dyad ' num2str(iDyad)];
end
if iDyad == 1
    ax.Title.String = ystr;
end

if iDyad < 5
    ax.XAxis.Visible = 'off';
end
end