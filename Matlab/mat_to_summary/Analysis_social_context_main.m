% Mains cript
close all
clear all
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/mat_to_summary/

file_pth                    = '/Users/fschneider/Desktop/social_context_study/';
num_dyads                   = 3;

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
    
    dyad_summary = [out_p1 out_p2];
    
    for iPlayer = 1:2
        for iCondition = 1:3
            out = get_summary(dyad_summary, exp_info{2}, iPlayer, iCondition);
                        
%             ax = subplot(1,3,iCondition); hold on
           
            if iPlayer == 1
                col = 'k';
            else
                col = 'r';
            end
            
            bonus_sum(iPlayer,iCondition,iDyad) = out{1}.bonus+out{2}.bonus;
            
%             plot(out{1}.coh_sorted.acc_mean, col)
%             plot(out{2}.coh_sorted.acc_mean, col)
%             ax.Title.String = out{1}.condition;
%             ax.YLim = [0 1];
%             grid on
        end
    end
end

% Compare social conditions
% 3x3 matrix: acc, ecc, lag for neutral, cooperation and competition

function out = get_summary(dyad_summary, dyad_id, iPlayer, iCondition)

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