% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/cbrewer/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

close all
clear all

source_pth = '/Users/fschneider/ownCloud/var_plot/';
load([source_pth '/solo_correlation.mat'])
load([source_pth '/solo_performance.mat'])
load([source_pth '/hh_dyad_pairwise_correlation.mat'])
load([source_pth 'hh_dyad_pairwise_performance.mat'])
load([source_pth '/hc_dyad_pairwise_correlation.mat'])
load([source_pth '/hc_dyad_pairwise_performance.mat'])
load([source_pth '/hh_dyad_performance.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Within-dyad effect size %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[acc_df, ecc_df, auc, score, hir, ply1_flag, ply2_flag, raw] = ...
    dyad_effect_size(solo_perf, dyad_pw_perf, dyad_perf);

%%

close all
dest_dir                    = '/Users/fschneider/Desktop/';

for iP = 1:2
    
    f = figure; hold on
    line([0.5 3.5],[.5 .5],'Color', 'k', 'LineStyle', ':', 'LineWidth',1.5)
    
    if iP == 1
        x = [acc_df.solo];
        thresh = std([auc.acc1 auc.acc2])/2; % 1/2 Std as threshold
        more = [auc.acc1(x > thresh) auc.acc2(x < -thresh)]; % AUC(more accurate player)
        less = [auc.acc1(x < -thresh) auc.acc2(x > thresh)]; % AUC(less accurate player)
        same = [auc.acc1(x > -thresh & x < thresh) auc.acc2(x > -thresh & x < thresh)]; % AUC(similar player)
        str = 'Accuracy';
        str2 = 'accurate';
    elseif iP == 2
        x = [ecc_df.solo];
        thresh = std([auc.ecc1 auc.ecc2])/2;
%         x = [acc_df.solo];
%         thresh = std([auc.acc1 auc.acc2])/2;
        more = [auc.ecc1(x > thresh) auc.ecc2(x < -thresh)]; % AUC(more confident player)
        less = [auc.ecc1(x < -thresh) auc.ecc2(x > thresh)]; % AUC(less confident player)
        same = [auc.ecc1(x > -thresh & x < thresh) auc.ecc2(x > -thresh & x < thresh)]; % AUC(similar player)
        str = 'Eccentricity';
        str2 = 'confident';
    end
    
    group = [zeros(1,length(more)) ones(1,length(same)) ones(1,length(less))+1];
    violinplot([more same less],group)
    title([str ' modulation if partner is...'])
    set(gca,'xticklabel',{['less ' str2],['similar'], ['more ' str2]}) % Labels swapped - to get partner perspective
    set(gca,'xticklabelrotation',0)
    set(gca,'fontsize',20)
    ylabel('AUC')
    
    print(f, [dest_dir 'auc_partner_grouping_' str], '-r500', '-dpng');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [acc_df, ecc_df, auc, score, hir, ply1_flag, ply2_flag,raw] = dyad_effect_size(in_solo, in_dyad, dyad_pooled)

for iSubj = 1:size(in_solo,2)
    solo_id{iSubj}          = in_solo{iSubj}.id;
end

for iDyad = 1:size(in_dyad,1)
    % Find solo data of subjects
    idx_ply1                = cellfun(@(x) strcmp(x,in_dyad{iDyad,1}.id),solo_id);
    idx_ply2                = cellfun(@(x) strcmp(x,in_dyad{iDyad,2}.id),solo_id);
    
    score.solo(iDyad,:)   	= [mean(in_solo{idx_ply1}.trg_all_score) mean(in_solo{idx_ply2}.trg_all_score)];
    score.dyad(iDyad,:)  	= [mean(in_dyad{iDyad,1}.trg_all_score) mean(in_dyad{iDyad,2}.trg_all_score)];
    
    hir.solo(iDyad,:)    	= [sum(in_solo{idx_ply1}.trg_all_outc)/length(in_solo{idx_ply1}.trg_all_outc) sum(in_solo{idx_ply2}.trg_all_outc)/length(in_solo{idx_ply2}.trg_all_outc)];
    hir.dyad(iDyad,:)     	= [sum(in_dyad{iDyad,1}.trg_all_outc)/length(in_dyad{iDyad,1}.trg_all_outc) sum(in_dyad{iDyad,2}.trg_all_outc)/length(in_dyad{iDyad,2}.trg_all_outc)];

    % Label subjects of interest
    % soi                     = ID_subj_negAUC(in_solo, dyad_pooled);
    % ply1_flag(iDyad)      	= logical(sum(cellfun(@(x) strcmp(x,in_dyad{iDyad,1}.id),soi)));
    % ply2_flag(iDyad)      	= logical(sum(cellfun(@(x) strcmp(x,in_dyad{iDyad,2}.id),soi)));
    ply1_flag(iDyad) = false;
    ply2_flag(iDyad) = false;
    
    % Raw performance
    raw.acc1(iDyad)         = mean(in_solo{idx_ply1}.macc_trg); % Mean(Median per coherence)
    raw.acc2(iDyad)         = mean(in_solo{idx_ply2}.macc_trg);
    raw.ecc1(iDyad)         = mean(in_solo{idx_ply1}.mecc_state);
    raw.ecc2(iDyad)         = mean(in_solo{idx_ply2}.mecc_state);
        
    raw.dacc1(iDyad)        = mean(in_dyad{iDyad,1}.macc_trg);
    raw.dacc2(iDyad)        = mean(in_dyad{iDyad,2}.macc_trg);
    raw.decc1(iDyad)        = mean(in_dyad{iDyad,1}.mecc_state);
    raw.decc2(iDyad)        = mean(in_dyad{iDyad,2}.mecc_state);
    
    % Performance difference
    acc_df.solo(iDyad)      = mean(in_solo{idx_ply1}.macc_trg) - mean(in_solo{idx_ply2}.macc_trg);
    ecc_df.solo(iDyad)     	= mean(in_solo{idx_ply1}.mecc_state) - mean(in_solo{idx_ply2}.mecc_state);
    acc_df.dyad(iDyad)      = mean(in_dyad{iDyad,1}.macc_trg) - mean(in_dyad{iDyad,2}.macc_trg);
    ecc_df.dyad(iDyad)     	= mean(in_dyad{iDyad,1}.mecc_state) - mean(in_dyad{iDyad,2}.mecc_state);
       
    % Effect size: Solo vs Dyadic
    [auc.ecc1(iDyad), auc.ecc2(iDyad)] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, 'ecc_state');
    [auc.acc1(iDyad), auc.acc2(iDyad)] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, 'acc_trg');

end
end

function [out_ply1, out_ply2] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, str)

for i = 1:length(in_solo{idx_ply1}.(str))
    out1(i)     = getAUROC(in_solo{idx_ply1}.(str){i},in_dyad{iDyad,1}.(str){i});
    out2(i)     = getAUROC(in_solo{idx_ply2}.(str){i},in_dyad{iDyad,2}.(str){i}); 
end

out_ply1 = mean(out1);
out_ply2 = mean(out2);

end

function [out] = getAUROC(in1, in2)

if size(in1,1) == 1
    in1             = in1';
    in2             = in2';
end

lab                 = [zeros(length(in1),1); ones(length(in2),1)];
[~,~,~,out]         = perfcurve(lab,[in1; in2],1);

% figure
% histogram(in1,20)
% hold on
% histogram(in2,20)
% title([median(in1) median(in2) out])

end