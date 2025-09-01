close all
clear all

% Adjust path
source_pth = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/var_plot/';
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
alignment_str = 'state';
% alignment_str = 'trg';

[acc_df, ecc_df, auc, score, hir, raw] = dyad_effect_size(solo_perf, dyad_pw_perf, alignment_str);

%%% PLOT PARAMS %%%
dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/FIG_pairwise_modulation/raw/';
rcol                       	= [.6 .1 .1];
scol                        = [.1 .1 .1];
hcol                        = [.1 .1 .1];
lw                         	= 1;
lb_fs                    	= 8;
sc_size                     = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Dyad-wise modulation %%% SUPPLEMENTARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                           = figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
ofs                         = .025;
nbin                        = 12;
dyad_pw_modulation          = computer_player_modulation(solo_perf,hc_dyad_pw_perf,dyad_pw_perf);

for iPlot = 1:4
    ax                      = subplot(2,2,iPlot); hold on
    
    ln                     	= line([0 1],[.5 .5]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
    
    ln                     	= line([.5 .5],[-1 1]);
    ln.Color              	= [0 0 0];
    ln.LineWidth           	= lw;
    ln.LineStyle           	= ':';
        
        
    if iPlot == 1
        %%% Subject-wise social modulation of accuracy and eccentricity %%%
        sc                      = scatter([auc.acc1 auc.acc2],[auc.ecc1 auc.ecc2]);
        sc.MarkerFaceColor      = [hcol];
        sc.MarkerEdgeColor      = 'none';
        sc.MarkerFaceAlpha      = .5;
        ax.XLabel.String        = 'Accuracy change [AUC]';
        ax.YLabel.String        = 'Eccentricity change [AUC]';
        ax.XLim                 = [.3 .6];
        ax.XTick                = [.3:.1:.6];
        ax.YLim                 = [0 1];
        
        [r,pv]               	= corrcoef([auc.acc1 auc.acc2],[auc.ecc1 auc.ecc2]);
        r                       = r(2)
        pv                      = pv(2)
        
        linearCoefficients      = polyfit([auc.acc1 auc.acc2],[auc.ecc1 auc.ecc2], 1);
        x_fit                   = linspace(-1, 1, 50);
        y_fit                   = polyval(linearCoefficients, x_fit);
        pl                      = plot(x_fit, y_fit, '-', 'LineWidth', 1);

        pl.Color                = rcol;
        tx                      = text(.5,.15,{['r=' num2str(round(r,2))];['p=' num2str(round(pv,2))]});
        tx.Color                = hcol;
        tx.FontSize             = 8;
        
    elseif iPlot == 3
        %%% Average dyad-wise accuracy modulation as function of solo eccentricity difference %%%       
        sc                      = scatter(mean([auc.acc1; auc.acc2]),abs(raw.ecc1 - raw.ecc2));
        sc.MarkerFaceColor      = [hcol];
        sc.MarkerEdgeColor      = 'none';
        sc.MarkerFaceAlpha      = .5;
        ax.YLabel.String        = 'Abs ecc dff [Solo]';
        ax.XLabel.String        = 'Avg. accuracy change [AUC]';
        ax.XLim                 = [.4 .75];
        ax.XTick                = [.4:.1:.7];
        ax.YLim                 = [0 .4];
        
        [r,pv]               	= corrcoef(mean([auc.acc1; auc.acc2]),abs(raw.ecc1 - raw.ecc2));
        r                       = r(2)
        pv                      = pv(2)
        
        linearCoefficients      = polyfit(mean([auc.acc1; auc.acc2]),abs(raw.ecc1 - raw.ecc2), 1);
        x_fit                   = linspace(-1, 1, 50);
        y_fit                   = polyval(linearCoefficients, x_fit);
        pl                      = plot(x_fit, y_fit, '-', 'LineWidth', 1);

        pl.Color                = rcol;
        tx                      = text(.5,.15,{['r=' num2str(round(r,2))];['p=' num2str(round(pv,2))]});
        tx.Color                = hcol;
        tx.FontSize             = 8;
            
    elseif iPlot == 2
        %%% Solo vs HC comparison
        ydat = mean(dyad_pw_modulation.auc_ecc_SC,2);
        xdat = mean(dyad_pw_modulation.auc_acc_SC,2);
        
        sc                      = scatter(xdat, ydat);
        sc.MarkerFaceColor      = [hcol];
        sc.MarkerEdgeColor      = 'none';
        sc.MarkerFaceAlpha      = .5;
        
        [r,pv] = corrcoef(xdat, ydat);
        r                       = r(2)
        pv                      = pv(2)
        
        linearCoefficients      = polyfit(xdat,ydat, 1);
        x_fit                   = linspace(-1, 1, 50);
        y_fit                   = polyval(linearCoefficients, x_fit);
        pl                      = plot(x_fit, y_fit, '-', 'LineWidth', 1);
        pl.Color                = rcol;
        tx                      = text(.6,.15,{['r=' num2str(round(r,2))];['p=' num2str(round(pv,2))]});
        tx.Color                = hcol;
        tx.FontSize             = 8;
        
        ax.YLabel.String        = 'Eccentricity change [AUC]';
        ax.XLabel.String        = 'Accuracy change [AUC]';
        ax.YLim                 = [0 1];
        ax.XLim                 = [.4 .75];
        ax.XTick                = [.4:.1:.7];
        
    elseif iPlot == 4
        %%% HH vs HC modulation comparison
        ydat = mean(dyad_pw_modulation.auc_ecc_HC,2);
        xdat = mean(dyad_pw_modulation.auc_acc_HC,2);
        
        sc                      = scatter(xdat, ydat);
        sc.MarkerFaceColor      = [hcol];
        sc.MarkerEdgeColor      = 'none';
        sc.MarkerFaceAlpha      = .5;
        
        [r,pv] = corrcoef(xdat, ydat);
        r                       = r(2)
        pv                      = pv(2)
        
        linearCoefficients      = polyfit(xdat,ydat, 1);
        x_fit                   = linspace(-1, 1, 50);
        y_fit                   = polyval(linearCoefficients, x_fit);
        pl                      = plot(x_fit, y_fit, '-', 'LineWidth', 1);
        pl.Color                = rcol;
        tx                      = text(.6,.15,{['r=' num2str(round(r,2))];['p=' num2str(round(pv,2))]});
        tx.Color                = hcol;
        tx.FontSize             = 8;
        
        ax.YLabel.String      	= 'AUC: Eccentricity';
        ax.XLabel.String       	= 'AUC: Accuracy';
        ax.YLim                 = [0 1];
        ax.XLim                 = [.4 .75];  
        ax.XTick                = [.4:.1:.7];
    end
    
    sc.SizeData             = sc_size;
    ax.Position(1)              = ax.Position(1) - ofs;
    ax.FontSize                 = lb_fs;
end

dest_dir = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/FIG_modulation_correlation/raw/';
print(f, [dest_dir '/FIG_auc_dyadwise_' alignment_str], '-r500', '-dpng');
print(f, [dest_dir '/FIG_auc_dyadwise_' alignment_str], '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [acc_df, ecc_df, auc, score, hir, raw] = dyad_effect_size(in_solo, in_dyad, str)

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
    
    % Raw performance
    raw.acc1(iDyad)         = mean(in_solo{idx_ply1}.(['macc_' str])); % Mean(Median per coherence)
    raw.acc2(iDyad)         = mean(in_solo{idx_ply2}.(['macc_' str]));
    raw.ecc1(iDyad)         = mean(in_solo{idx_ply1}.(['mecc_' str]));
    raw.ecc2(iDyad)         = mean(in_solo{idx_ply2}.(['mecc_' str]));
        
    raw.dacc1(iDyad)        = mean(in_dyad{iDyad,1}.(['macc_' str]));
    raw.dacc2(iDyad)        = mean(in_dyad{iDyad,2}.(['macc_' str]));
    raw.decc1(iDyad)        = mean(in_dyad{iDyad,1}.(['mecc_' str]));
    raw.decc2(iDyad)        = mean(in_dyad{iDyad,2}.(['mecc_' str]));
    
    % Performance difference
    acc_df.solo(iDyad)      = mean(in_solo{idx_ply1}.(['macc_' str])) - mean(in_solo{idx_ply2}.(['macc_' str]));
    ecc_df.solo(iDyad)     	= mean(in_solo{idx_ply1}.(['mecc_' str])) - mean(in_solo{idx_ply2}.(['mecc_' str]));
    acc_df.dyad(iDyad)      = mean(in_dyad{iDyad,1}.(['macc_' str])) - mean(in_dyad{iDyad,2}.(['macc_' str]));
    ecc_df.dyad(iDyad)     	= mean(in_dyad{iDyad,1}.(['mecc_' str])) - mean(in_dyad{iDyad,2}.(['mecc_' str]));
       
    % Effect size: Solo vs Dyadic
    [auc.ecc1(iDyad), auc.ecc2(iDyad)] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, ['ecc_' str]);
    [auc.acc1(iDyad), auc.acc2(iDyad)] = calcAUROC(in_solo, in_dyad, idx_ply1, idx_ply2, iDyad, ['acc_' str]);

end
end

function out = computer_player_modulation(solo_perf,hc_dyad_pw_perf,dyad_pw_perf)

tmp_dyad                        = [dyad_pw_perf(:,1); dyad_pw_perf(:,2)];
tmp_agent                       = hc_dyad_pw_perf(:,2);
snr                             = solo_perf{1}.carr;

% Extract subject sequence
for iSubj = 1:size(tmp_agent,1)
    if ~isempty(tmp_agent{iSubj})
        id_hc_dyad{iSubj}           	= tmp_agent{iSubj}.id;
    else
        id_hc_dyad{iSubj}           	= 'empty';
    end
end

for iSubj = 1:size(tmp_dyad,1)
    if ~isempty(tmp_dyad{iSubj})
        id_dyad{iSubj}           	= tmp_dyad{iSubj}.id;
    else
        id_dyad{iSubj}           	= 'empty';
    end
end

% Contrast data
scnt = 0;
dcnt = 0;

for iSub = 1:length(solo_perf)
        
    idx_pc                          = cellfun(@(x) strcmp(x,solo_perf{iSub}.id),id_hc_dyad);
    idx_dy                          = cellfun(@(x) strcmp(x,solo_perf{iSub}.id),id_dyad);
    
    if sum(idx_pc) == 0 || sum(idx_dy) == 0
        continue
    end
    
    scnt = scnt+1;
    
    for iCoh = 1:length(snr)
        out.auc_acc_SC(scnt,iCoh)   	= getAUROC(solo_perf{iSub}.acc_state{iCoh},tmp_agent{idx_pc}.acc_state{iCoh});
        out.auc_ecc_SC(scnt,iCoh)     	= getAUROC(solo_perf{iSub}.ecc_state{iCoh},tmp_agent{idx_pc}.ecc_state{iCoh});
    end
    
    sessions = find(idx_dy);
    for iSession = sessions
        dcnt = dcnt+1;

        for iCoh = 1:length(snr)
            
            out.auc_acc_SH(dcnt,iCoh)   	= getAUROC(solo_perf{iSub}.acc_state{iCoh},tmp_dyad{iSession}.acc_state{iCoh});
            out.auc_ecc_SH(dcnt,iCoh)     	= getAUROC(solo_perf{iSub}.ecc_state{iCoh},tmp_dyad{iSession}.ecc_state{iCoh});
            
            out.auc_acc_HC(dcnt,iCoh)       = getAUROC(tmp_dyad{iSession}.acc_state{iCoh},tmp_agent{idx_pc}.acc_state{iCoh});
            out.auc_ecc_HC(dcnt,iCoh)     	= getAUROC(tmp_dyad{iSession}.ecc_state{iCoh},tmp_agent{idx_pc}.ecc_state{iCoh});
        end
    end
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