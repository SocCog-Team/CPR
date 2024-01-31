load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/hh_dyad_performance.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUROC solo - dyadic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:size(dyad_perf,2)
    if isempty(dyad_perf{j})
        ply_id{j} = nan;
    else
        ply_id{j} = dyad_perf{j}.id;
    end
end

cnt = 0;
for iSubj = 1:size(dyad_perf,2)
    
    d.sp                  	= solo_perf{iSubj};                             % Performance data
    d.sc                  	= solo_cr{iSubj};                               % Correlation data
    
    sIdx                    = cellfun(@(x) strcmp(x,solo_perf{iSubj}.id),ply_id);
    
    if sum(sIdx) == 0
        continue
    end
    
    d.dp                  	= dyad_perf{sIdx};
    d.dc                  	= dyad_cr{sIdx};
         
    if isempty(d.dc)
        continue
    end
    
    cnt = cnt+1;

    for iCoh = 1:length(snr)
        bl_ecc(cnt,iCoh)     	= d.sp.mecc_state(iCoh);
        bl_acc(cnt,iCoh)     	= d.sp.macc_trg(iCoh);
        bl_scr(cnt,iCoh)     	= d.sp.trg_mscore(iCoh);
        bl_hir(cnt,iCoh)     	= d.sp.hir(iCoh);

        auc_acc(cnt,iCoh)     	= getAUROC(d.sp.acc_trg{iCoh},d.dp.acc_trg{iCoh});
        auc_ecc(cnt,iCoh)     	= getAUROC(d.sp.ecc_state{iCoh},d.dp.ecc_state{iCoh});
        auc_score(cnt,iCoh)    	= getAUROC(d.sp.trg_score{iCoh},d.dp.trg_score{iCoh});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
close all 

figure(1);hold on
figure(2);hold on

coh_col                     = cool(length(snr));

for iPlot = 1:4
    clear xdat
    if iPlot == 1
        xdat = bl_ecc;
        str = 'Eccentricity';
    elseif iPlot == 2
        xdat = bl_acc;
        str = 'Accuracy'; 
     elseif iPlot == 3
        xdat = bl_scr;
        str = 'Hit score';  
            elseif iPlot == 4
        xdat = bl_hir;
        str = 'Hit rate'; 
    end
    
    figure(1)
    ax = subplot(2,2,iPlot); hold on
    ln = line([0 1],[.5 .5], 'Color','k', 'LineStyle', ':', 'LineWidth',1);
    
    for iCoh = 1:7
        scatter(xdat(:,iCoh),auc_ecc(:,iCoh),'MarkerFaceColor', coh_col(iCoh,:),'MarkerEdgeColor','none')
    end
    
    xlabel({str ' [Solo]'});
    ylabel({'Eccentricity AUC'; '[Solo vs Dyadic]'});
    xlim([0 1])
    ylim([0 1])
    set(ax,'fontsize', 14)
    
    figure(2)
    ax = subplot(2,2,iPlot); hold on
    
    for iCoh = 1:7
        [~,idx_solo] = sort(xdat(:,iCoh));
        [~,idx_auc] = sort(auc_ecc(:,iCoh));
        mat = [(1:34)' idx_solo idx_auc];
        for j = 1:34
            xx(j) = find(mat(:,2) == j);
            yy(j) = find(mat(:,3) == j);
        end
        scatter(xx,yy,'MarkerFaceColor',coh_col(iCoh,:),'MarkerEdgeColor','none')
    end
    
    xlabel({str ' [Solo rank]'});
    ylabel('Rank AUC');
    set(ax,'fontsize', 14)
end

function [out] = getAUROC(in1, in2)

if size(in1,1) == 1
    in1         = in1';
    in2         = in2';
end

lab          	= [zeros(length(in1),1); ones(length(in2),1)];
[~,~,~,out]     = perfcurve(lab,[in1; in2],1);

end

