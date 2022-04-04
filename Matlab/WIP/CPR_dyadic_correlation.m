addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/Documents/MATLAB/cbrewer

clear all
close all

pth                     = '/Users/fschneider/Documents/CPR_psychophysics/';
fname                   = 'Subjects_summary.xlsx';
xd                      = readtable([pth fname],'Sheet','Dyads');
dte_format              = 'yyyymmdd';
% col                     = cbrewer('div', 'Spectral', length(snr), 'PCHIP');
alph                    = .2;

for iDyad = 1:13%xd.Dyad'
    % Specify directories
    summ_pth            = [pth 'Dyad' num2str(iDyad) '/summary/'];
    cd(summ_pth)
    
    % Identify players
    p1                  = xd.Block1_PSY4{iDyad};
    p2                  = xd.Block1_PSY3{iDyad};
    
    % Import tables for each block
    tmp11               = load([datestr(xd.Date(iDyad),dte_format) '_' p1 '_CPRdyadic_block1_tbl.mat']);
    tmp21               = load([datestr(xd.Date(iDyad),dte_format) '_' p2 '_CPRdyadic_block1_tbl.mat']);
    tmp12               = load([datestr(xd.Date(iDyad),dte_format) '_' p1 '_CPRdyadic_block2_tbl.mat']);
    tmp22               = load([datestr(xd.Date(iDyad),dte_format) '_' p2 '_CPRdyadic_block2_tbl.mat']);  
   
    % Concatenate tables for each player
    t1                  = [tmp11.t;tmp12.t];
    t2                  = [tmp21.t;tmp22.t];
    
    snr              	= unique(t1.ss_coh);
    nSamples           	= 30;
    cl                	= linspace(0,1,length(snr));
    
%     [acc11all,str11all] = getAcc(tmp11.t.js_dir,tmp11.t.js_str,tmp11.t.rdp_dir,30);
%     [acc12all,str12all] = getAcc(tmp12.t.js_dir,tmp12.t.js_str,tmp12.t.rdp_dir,30);
%     [acc21all,str21all] = getAcc(tmp21.t.js_dir,tmp21.t.js_str,tmp21.t.rdp_dir,30);
%     [acc22all,str22all] = getAcc(tmp22.t.js_dir,tmp22.t.js_str,tmp22.t.rdp_dir,30);
%     [h1s(iDyad),p1s(iDyad)]=ranksum(str11all,str12all);
%     [h1a(iDyad),p1a(iDyad)]=ranksum(acc11all,acc12all);
%     [h(iDyad),p(iDyad)]=ranksum(str21all,str22all);
%     [h(iDyad),p(iDyad)]=ranksum(acc21all,acc22all);
%     
%     for iCoh = 1:length(snr)
%         clear cIdx
%         
%         cIdx1                = tmp11.t.ss_coh == snr(iCoh);
%         cIdx2                = tmp12.t.ss_coh == snr(iCoh);
%         cIdx1(cellfun(@length,tmp11.t.js_str) < 100) = false;
%         cIdx2(cellfun(@length,tmp21.t.js_str) < 100) = false;
% 
%         % Check consistency
%         [acc11{iCoh},str11{iCoh}] 	= getAcc(tmp11.t.js_dir(cIdx1),tmp11.t.js_str(cIdx1),tmp11.t.rdp_dir(cIdx1),30);
%         [acc12{iCoh},str12{iCoh}]   = getAcc(tmp12.t.js_dir(cIdx2),tmp12.t.js_str(cIdx2),tmp12.t.rdp_dir(cIdx2),30);
%         [acc21{iCoh},str21{iCoh}]   = getAcc(tmp21.t.js_dir(cIdx1),tmp21.t.js_str(cIdx1),tmp21.t.rdp_dir(cIdx1),30);
%         [acc22{iCoh},str22{iCoh}]   = getAcc(tmp22.t.js_dir(cIdx2),tmp22.t.js_str(cIdx2),tmp22.t.rdp_dir(cIdx2),30);
%     end
    
%     figure;hold on
%     plot(cellfun(@median, str11) - cellfun(@median, str12))
%     plot(cellfun(@median, str21) - cellfun(@median, str22))
%     plot(cellfun(@median, acc11) - cellfun(@median, acc12))
%     plot(cellfun(@median, acc21) - cellfun(@median, acc22))
    
% pa1(iDyad,:) = cellfun(@(x,y) ttest2(x,y),acc11,acc12) < .05;
% pa2(iDyad,:) = cellfun(@(x,y) ttest2(x,y),acc21,acc22) < .05;
% ps1(iDyad,:) = cellfun(@(x,y) ttest2(x,y),str11,str12) < .05;
% ps2(iDyad,:) = cellfun(@(x,y) ttest2(x,y),str21,str22) < .05;

    for iCoh = 1:length(snr)
        clear cIdx rdp_dir t.js_dir
        
        cIdx                = t1.ss_coh == snr(iCoh);
        cIdx(cellfun(@length,t1.js_str) < 100) = false;

        % Joystick displacement
        js1_str             = t1.js_str(cIdx);
        js2_str             = t2.js_str(cIdx);
        
        % Joystick accuracy
        rdp_dir             = t1.rdp_dir(cIdx);
        js1_dir             = t1.js_dir(cIdx);
        js2_dir             = t2.js_dir(cIdx);
        
        [acc1{iCoh},str1{iCoh}]     = getAcc(js1_dir,js1_str,rdp_dir,nSamples);
        [acc2{iCoh},str2{iCoh}]     = getAcc(js2_dir,js2_str,rdp_dir,nSamples);
  
        tmp1 = corrcoef(acc1{iCoh},str1{iCoh});
        R1(iDyad,iCoh) = tmp1(1,2);
        tmp2 = corrcoef(acc2{iCoh},str2{iCoh});
        R2(iDyad,iCoh) = tmp2(1,2);
        
        tmp12 = corrcoef(acc1{iCoh},str2{iCoh});
        R12(iDyad,iCoh) = tmp12(1,2);
        
        tmp21 = corrcoef(acc2{iCoh},str1{iCoh});
        R21(iDyad,iCoh) = tmp21(1,2);
    end
    
    
    figure('units','normalized','outerposition',[0 0 1 1])
    ax = subplot(2,4,1);hold on
    for iCoh = 1:length(snr)
        sc = scatter(acc1{iCoh},str1{iCoh},'MarkerFaceColor', [cl(iCoh) 0 0],'MarkerEdgeColor', [cl(iCoh) 0 0]);
        sc.MarkerFaceAlpha = alph;
        sc.MarkerEdgeAlpha = alph;
        coefs1(iCoh,:) = polyfit(acc1{iCoh}, str1{iCoh}, 1);
        plot([0 1],polyval(coefs1(iCoh,:),0:1),'Color',[cl(iCoh) 0 0 alph*3],'LineWidth',2)
    end
    ax.XLabel.String = 'Player1 JS Accuracy';
    ax.YLabel.String = 'Player1 JS Displacement';
    ax.Title.String = 'JS response_within_sbj';
    ax.FontSize = 16;
    ax.Title.Interpreter = 'none';
    
    ax = subplot(2,4,2);hold on
    for iCoh = 1:length(snr)
        sc = scatter(acc2{iCoh},str2{iCoh},'MarkerFaceColor', [cl(iCoh) 0 0],'MarkerEdgeColor', [cl(iCoh) 0 0]);
        sc.MarkerFaceAlpha = alph;
        sc.MarkerEdgeAlpha = alph;
        coefs2(iCoh,:) = polyfit(acc2{iCoh}, str2{iCoh}, 1);
        plot([0 1],polyval(coefs2(iCoh,:),0:1),'Color',[cl(iCoh) 0 0 alph*3],'LineWidth',2)
    end
    ax.XLabel.String = 'Player2 JS Accuracy';
    ax.YLabel.String = 'Player2 JS Displacement';
    ax.Title.String = 'JS response_within_sbj';
    ax.FontSize = 16;
    ax.Title.Interpreter = 'none';
        
    ax = subplot(2,4,5);hold on
    for iCoh = 1:length(snr)
        sc = scatter(acc1{iCoh},str2{iCoh},'MarkerFaceColor', [cl(iCoh) 0 0],'MarkerEdgeColor', [cl(iCoh) 0 0]);
        sc.MarkerFaceAlpha = alph;
        sc.MarkerEdgeAlpha = alph;
        coefs12(iCoh,:) = polyfit(acc1{iCoh}, str2{iCoh}, 1);
        plot([0 1],polyval(coefs12(iCoh,:),0:1),'Color',[cl(iCoh) 0 0 alph*3],'LineWidth',2)
    end
    ax.XLabel.String = 'Player1 JS Accuracy';
    ax.YLabel.String = 'Player2 JS Displacement';
    ax.Title.String = 'JS response_between_sbj';
    ax.FontSize = 16;
    ax.Title.Interpreter = 'none';
    
    ax = subplot(2,4,6);hold on
    for iCoh = 1:length(snr)
        sc = scatter(acc2{iCoh},str1{iCoh},'MarkerFaceColor', [cl(iCoh) 0 0],'MarkerEdgeColor', [cl(iCoh) 0 0]);
        sc.MarkerFaceAlpha = alph;
        sc.MarkerEdgeAlpha = alph;
        coefs21(iCoh,:) = polyfit(acc2{iCoh}, str1{iCoh}, 1);
        plot([0 1],polyval(coefs21(iCoh,:),0:1),'Color',[cl(iCoh) 0 0 alph*3],'LineWidth',2)
    end
    ax.XLabel.String = 'Player2 JS Accuracy';
    ax.YLabel.String = 'Player1 JS Displacement';
    ax.Title.String = 'JS response_between_sbj';
    ax.FontSize = 16;
    ax.Title.Interpreter = 'none';
   
    ax = subplot(2,4,3);hold on
    p1=plot(coefs1(:,1),'Color', [1 0 0], 'Linewidth',2);
    p2 = plot(coefs2(:,1),'Color', [0 0 0], 'Linewidth',2);
    ax.XLabel.String = 'Coherence';
    ax.XTick= 1:length(snr);
    ax.XTickLabel= snr;
    ax.YLabel.String = 'Slope Accuracy - Displacement';
    ax.Title.String = 'Slope_within';
    ax.FontSize = 16;
    ax.Title.Interpreter = 'none';
    lg = legend([p1,p2], 'Player1', 'Player2');
    
    ax = subplot(2,4,4);hold on
    pr1 = plot(R1(iDyad,:)-R12(iDyad,:),'Color', [1 0 0], 'Linewidth',2);
    pr2 = plot(R2(iDyad,:)-R21(iDyad,:),'Color', [0 0 0], 'Linewidth',2);
    ax.XLabel.String = 'Coherence';
    ax.XTick= 1:length(snr);
    ax.XTickLabel= snr;
    ax.YLabel.String = 'Diff Corr_Coef';
    ax.Title.String = 'Correlation: within minus between';
    ax.FontSize = 16;
    ax.Title.Interpreter = 'none';
    lg = legend([pr1,pr2], 'Player1', 'Player2');
    
    ax = subplot(2,4,7);hold on
    p12 = plot(coefs12(:,1),'Color', [1 0 0], 'Linewidth',2);
    p21 = plot(coefs21(:,1),'Color', [0 0 0], 'Linewidth',2);
    ax.XLabel.String = 'Coherence';
    ax.XTick= 1:length(snr);
    ax.XTickLabel= snr;
    ax.YLabel.String = 'Slope Accuracy - Displacement';
    ax.Title.String = 'Slope_between';
    ax.FontSize = 16;
    ax.Title.Interpreter = 'none';
    lg = legend([p12,p21], 'Player1_acc vs Player2_dis', 'Player2_acc vs Player1_dis');
    lg.Interpreter = 'none';

end

%% PLOT 

rr = [R1];
figure;hold on
boxplot(rr)
for i=1:17
    sc(i) = scatter(repmat(i,length(rr(:,i)),1),rr(:,i));
    sc(i).MarkerFaceColor = [0 0 0];
    sc(i).MarkerEdgeColor = 'none';
end
ax = gca;
ax.XLabel.String = 'Coherence';
ax.XTick= 1:length(snr);
ax.XTickLabel= snr;
ax.YLabel.String = 'Correlation coef';
ax.Title.String = 'Player1_correlation_acc_displ';
ax.Title.Interpreter = 'none';
ax.FontSize = 16;

rrr = [R12];
figure;hold on
boxplot(rrr)
for i=1:17
    sc(i) = scatter(repmat(i,length(rrr(:,i)),1),rrr(:,i));
    sc(i).MarkerFaceColor = [0 0 0];
    sc(i).MarkerEdgeColor = 'none';
end
ax = gca;
ax.XLabel.String = 'Coherence';
ax.XTick= 1:length(snr);
ax.XTickLabel= snr;
ax.YLabel.String = 'Correlation coef';
ax.Title.String = 'Correlation_player1_acc_player2_displ';
ax.Title.Interpreter = 'none';
ax.FontSize = 16;

figure;hold on
for i = 1:size(R1,1)
plot(R1(i,:)-R12(i,:),'Color', [0 0 0 .5], 'Linewidth',2)
end

for i = 1:size(R2,1)
plot(R2(i,:)-R21(i,:),'Color', [0 0 0 .5], 'Linewidth',2)
end

ax = gca;
ax.XLabel.String = 'Coherence';
ax.XTick= 1:length(snr);
ax.XTickLabel= snr;
ax.YLabel.String = 'Correlation coef';
ax.Title.String = 'Difference_correlation: Within-Between';
ax.Title.Interpreter = 'none';
ax.FontSize = 16;

plot(median([R1-R12;R2-R21]),'Color', [1 0 0], 'Linewidth',2)

%% Helper

function [out_acc,out_str] = getAcc(js_dir,js_str,rdp_dir,nSamples)
for iState = 1:length(rdp_dir)
    clear js_dev
    js_dev           	= rad2deg(circ_dist(deg2rad(js_dir{iState}),deg2rad(rdp_dir(iState))));         % Minimum RDP-Joystick difference
%     out_acc(iState)  	= nanmean(abs(1 - abs(js_dev(end-nSamples:end)) / 180));                        % Joystick accuracy
%     out_str(iState)  	= nanmean(js_str{iState}(end-nSamples:end));
    out_acc(iState)  	= nanmean(abs(1 - abs(js_dev) / 180));                       
    out_str(iState)  	= nanmean(js_str{iState});
end
end

% FILTER
% only hit states
% only right after change
% -> reduce variability
% FINDING: Differences between subjects
% What produces this type of effect? CLUSTER SUBJECTS! Confusion/Detrimental performance vs Integration of additional
% information - better
% Correelation between str and acc fro different classes
% look also at agent vs dyadic cluster according to psychometric curve -> also modulated by signal2noise
% follower vs leader vs independent - might change with SNR
% decision times for each subject

