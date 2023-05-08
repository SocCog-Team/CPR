% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/WIP/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

% Import subject summary table
pth                         = '/Users/fschneider/Documents/CPR_psychophysics/';
x                           = readtable([pth 'Subjects_summary.xlsx']);
sbj_lst                     = x.Abbreviation;
sbj_lst(cellfun(@isempty,sbj_lst)) = [];
nSample                     = 30;
s_cnt                       = 0;
nLag                        = 150;

% Extract file names of experiments
c                           = 0;
cc                          = 0;
for iSubj = 1:length(sbj_lst)
    data_pth                = [pth sbj_lst{iSubj} '/summary/'];
    
    load([pth sbj_lst{iSubj} '/summary/' sbj_lst{iSubj} '_RT.mat'])
    mrt(iSubj)              = mean(rt.dat);
    
    if isdir(data_pth)
        cd(data_pth)
        mat_files        	= dir('*.mat');
        
        for iFile = 1:length(mat_files)
            if contains(mat_files(iFile).name,'CPRsolo')
                c               = c+1;
                fname_solo{c} 	= mat_files(iFile).name;
            elseif contains(mat_files(iFile).name,'CPRagent')
                cc              = cc+1;
                fname_agnt{cc} 	= mat_files(iFile).name;
            end
        end
    end
end

for iExp = 1:2
    
    if iExp == 1
        fname = fname_solo;
    elseif iExp == 2
        fname = fname_agnt;
    end
    
    for iSub = 1:length(sbj_lst)
        
        lg{iExp,iSub}                   = [];
        avg_sxc{iExp,iSub}              = [];
        
        fname_crop                      = cellfun(@(x) x(1:12), fname, 'UniformOutput', false);
        fidx                            = find(cellfun(@(x) contains(x,lower(sbj_lst{iSub})),fname_crop));
        cc                              = 0;
        t                               = [];
        all                             = [];
        
        for iBlock = 1:length(fidx)
            clear tbl
            tbl = load([pth sbj_lst{iSub} '/summary/' fname{fidx(iBlock)}]);
            
            if sum(unique(tbl.t.rdp_coh) < .1) > 2
                continue
            end
            
            t = [t; tbl.t];
            
            % Extract trials for correlation analysis
            for iTrl = 1:tbl.t.trl_no(end)
                tidx                    = tbl.t.trl_no == iTrl;
                tt                      = tbl.t(tidx,:);
                cc                      = cc+1;
                tmp.frme_ts{cc}         = [];
                tmp.rdp_dir{cc}         = [];
                tmp.rdp_coh{cc}     	= [];
                tmp.js_dir{cc}      	= [];
                tmp.js_str{cc}      	= [];
                tmp.refresh{cc}         = [];
                
                for iState = 1:size(tt,1)
                    tmp.frme_ts{cc}  	= [tmp.frme_ts{cc} tt.frme_ts{iState}];
                    tmp.rdp_dir{cc}  	= [tmp.rdp_dir{cc} repmat(tt.rdp_dir(iState),1,length(tt.frme_ts{iState}))];
                    tmp.rdp_coh{cc}  	= [tmp.rdp_coh{cc} repmat(tt.rdp_coh(iState),1,length(tt.frme_ts{iState}))];
                    tmp.js_dir{cc}  	= [tmp.js_dir{cc} tt.js_dir{iState}];
                    tmp.js_str{cc}  	= [tmp.js_str{cc} tt.js_str{iState}];
                    tmp.refresh{cc}     = [tmp.refresh{cc} median(diff(tmp.frme_ts{cc}))];
                end
            end
            
            [cr,ps]                     = CPR_correlation_analysis_WIP(tmp, nLag, false);
            avg_sxc{iExp,iSub}          = [avg_sxc{iExp,iSub}; mean(cr.sxc)];
            max_pos                  	= find(mean(cr.sxc) == max(mean(cr.sxc)),1,'first');
            lag                         = ((max_pos - nLag) * median(tmp.refresh{cc})) / 1e3;
            lg{iExp,iSub}               = [lg{iExp,iSub} lag];
            id{iSub}                    = sbj_lst{iSub} ;
        end
        
        if isempty(t)
            continue
        end
        
        snr                     = unique(t.rdp_coh);
        
        %         % Target score
        %         clear score score_hi score_coh
        %         c                           = 0;
        %         for iState = 1:size(t.trg_ts,1)
        %             for iTrg = 1:length(t.trg_ts{iState})
        %                 c                   = c +1;
        %                 score_cum(c)        = t.trg_score{iState}(iTrg);
        %                 score_coh(c)        = t.rdp_coh(iState);
        %                 score_hi(c)         = t.trg_hit{iState}(iTrg);
        %             end
        %         end
        %
        %         tmp_score                   = score_cum(~isnan(score_cum));
        %         tscore(~isnan(score_cum))   = [0 diff(tmp_score)];
        %         tscore(tscore < 0)          = nan;
        %         highscore(iExp,iSub)        = tmp_score(end); %% across all sessions?
        
        for iCoh = 1:length(snr)
            clear cIdx rdp_dir t.js_dir
            
            cIdx = t.rdp_coh == snr(iCoh);
            cIdx(cellfun(@length,t.js_str) < 100) = false;
            
            % Hit rate
            nhi                     = sum(cellfun(@sum,t.trg_hit(cIdx)));
            ntrg                    = sum(cellfun(@numel,t.trg_hit(cIdx)));
            hir(iExp,iSub,iCoh)     = nhi / ntrg;
            
            % Target score
            %             trg_score(iExp,iSub,iCoh) 	= nanmean(tscore(score_coh  == snr(iCoh) & score_hi == true));
            
            % Joystick displacement
            mstr(iExp,iSub, iCoh)       = nanmedian(cellfun(@(x) nanmedian(x(end-nSample:end)), t.js_str(cIdx)));
            
            % Joystick accuracy
            rdp_dir                 = t.rdp_dir(cIdx);
            js_dir                  = t.js_dir(cIdx);
            
            for iState = 1:length(rdp_dir)
                clear js_dev
                js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iState}(end-nSample:end)),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
                js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));                           % Joystick accuracy
            end
            
            carr(iExp,iSub,snr == snr(iCoh))   	= snr(iCoh);
            macc(iExp,iSub,snr == snr(iCoh))    = nanmedian(js_acc);
            %             sxc{iExp,iSub,snr == snr(iCoh)}   	= cr.sxc(cr.coh == snr(iCoh),:);
            %             mR(iExp,iSub,snr == snr(iCoh))     	= mean(cr.maxR(cr.coh == snr(iCoh)));
            %             mcc(iExp,iSub,snr == snr(iCoh))   	= mean(cr.cc(cr.coh == snr(iCoh)));
            %             mpk(iExp,iSub,snr == snr(iCoh))   	= mean(cr.posPk(cr.coh == snr(iCoh)));
            %             spk(iExp,iSub,snr == snr(iCoh))    	= std(cr.posPk(cr.coh == snr(iCoh)));
        end
    end
end

%% Average lag

figure
df                          = cellfun(@nanmean,lg)';
mrt(sum(isnan(df),2)> 0)    = [];
df(sum(isnan(df),2)> 0,:)   = [];
% cl                          = cbrewer('qual', 'Paired', size(df,1), 'PCHIP');
cl                          = jet(size(df,1));

f = figure; hold on
vl                          = violinplot([mrt' df]);
for v = 1:3
vl(v).ShowData = 0;
vl(v).WhiskerPlot.Color     = 'none';
vl(v).BoxPlot.FaceColor     = 'none';
vl(v).BoxPlot.EdgeColor     = 'none';
vl(v).MedianPlot.Marker     = 'none';
vl(v).ViolinColor           = [.5 .5 .5];
end

for iSubj = 1:size(df,1)    
    plot([1 2 3],[mrt(iSubj) df(iSubj,:)], 'Color', [cl(iSubj,:) .5], 'LineWidth',2)
    scatter(1, mrt(iSubj),'filled', 'CData', cl(iSubj,:), 'SizeData', 100, 'MarkerFaceAlpha', 1);
    scatter(2, df(iSubj,1),'filled', 'CData', cl(iSubj,:), 'SizeData', 100, 'MarkerFaceAlpha', 1);
    scatter(3, df(iSubj,2),'filled', 'CData', cl(iSubj,:), 'SizeData', 100, 'MarkerFaceAlpha', 1)
end

ax                          = gca;
ax.XTick                    = [1 2 3];
ax.XTickLabel               = {'Reaction Time','CPR_solo', 'CPR_computer'};
ax.XLim                     = [.5 3.5];
ax.YTick                    = [200 450 700];
ax.YLim                     = [200 700];
ax.YLabel.String            = 'Response lag [ms]';
ax.TickLabelInterpreter     = 'none';
ax.FontSize                 = 28;
ax.XTickLabelRotation       = 13;
ax.FontWeight               = 'bold';

print(f, ['/Users/fschneider/Desktop/lag_diff'], '-r300', '-dpng');

%%

for iPlot = 1:3
    clear dat_df
    
    if iPlot == 1
        dat         = macc;
        ystr        = {'Accuracy difference'; '[Dyad.computer - Solo]'};
    elseif iPlot == 2
        dat         = mstr;
        ystr        = {'Eccentricity difference'; '[Dyad.computer - Solo]'};
    elseif iPlot == 3
        dat         = hir;
        ystr        = {'Hit rate difference'; '[Dyad.computer - Solo]'};
    end
    
    c               = 0;
    for iSub = 1:size(dat,2)
        
        if sum(dat(2,iSub,:)) == 0
            continue
        end
        
        c                           = c+1;
        dat_df(c,:)                 = squeeze(dat(2,iSub,:)-dat(1,iSub,:));
    end
    
    f = figure;
    
    ln = line([0 8],[0 0],'LineWidth',2,'LineStyle',':','Color',[.3 .3 .3]);
    
    vl                           	= violinplot(dat_df);
    for iV = 1:size(vl,2)
        vl(iV).ViolinColor       	= [.5 .5 .5];
        vl(iV).ShowData            	= 0;
        vl(iV).WhiskerPlot.Color     = 'none';
        vl(iV).BoxPlot.FaceColor     = 'none';
        vl(iV).BoxPlot.EdgeColor     = 'none';
        vl(iV).MedianPlot.Marker     = 'none';
        vl(iV).ViolinColor           = [.5 .5 .5];
        
        for iDat = 1:size(dat_df,1)
            sc = scatter(iV, dat_df(iDat,iV),'filled', 'CData', cl(iDat,:), 'SizeData', 100, 'MarkerFaceAlpha', .6);
        end
    end
    
    ax                              = gca;
    ax.FontSize                     = 22;
    ax.YLabel.String                = ystr;
    ax.XLabel.String                = 'Coherence [%]';
    ax.XTickLabel                   = round(snr,2)*100;
    ax.TickLabelInterpreter         = 'none';

    print(f, ['/Users/fschneider/Desktop/agnt_diff' num2str(iPlot)], '-r300', '-dpng');

end


%%



% %% PLOT DIFF
%
% lb = {'Hit Rate',...
%     'Hit Score',...
%     'Avg JS Displacement',...
%     'Avg JS Accuracy',...
%     'Avg CirCorr Coeff',...
%     'Avg XCorr Peak Position',...
%     'Std XCorr Peaks',...
%     'Avg XCorr Coeff'};
%
% f                       = figure('units','normalized','outerposition',[0 0 1 1]);
% hold on
%
% for iDat = 1:8
%
%     if iDat == 1
%         d = hir;
%     elseif iDat == 2
%         d = trg_score;
%     elseif iDat == 3
%         d = mstr;
%     elseif iDat == 4
%         d = macc;
%     elseif iDat == 5
%         d = mcc;
%     elseif iDat == 6
%         d = mpk;
%     elseif iDat == 7
%         d = spk;
%     elseif iDat == 8
%         d = mR;
%     end
%
%     xtick                   = [1 9 length(snr)];
%     ax                      = subplot(2,4,iDat);
%     im                      = imagesc(d(:,:,1)-d(:,:,2));
%     ax.Title.String         = lb{iDat};
%     ax.Title.Interpreter    = 'none';
%     ax.FontSize             = 16;
%     ax.YTick                = 1:size(fname,2);
%     ax.YTickLabel           = cellfun(@(x) x(10:12),fname,'UniformOutput',false);
%     ax.XTick                = xtick;
%     ax.XTickLabel           = snr(xtick);
%     ax.XLabel.String        = 'Coherence';
%     cb                      = colorbar;
%
%     if iDat == 1
%         caxis([-.25 .25])
%     elseif iDat == 2
%         caxis([-.301 .301])
%     elseif iDat == 3
%         caxis([-.35 .35])
%     elseif iDat == 4
%         caxis([-.1 .1])
%     elseif iDat == 5
%         caxis([-.5 .5])
%     elseif iDat == 6
%         caxis([-60 60])
%     elseif iDat == 7
%         caxis([-50 50])
%     elseif iDat == 8
%         caxis([-1000 1000])
%     end
%
%     colormap(cbrewer('div', 'RdGy', 256, 'PCHIP'))
%     axis square;
%
% end
%
%
% print(f, ['/Users/fschneider/Desktop/solo_minus_agnt'], '-r300', '-dpng');
%
% %% PLOT LAG
%
% for j = 1:size(fname,2)
%     solo                = [];
%     agnt                = [];
%
%     for i = 1:length(snr)
%         if size(sxc{j,i,1},1) == 1
%             solo        = [solo; sxc{j,i,1}];
%         else
%             solo       	= [solo; mean(sxc{j,i,1})];
%         end
%
%         agnt            = [agnt; mean(sxc{j,i,2})];
%     end
%
%     all{j}              = solo;
%     ms(j,:)             = mean(solo);
%     ma(j,:)             = mean(agnt);
%
%     p(j,1)              = find(mean(solo) == max(mean(solo)));
%     p(j,2)              = find(mean(agnt) == max(mean(agnt)));
% end
%
% c                       = cbrewer('qual', 'Dark2', size(fname,2), 'PCHIP');
% alph                    = .75;
%
% f                       = figure('units','normalized','outerposition',[0 0 1 1]);
% ax                      = subplot(2,3,1);
% im                      = imagesc(solo);
% ax.YTick                = 1:length(snr);
% ax.YTickLabel           = snr;
% ax.YLabel.String        = 'Coherence';
% ax.XTick                = [1 151 301];
% ax.XTickLabel           = round([-150 0 150] * 8.333);
% ax.XLabel.String        = 'Lag [ms]';
% ax.Title.String         = ['Exmpl Subj: ' fname_solo{j}(10:12) '_solo'];
% ax.FontSize             = 20;
% ax.Title.Interpreter    = 'none';
% cb                      = colorbar;
% cb.Label.String         = 'Xcorr coeff';
%
% caxis([0 750])
% colormap(gray)
%
% ax                      = subplot(2,3,2);
% hold on
% for i = 1:size(ms,1)
%     pl                  = plot(ms(i,:),'LineWidth',2,'Color',[c(i,:) alph]);
% end
% ax.XTick                = [1 151 301];
% ax.XTickLabel           = round([-150 0 150] * 8.333);
% ax.XLabel.String        = 'Lag [ms]';
% ax.YLabel.String        = 'XCorr Coeff';
% ax.FontSize         	= 20;
% ax.Title.String         = 'Average/Subj: Solo';
% ax.Position(1)          = ax.Position(1)+.02;
%
% ax                      = subplot(2,3,4);
% im                      = imagesc(agnt);
% ax.YTick                = 1:length(snr);
% ax.YTickLabel           = snr;
% ax.YLabel.String        = 'Coherence';
% ax.XTick                = [1 151 301];
% ax.XTickLabel           = round([-150 0 150] * 8.333);
% ax.XLabel.String        = 'Lag [ms]';
% ax.Title.String         = ['Exmpl Subj: ' fname_agent{j}(10:12) '_computer'];
% ax.Title.Interpreter    = 'none';
% ax.FontSize             = 20;
% cb                      = colorbar;
% cb.Label.String         = 'Xcorr coeff';
%
% caxis([0 750])
% colormap(gray)
%
% ax                      = subplot(2,3,5);
% hold on
% for i = 1:size(ma,1)
%     pl                  = plot(ma(i,:),'LineWidth',2,'Color',[c(i,:) alph]);
% end
% ax.XTick                = [1 151 301];
% ax.XTickLabel           = round([-150 0 150] * 8.333);
% ax.XLabel.String        = 'Lag [ms]';
% ax.YLabel.String        = 'XCorr Coeff';
% ax.FontSize             = 20;
% ax.Title.String       	= 'Average/Subj: Computer';
% ax.Position(1)          = ax.Position(1)+.02;
%
% dat                     = ((p(:,1)-p(:,2))*8.3333);
% ax                  	= subplot(2,3,[3 6]);
% bx                      = boxplot(dat, 'Colors', 'k');
%
% set(bx(end,:),'Visible','off')
% set(bx, {'linew'},{3})
% hold on
%
% for i = 1:length(dat)
%     sc               	= scatter(1,((p(i,1)-p(i,2))*8.3333));
%     sc.SizeData         = 100;
%     sc.MarkerFaceColor  = c(i,:);
%     sc.MarkerFaceAlpha  = 1;
%     sc.MarkerEdgeColor  = 'none';
% end
%
% ax.YLabel.String        = 'Peak difference [ms]';
% ax.XTick                = [];
% ax.Box                  = 'off';
% ax.FontSize             = 20;
%
% lg                      = legend(cellfun(@(x) x(10:12),fname,'UniformOutput',false));
% lg.Location             = 'east';
% lg.Box                  = 'off';
%
% print(f, ['/Users/fschneider/Desktop/Pk_diff'], '-r300', '-dpng');
%
% %% Solo plot
%
% f = figure('units','normalized','outerposition',[0 0 1 1]);
% for s = 1:length(all)
%     ax                      = subplot(2,6,s);
%     im                      = imagesc(all{s}./max(max(all{s})));
%     cm                      = colormap(gray(256));
%     ax.XTick                = [1 151 301];
%     ax.XTickLabel           = round([-150 0 150] * 8.333);
%     ax.XLabel.String        = 'Lag [ms]';
%     ax.YTick                = 1:length(snr);
%     ax.YTickLabel           = snr;
%     ax.YLabel.String        = 'Coherence';
%     ax.Title.String         = [fname_solo{s}(10:12)];
%     ax.Title.Interpreter    = 'none';
%     ax.FontSize             = 16;
% end
% print(f, ['/Users/fschneider/Desktop/xCorr_coh_subj'], '-r300', '-dpng');
%
% %%
%
% f                       = figure;
% dat                     = (p(:,1)*8.3333);
% bx                      = boxplot(dat, 'Colors', 'k');
%
% set(bx(end,:),'Visible','off')
% set(bx, {'linew'},{3})
% hold on
%
% for i = 1:length(dat)
%     sc               	= scatter(1+randi([-3 3])/100,(p(i,1)*8.3333));
%     sc.SizeData         = 100;
%     sc.MarkerFaceColor  = c(i,:);
%     sc.MarkerFaceAlpha  = 1;
%     sc.MarkerEdgeColor  = 'none';
% end
%
% ax                      = gca;
% ax.YLabel.String        = 'Response delay [ms]';
% ax.XTick                = [];
% ax.Box                  = 'off';
% ax.FontSize             = 20;
%
% lg                      = legend(cellfun(@(x) x(10:12),fname,'UniformOutput',false));
% lg.Location             = 'east';
% lg.Box                  = 'off';
% print(f, ['/Users/fschneider/Desktop/Resp_delay'], '-r300', '-dpng');
%
% %%
%
% pool = [];
% for ii = 1:length(all)
%     pool = [pool all{ii}];
%     df(:,ii) = max(all{ii},[],2) ./ mean(all{ii},2);
% end
%
% f                       = figure;
% ax                      = gca;
% im                      = imagesc(df);
% cm                      = colormap(gray(256));
% ax.XLabel.String        = 'Subject';
% ax.YTick                = 1:length(snr);
% ax.YTickLabel           = snr;
% ax.YLabel.String        = 'Coherence';
% ax.XTick                = 1:length(fname);
% ax.XTickLabel           = cellfun(@(x) x(10:12),fname,'UniformOutput',false);
% ax.Title.String         = 'XCorr coeff: Max / Mean';
% ax.Title.Interpreter    = 'none';
% ax.FontSize             = 16;
% cb                      = colorbar;
% cb.Label.String         = 'Max/Mean ratio';