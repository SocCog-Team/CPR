function [summ,crr] = CPR_plot_performance_summary(d,idx,tbl,rew_str,sbj,fname)

for iSubj = 1:size(tbl,2)
    
    f                           = figure('Units', 'normalized', 'Position', [0 0 .8 1]); set(gcf,'color', [1 1 1]);
    t                           = [];
    
    if iSubj == 1
        t                       = tbl{1};
    elseif iSubj == 2
        t                       = tbl{2};
    end
        
    trg_tbl                   	= t(logical(t.trg_shown),:);
    thi                         = [trg_tbl.trg_hit{:}];
    tts                         = [trg_tbl.trg_ts{:}];
    HIidx                       = thi(~isnan(tts));
    HIr                         = sum(HIidx) / length(HIidx);
    MIr                         = sum(~HIidx) / length(HIidx);
    trg_ts                      = tts(~isnan(tts));
    trg_ts                      = (trg_ts - t.frme_ts{1}(1)) ./ 1e6;
    rew                         = d.value(idx.reward);
    exIdx                       = strcmp(cellfun(@class, rew, 'UniformOutput', false),'double');
    rew                         = cell2mat(rew(exIdx));
    rew_ts                      = d.time(idx.reward);
    rew_ts                      = rew_ts(exIdx);
    rew_ts                      = (rew_ts - t.frme_ts{1}(1)) ./ 1e6;
    
    %%% Performance pie chart
    s                           = subplot(3,2,1);
    p                           = pie([HIr,MIr], [1,1]);
    pPatch                      = findobj(p,'Type','patch');
    pText                       = findobj(p,'Type','text');
    percentValues               = get(pText,'String');
    txt                         = {'HIr: ';'MIr: '};
    combinedtxt                 = strcat(txt,percentValues);
    s.Title.String              = {'Overall', 'Performance'};
    s.Title.FontSize            = 14;
    s.Title.Position(1)         = -1;
    s.Position(1)               = 0;
    col                         = {[1 0 0],[0 0 0]};
    
    for i = 1:size(pText,1)
        pText(i).String         = combinedtxt(i);
        pText(i).FontSize       = 14;
        pPatch(i).FaceColor     = col{i};
    end
    
    % Condition-wise performance
    trg_shown                   = cellfun(@(x) (~isnan(x)),trg_tbl.trg_ts,'uni',false);
    tcoh_tmp                    = cellfun(@times,trg_shown,num2cell(trg_tbl.rdp_coh),'uni',false);
    tcoh                        = [tcoh_tmp{:}];
    tcoh(isnan(tts))            = [];
    clvl                        = unique(trg_tbl.rdp_coh);
    for iCoh = 1:length(clvl)
        cindx               	= tcoh == clvl(iCoh);
        hir(iCoh)               = sum(HIidx(cindx)) / length(HIidx(cindx));
        mir(iCoh)               = sum(~HIidx(cindx)) / length(HIidx(cindx));
    end
    
    ax                          = axes('Position',[.345 .75 .15 .15]); hold on
    bp                          = bar([hir',mir'], 'stacked');
    cl                          = [1 0];
    
    for iOutc = 1:2
        bp(iOutc).FaceColor    	= [cl(iOutc) 0 0];
        bp(iOutc).EdgeColor     = [cl(iOutc) 0 0];
    end
    
    ax.XTick                    = [1:length(clvl)];
    ax.XLim                     = [0.5 length(clvl)+.5];
    ax.YLabel.String            = 'Rate';
    ax.XLabel.String            = 'Coherence level';
    ax.XTickLabel               = {round(clvl,2)};
    ax.FontSize                 = 12;
    
    %%% Performance over time
    col                         = {[1 0 0], [0 0 0], [.5 .5 .5]};             	% Color specs
    
    s                           = subplot(3,2,2); hold on
    p                           = stairs(rew_ts,rew);
    p.LineWidth                 = 2;
    p.LineStyle                 = '-';
    p.Color                     = [1 1 1];
    s.YLim                      = [0 rew(end)];
    s.YLabel.String             = ['Cumulative reward ' rew_str];
    
    x                           = [p.XData(1),repelem(p.XData(2:end),2)];       % Fill area
    y                           = [repelem(p.YData(1:end-1),2),p.YData(end)];
    fl                          = fill([x,fliplr(x)],[y,0*ones(size(y))], col{3});
    fl.FaceAlpha                = .5;
    fl.EdgeAlpha                = .5;
    fl.FaceColor                = col{3};
    fl.EdgeColor                = col{3};
    
    yyaxis right
    p2                          = plot(trg_ts,movmean(HIidx,20));
    p2.LineWidth                = 2;
    p2.Color                    = col{1};
    p3                          = plot(trg_ts,movmean(~HIidx,20));
    p3.LineWidth                = 2;
    p3.LineStyle                = ':';
    p3.Color                    = col{2};
    s.YLim                      = [0 1];
    s.XLim                      = [1 rew_ts(end)];
    s.Title.String              = 'Performance [movmean, 5 targets]';
    s.XLabel.String             = 'Time [s]';
    s.YLabel.String             = 'Rate';
    s.FontSize                  = 14;
    s.Box                       = 'off';
    s.YAxis(1).Color            = col{3};
    s.YAxis(2).Color            = col{1};
    l                           = legend([p2,p3,fl],{'HIr', 'MIr', rew_str});
    %     l                           = legend([p2,p3],{'HIr', 'MIr'});
    l.Location                  = 'west';
    
    %% ANALYSIS OF TIME WINDOW
    
    cohPool                     = unique(t.rdp_coh);                           	% Tested coherence levels
    out                         = CPR_time_window_analysis(t,29);
    summ{iSubj}                 = out;
    
    %%% PLOT STEADY STATE DATA %%%
    bx                          = [];
    by                          = [];
    bbx                         = [];
    bby                         = [];
    
    snr                         = unique(cohPool);
    arr_str                  	= cell(1,length(snr));
    arr_acc                     = cell(1,length(snr));
    
    for iCoh = 1:size(snr,1)
        arr_str{iCoh}         	= [arr_str{iCoh} out.str_mean_dist{iCoh}'];
        arr_acc{iCoh}         	= [arr_acc{iCoh} out.acc_dist{iCoh}'];
        
        by                      = [by; arr_str{iCoh}'];
        bx                      = [bx; repmat(iCoh,length(arr_str{iCoh}),1)];
        
        bby                 	= [bby; arr_acc{iCoh}'];
        bbx                  	= [bbx; repmat(iCoh,length(arr_acc{iCoh}),1)];
    end
    
    rsnr = round(snr,2);
    
    ax                          = subplot(3,3,4); hold on
    vs                        	= violinplot(bby,bbx);
    cl                        	= linspace(0,1,size(vs,2));
    for iSub = 1:size(vs,2)
        vs(iSub).ViolinColor    = [cl(iSub) 0 0];
    end
    ax.YLabel.String            = 'Tracking accuracy';
    ax.XLabel.String            = 'Coherence level';
    ax.XTickLabel               = {rsnr};
    ax.Title.String             = 'Avg motion tracking accuracy [states]';
    ax.XColor                   = [0 0 0];
    ax.YColor                   = [0 0 0];
    ax.FontSize                 = 14;
    box off
    
    % [P,ANOVATAB,STATS]          = kruskalwallis(bby,bbx,'off');
    %
    % if P > .05
    %     ax.Title.String     	= {'Avg motion tracking accuracy [states]','n.s.'};
    % end
    
    ax                          = subplot(3,3,6); hold on
    vs                        	= violinplot(by,bx);
    cl                        	= linspace(0,1,size(vs,2));
    for iSub = 1:size(vs,2)
        vs(iSub).ViolinColor    = [cl(iSub) 0 0];
    end
    ax.YLabel.String            = 'Joystick strength [norm]';
    ax.XLabel.String            = 'Coherence level';
    ax.XTickLabel               = {rsnr};
    ax.Title.String             = 'Avg joystick eccentricity [states]';
    ax.XColor                   = [0 0 0];
    ax.YColor                   = [0 0 0];
    ax.FontSize                 = 14;
    box off
    
    
    %%% PLOT TRIAL DATA %%%
    for iTrl = 1:t.trl_no(end)
        tidx = t.trl_no == iTrl;
        tt = t(tidx,:);
        
        tmp.frme_ts{iTrl}    	= [];
        tmp.rdp_dir{iTrl}    	= [];
        tmp.rdp_coh{iTrl}     	= [];
        tmp.js_dir{iTrl}      	= [];
        tmp.js_str{iTrl}      	= [];
        tmp.refresh{iTrl}      	= [];
        
        for iState = 1:size(tt,1)
            tmp.frme_ts{iTrl}  	= [tmp.frme_ts{iTrl} tt.frme_ts{iState}];
            tmp.rdp_dir{iTrl}  	= [tmp.rdp_dir{iTrl} repmat(tt.rdp_dir(iState),1,length(tt.frme_ts{iState}))];
            tmp.rdp_coh{iTrl}  	= [tmp.rdp_coh{iTrl} repmat(tt.rdp_coh(iState),1,length(tt.frme_ts{iState}))];
            tmp.js_dir{iTrl}  	= [tmp.js_dir{iTrl} tt.js_dir{iState}];
            tmp.js_str{iTrl}  	= [tmp.js_str{iTrl} tt.js_str{iState}];  
            tmp.refresh{iTrl} 	= [tmp.refresh{iTrl} median(diff(tmp.frme_ts{iTrl}))];
        end
        
        coh_block = unique(tmp.rdp_coh{iTrl});
        
      % For each coherence block...
        for iBlock = 1:size(coh_block,2)
            cohIdx           	= tmp.rdp_coh{iTrl} == coh_block(iBlock);         % Index of coherence block
            cohID(iTrl,iBlock) 	= coh_block(iBlock);                        % Coherence ID
            mStr(iTrl,iBlock)	= nanmean(tmp.js_str{iTrl}(cohIdx));                   % Average displacement for given stimulus condition
            sdStr(iTrl,iBlock)	= nanstd(tmp.js_str{iTrl}(cohIdx));                    % Standard deviation
        end
    end
    
    by                          = [];
    bx                          = [];
    bby                         = [];
    bbx                         = [];
    
    % ADD IF STATEMENT: CHECK INPUT DIMENSIONS - IS VECTOR OR MATRIX?
    for iCoh = 1:length(clvl)
        cidx                    = [];
        cidx                    = cohID == clvl(iCoh);
        
        by                      = [by; sdStr(cidx)];
        bx                      = [bx; repmat(iCoh,length(sdStr(cidx)),1)];
        
        bby                     = [bby; mStr(cidx)];
        bbx                     = [bbx; repmat(iCoh,length(mStr(cidx)),1)];
    end
    
    % ax                     	  = subplot(3,3,5); hold on
    % vs                          = violinplot(by,bx);
    % cl                       	  = linspace(0,1,size(vs,2));
    % for iSub = 1:size(vs,2)
    %     vs(iSub).ViolinColor    = [cl(iSub) 0 0];
    % end
    % ax.YLabel.String            = 'JS strength variability';
    % ax.XLabel.String            = 'Coherence level';
    % ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
    % ax.Title.String             = 'Trial-wise strength variability';
    % ax.XColor                   = [0 0 0];
    % ax.YColor                   = [0 0 0];
    % ax.FontSize                 = 14;
    % box off
    
    ax                        	= subplot(3,3,5); hold on
    vs                        	= violinplot(bby,bbx);
    cl                        	= linspace(0,1,size(vs,2));
    for iSub = 1:size(vs,2)
        vs(iSub).ViolinColor    = [cl(iSub) 0 0];
    end
    ax.YLabel.String            = 'Joystick strength [norm]';
    ax.XLabel.String            = 'Coherence level';
    ax.XTickLabel               = {rsnr};
    ax.Title.String             = 'Avg joystick eccentricity [block]';
    ax.XColor                   = [0 0 0];
    ax.YColor                   = [0 0 0];
    ax.FontSize                 = 14;
    box off
    
    %% CORRELATION ANALYSIS
    
    ax                              = subplot(3,3,7); hold on
    nLag                            = 150;
    [cr,ps]                         = CPR_correlation_analysis_WIP(tmp, nLag, true);
    crr{iSubj}                      = cr;
    cl                              = linspace(0,1,length(snr));

    ax.XTick                        = [1 nLag/2 nLag];
    ax.XLim                         = [1 nLag];
    ax.XTickLabel                   = {['-' num2str(nLag)],['-' num2str(nLag/2)],'0'};
    ax.XLabel.String                = 'Lag';
    ax.YLabel.String                = 'XCorr coeff';
    ax.FontSize                     = 14;
    ax.Title.String                 = 'XCorr(RDP, JS)';
    ax.Title.Interpreter            = 'none';
    yyaxis right
    
    for iCoh = 1:length(snr)
        vec                         = nanmean(cr.sxc(cr.coh == snr(iCoh),:));
        svec                        = smoothdata(vec,'gaussian',30);
        pm                         	= plot(svec,'Color',[cl(iCoh) 0 0 .8], 'Marker','none', 'LineStyle', '-', 'LineWidth',3);
        %         pm(iCoh)                 	= plot(svec,'Color',[cl(iCoh) 0 0 .8], 'Marker','none', 'LineStyle', '-', 'LineWidth',3);
    end
    
    ax.YAxis(2).Color               = [cl(iCoh) 0 0];
    lg                              = legend([ps, pm], {sprintf('All Trials\n[smoothed]'),sprintf('Mean\n[smoothed]')});
    lg.Position(1)                  = .13;
    
    %     lg                              = legend([ps, pm(1),pm(2),pm(3),pm(4),pm(5)], ...
    %         {'Trials',['Mean ' num2str(round(snr(1),2))],['Mean ' num2str(round(snr(2),2))],...
    %         ['Mean ' num2str(round(snr(3),2))],['Mean ' num2str(round(snr(4),2))],...
    %         ['Mean ' num2str(round(snr(5),2))]},'Location', 'northwest');
    
    %%% Coherence %%%
    bx                              = [];
    by                              = [];
    bbx                             = [];
    bby                             = [];
    
    snr                             = unique(cr.coh);
    ccorr                           = cell(1,length(snr));
    aucPeak                         = cell(1,length(snr));
    xcLag                           = cell(1,length(snr));
    
    for iCoh = 1:size(snr,2)
        cIdx                        = cr.coh == snr(iCoh);
        ccorr{iCoh}                 = [ccorr{iCoh} cr.cc(cIdx)'];
        aucPeak{iCoh}               = [aucPeak{iCoh} cr.auPk(cIdx)];
        xcLag{iCoh}                 = [xcLag{iCoh} cr.posPk(cIdx)'-nLag];
        
        by                          = [by; ccorr{iCoh}];
        bx                          = [bx; repmat(iCoh,length(ccorr{iCoh}),1)];
        
        %     by                          = [by; aucPeak{iCoh}];
        %     bx                          = [bx; repmat(iCoh,length(aucPeak{iCoh}),1)];
        
        bby                         = [bby; xcLag{iCoh}];
        bbx                         = [bbx; repmat(iCoh,length(xcLag{iCoh}),1)];
    end
    
    ax                              = subplot(3,3,9); hold on
    vs                              = violinplot(by,bx);
    cl                              = linspace(0,1,size(vs,2));
    for iSub = 1:size(vs,2)
        vs(iSub).ViolinColor        = [cl(iSub) 0 0];
    end
    ax.YLabel.String                = 'Area under peak';
    ax.YLabel.String                = 'Corr coeff';
    ax.XLabel.String                = 'Coherence level';
    ax.XTickLabel                   = {rsnr};
    ax.Title.String                 = 'Circular correlation';
    ax.XColor                       = [0 0 0];
    ax.YColor                       = [0 0 0];
    ax.FontSize                     = 14;
    box off
    
    ax                              = subplot(3,3,8); hold on
    vs                              = violinplot(bby,bbx);
    cl                              = linspace(0,1,size(vs,2));
    for iSub = 1:size(vs,2)
        vs(iSub).ViolinColor        = [cl(iSub) 0 0];
    end
    ax.YLabel.String                = 'Lag';
    ax.XLabel.String                = 'Coherence level';
    ax.XTickLabel                   = {rsnr};
    ax.Title.String                 = 'XCorr peak lag';
    ax.XColor                       = [0 0 0];
    ax.YColor                       = [0 0 0];
    ax.FontSize                     = 14;
    box off
    
    if strcmp(sbj{iSubj}, 'cla') || strcmp(sbj{iSubj}, 'nil')
        %         dest_dir = '/Users/fschneider/Documents/MWorks/Plots/';
        %         print(f, [dest_dir 'summary_' sbj{iSubj} '_' fname], '-r300', '-dpng');
        dest_dir = '/Users/fschneider/Desktop/';
        print(f, [dest_dir 'summary_' sbj{iSubj} '_' fname], '-r300', '-dpng');
    else
        if strcmp(sbj{iSubj}, 'agnt')
            dest_dir = ['/Volumes/T7_Shield/CPR_psychophysics/' sbj{1} '/summary/'];
        else
            dest_dir = ['/Volumes/T7_Shield/CPR_psychophysics/' sbj{iSubj} '/summary/'];
        end
        
        if iscell(fname)
            fid = strsplit(fname{iSubj},'_');
        else
            fid = strsplit(fname,'_');
        end
        print(f, [dest_dir fid{1}  '_' sbj{iSubj} '_' fid{3} '_' fid{4}], '-r300', '-dpng');
    end
end
end