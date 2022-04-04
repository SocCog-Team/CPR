addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/Documents/MATLAB/cbrewer

close all
clear all

fname = {'20220215_ana',...
    '20220222_anm',...
    '20220223_cah',...
    '20220215_dem',...
    '20220126_dih',...
    '20220201_mam',...
    '20220128_nak',...
    '20220214_sol',...
    '20220225_anb',...
    '20220325_ant',...
    '20220312_seo',...
    '20220322_las'};%    '20220216_tos',...


pth                 = '/Users/fschneider/Documents/CPR_psychophysics/';
nSample             = 30;

for iSub = 1:size(fname,2)
    clear tmp1 tmp2 t1 t2 t
    tmp1                    = load([pth fname{iSub}(end-2:end) '/summary/' fname{iSub} '_CPRsolo_block1_tbl.mat']);
    tmp2                    = load([pth fname{iSub}(end-2:end) '/summary/' fname{iSub} '_CPRsolo_block2_tbl.mat']);
    t                       = [tmp1.t; tmp2.t];
    
    snr                     = unique(t.ss_coh);
    id{iSub}                = fname{iSub}(end-2:end);
    
    % Plot target reward distribution
    ff                      = CPR_reward_distribution(t);
    print(ff, [pth id{iSub} '/summary/' 'summary_' id{iSub} '_' fname{iSub}(1:8) '_CPRsolo_rewMatrix'], '-r300', '-dpng');

    % Correlation analysis - for coherence blocks
    clear tmp tIdx tbl_trl cr ps
    
    tmp                    	= cellfun(@(x) size(x), t.trl_rdp_dir, 'UniformOutput', false);
    tIdx                   	= cellfun(@(x) x(2)>1, tmp, 'UniformOutput', false);
    tIdx                  	= cell2mat(tIdx);
    tbl_trl                	= t(tIdx,:);
    nLag                 	= 150;
    [cr,ps]               	= CPR_correlation_analysis_WIP(tbl_trl, nLag, false);
    
    % Target score
    clear score score_hi score_coh
    c                    	= 0;
    for i = 1:size(t.trg_ts,1)
        for iTrg = 1:length(t.trg_ts{i})
            
            [~,idx]         = min(abs(t.frme_ts{i} - t.trg_ts{i}(iTrg)));
            str             = t.js_str{i}(idx);
            acc             = abs(1 - abs(t.rdp_dir(i) - t.js_dir{i}(idx)) / 180);
            
            c = c +1;
            score(c)        = str * acc;
            score_coh(c)    = t.ss_coh(i);
            score_hi(c)     = t.trg_hit{i}(iTrg);            
        end
    end

    for iCoh = 1:length(snr)
        clear cIdx rdp_dir t.js_dir
        
        cIdx = t.ss_coh == snr(iCoh);
        cIdx(cellfun(@length,t.js_str) < 100) = false;
        
        % Hit rate
        nhi                     = sum(cellfun(@sum,t.trg_hit(cIdx)));
        ntrg                    = sum(cellfun(@numel,t.trg_hit(cIdx)));
        hir(iSub,iCoh)          = nhi / ntrg;
        
        % Target score
        trg_score(iSub,iCoh) 	= nanmean(score(score_coh  == snr(iCoh) & score_hi == true));
        
        % Joystick displacement
        mstr(iSub, iCoh)        = nanmedian(cellfun(@(x) nanmedian(x(end-nSample:end)), t.js_str(cIdx)));
        
        
        % Joystick accuracy
        rdp_dir                 = t.rdp_dir(cIdx);
        js_dir                  = t.js_dir(cIdx);
        
        for iState = 1:length(rdp_dir)
            clear js_dev
            js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iState}(end-nSample:end)),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
            js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));                           % Joystick accuracy
        end
        
        macc(iSub, iCoh)        = nanmedian(js_acc);
        sxc{iSub,iCoh}      	= cr.sxc(cr.coh == snr(iCoh),:);
        mR(iSub,iCoh)         	= mean(cr.maxR(cr.coh == snr(iCoh)));
        mcc(iSub,iCoh)         	= mean(cr.cc(cr.coh == snr(iCoh)));
        mpk(iSub,iCoh)       	= mean(cr.posPk(cr.coh == snr(iCoh)));
        spk(iSub,iCoh)        	= std(cr.posPk(cr.coh == snr(iCoh)));
    end
end


%% Plot

lb = {'Hit rate',...
    'Avg accuracy',...
    'Avg circular correlation coeff',...
    'Avg peak position',...
    'Avg peak variability',...
    'Avg joystick displacement [raw]',...
    'Avg joystick displacement [norm]',...
    'Avg target score [norm]'};

str = {'hir','macc','mcc','mpk','spk','mstr','mstr_norm','trg_score'};

for iPlot = 1:length(str)
    
    if iPlot == 1
        var     = hir;
    elseif iPlot == 2
        var     = macc;
    elseif iPlot == 3
        var     = mcc;
    elseif iPlot == 4
        var     = mpk-nLag;
    elseif iPlot == 5
        var     = spk;
    elseif iPlot == 6
        var     = mstr;
    elseif iPlot == 7
        var     = mstr ./ mean(mstr,2);
    elseif iPlot == 8
        var     = trg_score ./ mean(trg_score,2);
    end
    
    f                   = figure;
    imsc                = imagesc(var);
    ax                  = gca;
    ax.YTick            = [1:size(id,2)];
    ax.YTickLabel       = id;
    ax.XLabel.String    = 'Coherence';
    ax.XTick            = 1:length(snr);
    ax.XTickLabel       = snr;
    ax.FontSize         = 16;
    cb                  = colorbar;
    cb.Label.String     = lb{iPlot};
    
    if iPlot == 2
        caxis([.75 1])
    end
    
    axis image
    cmap=cbrewer('seq', 'OrRd', 256, 'PCHIP');
    colormap(flipud(cmap))
    print(f, ['/Users/fschneider/Desktop/' str{iPlot}], '-r300', '-dpng');
end

%% Crosscorrelation

alph = linspace(.1,.7,size(sxc,2));
cmap=cbrewer('qual', 'Dark2', size(sxc,1), 'PCHIP');

for iSub = 1:size(sxc,1)
    f = figure;
    hold on
    clear tmp
    tmp = cellfun(@mean,sxc(iSub,:),'UniformOutput',false);
    for iCoh = 1:size(sxc,2)
        pl = plot(tmp{iCoh});
        pl.LineWidth = 2;
        pl.Color = [cmap(iSub,:) alph(iCoh)];
    end
    
    ax                  = gca;
    ax.YLabel.String    = 'XCorr Coeff';
    ax.XLabel.String    = 'Lag [ms]';
    ax.XLim             = [0 301];
    ax.XTickLabel       = round((cellfun(@str2num, ax.XTickLabel)-150) * 8.33);
    ax.FontSize         = 16;
    ax.Title.String     = id{iSub};
    print(f, ['/Users/fschneider/Desktop/Sub' num2str(iSub)], '-r300', '-dpng');
end