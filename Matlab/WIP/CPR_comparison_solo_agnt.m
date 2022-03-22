addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/Documents/MATLAB/cbrewer

close all
clear all

fname_solo = {'20220215_ana',...
    '20220214_sol',...
    '20220126_dih',...
    '20220201_mam',...
    '20220128_nak',...
    '20220216_tos',...
    '20220223_cah',...
    '20220222_anm'}; % SOLO DATA         



fname_agent = {'20220127_ana',...
    '20220204_sol',...  
    '20220131_dih',...
    '20220228_mam',...
    '20220224_nak',...
    '20220201_tos',...
    '20220301_cah',...
    '20220301_anm'}; % AGNT DATA       


pth                 = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/';
nSample             = 30;

for iCond = 1:2
    
    if iCond == 1
        fname = fname_solo;
    elseif iCond == 2
        fname = fname_agent;
    end
    
    for iSub = 1:size(fname,2)
        disp(['Condition_' num2str(iCond) '_Subject_' num2str(iSub)])
        clear t
        load([pth fname{iSub}(end-2:end) '/summary/' fname{iSub} '_cpr_tbl.mat'])
        snr                     = unique(t.ss_coh);
        id{iSub}                = fname{iSub}(end-2:end);
        
        % Correlation analysis - for coherence blocks
        clear tmp tIdx tbl_trl cr ps
        
        tmp                    	= cellfun(@(x) size(x), t.trl_rdp_dir, 'UniformOutput', false);
        tIdx                   	= cellfun(@(x) x(2)>1, tmp, 'UniformOutput', false);
        tIdx                  	= cell2mat(tIdx);
        tbl_trl                	= t(tIdx,:);
        nLag                 	= 150;
        [cr,ps]               	= CPR_correlation_analysis_WIP(tbl_trl, nLag, true);
        
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
            hir(iSub,iCoh,iCond)    = nhi / ntrg;
            
            % Target score
            trg_score(iSub,iCoh,iCond) 	= nanmean(score(score_coh  == snr(iCoh) & score_hi == true));
            
            % Joystick displacement
            mstr(iSub,iCoh,iCond) 	= nanmedian(cellfun(@(x) nanmean(x(end-nSample:end)), t.js_str(cIdx)));
            
            % Joystick accuracy
            rdp_dir                 = t.rdp_dir(cIdx);
            js_dir                  = t.js_dir(cIdx);
            
            for iState = 1:length(rdp_dir)
                clear js_dev
                js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iState}(end-nSample:end)),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
                js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));                           % Joystick accuracy
            end
            
            macc(iSub, iCoh,iCond)	= nanmedian(js_acc);
            sxc{iSub,iCoh,iCond} 	= cr.sxc(cr.coh == snr(iCoh),:);
            mR(iSub,iCoh,iCond)    	= mean(cr.maxR(cr.coh == snr(iCoh)));
            mcc(iSub,iCoh,iCond)   	= mean(cr.cc(cr.coh == snr(iCoh)));
            mpk(iSub,iCoh,iCond)  	= mean(cr.posPk(cr.coh == snr(iCoh)));
            spk(iSub,iCoh,iCond)  	= std(cr.posPk(cr.coh == snr(iCoh)));
        end
    end
end

%% PLOT DIFF

lb = {'Hit Rate',...
    'Hit Score',...
    'Avg JS Displacement',...
    'Avg JS Accuracy',...
    'Avg CirCorr Coeff',...
    'Avg XCorr Peak Position',...
    'Std XCorr Peaks',...
    'Avg XCorr Coeff'};

f                       = figure('units','normalized','outerposition',[0 0 1 1]);
hold on

for iDat = 1:8
    
    if iDat == 1
        d = hir;
    elseif iDat == 2
        d = trg_score;
    elseif iDat == 3
        d = mstr;
    elseif iDat == 4
        d = macc;
    elseif iDat == 5
        d = mcc;
    elseif iDat == 6
        d = mpk;
    elseif iDat == 7
        d = spk;
    elseif iDat == 8
        d = mR;
    end
    
    xtick                   = [1 9 length(snr)];
    ax                      = subplot(2,4,iDat);
    im                      = imagesc(d(:,:,1)-d(:,:,2));
    ax.Title.String         = lb{iDat};
    ax.Title.Interpreter    = 'none';
    ax.FontSize             = 16;
    ax.YTick                = 1:size(fname,2);
    ax.YTickLabel           = cellfun(@(x) x(end-2:end),fname,'UniformOutput',false);
    ax.XTick                = xtick;
    ax.XTickLabel           = snr(xtick);
    ax.XLabel.String        = 'Coherence';
    cb                      = colorbar;
    
    if iDat == 1
        caxis([-.25 .25])
    elseif iDat == 2
        caxis([-.301 .301])
    elseif iDat == 3
        caxis([-.35 .35])
    elseif iDat == 4
        caxis([-.2 .2])
    elseif iDat == 5
        caxis([-.5 .5])
    elseif iDat == 6
        caxis([-60 60])
    elseif iDat == 7
        caxis([-50 50])
    elseif iDat == 8
        caxis([-1000 1000])
    end
    
    colormap(cbrewer('div', 'RdGy', 256, 'PCHIP'))  
    axis square;

end


print(f, ['/Users/fschneider/Desktop/solo_minus_agnt'], '-r300', '-dpng');

%% PLOT LAG

for j = 1:size(fname,2)
    solo                = [];
    agnt                = [];
    
    for i = 1:length(snr)
        solo            = [solo; mean(sxc{j,i,1})];
        agnt            = [agnt; mean(sxc{j,i,2})];
    end
    
    ms(j,:)             = mean(solo);
    ma(j,:)             = mean(agnt);
    
    p(j,1)              = find(mean(solo) == max(mean(solo)));
    p(j,2)              = find(mean(agnt) == max(mean(agnt)));
end

c                       = cbrewer('qual', 'Dark2', size(fname,2), 'PCHIP');
alph                    = .75;

f                       = figure('units','normalized','outerposition',[0 0 1 1]);
ax                      = subplot(2,3,1);
im                      = imagesc(solo);
ax.YTick                = 1:length(snr);
ax.YTickLabel           = snr;
ax.YLabel.String        = 'Coherence';
ax.XTick                = [1 151 301];
ax.XTickLabel           = round([-150 0 150] * 8.333);
ax.XLabel.String        = 'Lag [ms]';
ax.Title.String         = ['Exmpl Subj: ' fname_solo{j}(end-2:end) '_solo'];
ax.FontSize             = 20;
ax.Title.Interpreter    = 'none';
cb                      = colorbar;
cb.Label.String         = 'Xcorr coeff';

caxis([0 750])
colormap(gray)

ax                      = subplot(2,3,2);
hold on
for i = 1:size(ms,1)
pl                      = plot(ms(i,:),'LineWidth',2,'Color',[c(i,:) alph]);
end
ax.XTick                = [1 151 301];
ax.XTickLabel           = round([-150 0 150] * 8.333);
ax.XLabel.String        = 'Lag [ms]';
ax.YLabel.String        = 'XCorr Coeff';
ax.FontSize         	= 20;
ax.Title.String         = 'Average/Subj: Solo';
ax.Position(1)          = ax.Position(1)+.02;

ax                      = subplot(2,3,4);
im                      = imagesc(agnt);
ax.YTick                = 1:length(snr);
ax.YTickLabel           = snr;
ax.YLabel.String        = 'Coherence';
ax.XTick                = [1 151 301];
ax.XTickLabel           = round([-150 0 150] * 8.333);
ax.XLabel.String        = 'Lag [ms]';
ax.Title.String         = ['Exmpl Subj: ' fname_agent{j}(end-2:end) '_computer'];
ax.Title.Interpreter    = 'none';
ax.FontSize             = 20;
cb                      = colorbar;
cb.Label.String         = 'Xcorr coeff';

caxis([0 750])
colormap(gray)

ax                      = subplot(2,3,5);
hold on
for i = 1:size(ma,1)
pl                      = plot(ma(i,:),'LineWidth',2,'Color',[c(i,:) alph]);
end
ax.XTick                = [1 151 301];
ax.XTickLabel           = round([-150 0 150] * 8.333);
ax.XLabel.String        = 'Lag [ms]';
ax.YLabel.String        = 'XCorr Coeff';
ax.FontSize             = 20;
ax.Title.String       	= 'Average/Subj: Computer';
ax.Position(1)          = ax.Position(1)+.02;

dat                     = ((p(:,1)-p(:,2))*8.3333);
ax                  	= subplot(2,3,[3 6]);
bx                      = boxplot(dat, 'Colors', 'k');

set(bx(end,:),'Visible','off')
set(bx, {'linew'},{3})
hold on

for i = 1:length(dat)
    sc               	= scatter(1,((p(i,1)-p(i,2))*8.3333));
    sc.SizeData         = 100;
    sc.MarkerFaceColor  = c(i,:);
    sc.MarkerFaceAlpha  = 1;
    sc.MarkerEdgeColor  = 'none';
end

ax.YLabel.String        = 'Peak difference [ms]';
ax.XTick                = [];
ax.Box                  = 'off';
ax.FontSize             = 20;

lg                      = legend(cellfun(@(x) x(end-2:end),fname,'UniformOutput',false));
lg.Location             = 'east';
lg.Box                  = 'off';

print(f, ['/Users/fschneider/Desktop/Pk_diff'], '-r300', '-dpng');
