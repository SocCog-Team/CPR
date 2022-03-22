addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab

clear all 

tbl1 = load('/Users/fschneider/Desktop/pilot/dyadic/20211217_nes_cpr_tbl');
tbl2 = load('/Users/fschneider/Desktop/pilot/dyadic/20211217_ann_cpr_tbl');
clear t

tbl1 = tbl1.t;
tbl2 = tbl2.t;
nLag = 150;

% Check number of samples
tmp = cellfun(@size, tbl2.trl_frme_ts, 'UniformOutput', false);
idx = cell2mat(cellfun(@(x) x(2)>1, tmp, 'UniformOutput', false));

for iSubj = 1:2
    count = 0;
    
    if iSubj == 1
        t = tbl1;
    else
        t = tbl2;
    end
    
    trl_lst = find(idx);
    for i = 1:length(trl_lst)
        clear js_dir rdp_dir js_dev js_dff rdp_dff rdp_coh
        iRow                    = trl_lst(i);
        js_dir                  = t.trl_js_dir{iRow};
        rdp_dir                 = t.trl_rdp_dir{iRow};
        rdp_coh                 = t.trl_rdp_coh{iRow};
        ts                      = t.trl_frme_ts{iRow};
        js_dev                  = rad2deg(circ_dist(deg2rad(js_dir),deg2rad(rdp_dir)));  % Minimum RDP-Joystick difference
        
        for iSmple = 1:size(t.trl_frme_ts{iRow},2)-1
            js_dff(iSmple)     	= rad2deg(circ_dist(deg2rad(js_dir(iSmple)),deg2rad(js_dir(iSmple+1))));
            rdp_dff(iSmple)     = rad2deg(circ_dist(deg2rad(rdp_dir(iSmple)),deg2rad(rdp_dir(iSmple+1))));
        end
        
        ex                      = isnan(js_dff) | isnan(rdp_dff);
        js_dff(ex)              = [];
        rdp_dff(ex)             = [];
        jdff{i,iSubj}           = js_dff;
        
        ex                      = isnan(js_dir) | isnan(js_dev) | isnan(rdp_dir);
        js_dir(ex)              = [];
        js_dev(ex)              = [];
        rdp_dir(ex)             = [];
        rdp_coh(ex)             = [];
        ts(ex)                  = [];
        
        % Extract coherence chunks
        cindx                	= diff(rdp_coh) ~=0;
        cid                   	= rdp_coh([true cindx]);
        cts                     = ts([true cindx]);
        
        for iCoh = 1:size(cts,2)
            
            % Build coherence index
            if iCoh < size(cts,2)
                cIdx                = ts >= cts(iCoh) & ts < cts(iCoh+1);
            else
                cIdx                = ts >= cts(iCoh) & ts <= ts(end);
            end
            
            if sum(cIdx) == 0
                continue
            end
            
            % Correlation analysis
            count                   = count + 1;
            coh(count,iSubj)        = cid(iCoh);
            [xc(count,iSubj,:),lags]= xcorr(abs(rdp_dff(cIdx(1:end-1))),abs(js_dff(cIdx(1:end-1))),nLag,'coeff');
            sxc(count,iSubj,:)   	= smoothdata(xc(count,iSubj,:) ,'gaussian',20);         % Smooth correlation output with gaussian kernel
            maxR(count,iSubj)      	= max(sxc(count,iSubj,1:nLag+1));                        % Max cross-correlation coefficient
            posPk(count,iSubj)    	= find(sxc(count,iSubj,1:nLag+1) == maxR(count,iSubj));
            cc(count,iSubj)        	= circ_corrcc(deg2rad(js_dir(cIdx)),deg2rad(rdp_dir(cIdx)));	% Circular correlation between stimulus and joystick direction
            acc(count,iSubj)       	= abs(1 - median(abs(js_dev(cIdx))) / 180);
        end
    end
end

snr = unique(coh);
for iCoh = 1:size(snr,1)
    ind             = coh(:,1) == snr(iCoh);
    macc{iCoh}      = acc(ind,:);
    msxc1(iCoh,:)   = mean(squeeze(sxc(ind,1,:)));
    msxc2(iCoh,:)   = mean(squeeze(sxc(ind,2,:)));
    mpk{iCoh}       = posPk(ind,:);
end

%% Accuracy

figure
dat1 =[];
dat2 =[];
cat = [];
for i = 1:length(macc)
    dat1 = [dat1; macc{i}(:,1)];
    dat2 = [dat2; macc{i}(:,2)];
    cat = [cat; repmat(i,[1,length(macc{i}(:,1))])'];
end

v1 = violinplot(dat1,cat);
v2 = violinplot(dat2,cat);
for i = 1:length(v1)
    v1(i).ViolinColor = [1 0 0];
    v2(i).ViolinColor = [0 0 1];
    v1(i).ViolinAlpha = .1;
    v2(i).ViolinAlpha = .1;
end

title('Response accuracy')
xlim([0 length(macc)+1])
xlabel('Coherence')
ylabel('Accuracy [norm]')
set(gca,'FontSize',16)
set(gca,'XTickLabels',snr)

%% XC Peak distribution

figure
dat1 =[];
dat2 =[];
cat = [];
for i = 1:length(mpk)
    dat1 = [dat1; mpk{i}(:,1)];
    dat2 = [dat2; mpk{i}(:,2)];
    cat = [cat; repmat(i,[1,length(mpk{i}(:,1))])'];
end

v1 = violinplot((dat1-nLag)*1e3*1/120,cat);
v2 = violinplot((dat2-nLag)*1e3*1/120,cat);
for i = 1:length(v1)
    v1(i).ViolinColor = [1 0 0];
    v2(i).ViolinColor = [0 0 1];
    v1(i).ViolinAlpha = .1;
    v2(i).ViolinAlpha = .1;
end

title('xCorr peak distribution')
xlim([0 length(mpk)+1])
xlabel('Coherence')
ylabel('xCorr peak [ms]')
set(gca,'FontSize',16)
set(gca,'XTickLabels',snr)


%% XC curves

figure;hold on
alp = linspace(0,1,17);
for iCoh = 1:size(msxc1,1)
    plot(lags*1e3*1/120,msxc1(iCoh,:),'Color',[1 0 0 alp(iCoh)], 'LineWidth',1)
    plot(lags*1e3*1/120,msxc2(iCoh,:),'Color',[0 0 1 alp(iCoh)], 'LineWidth',1)
end

title('Crosscorrelation RDP-JS')
xlim([-1250 0])
xlabel('Lag [ms]')
ylabel('xCorr coeff')
set(gca,'FontSize',20)

%% XC between subjects
clearvars -except tbl1 tbl2 nLag trl_lst

count = 0;

for i = 1:length(trl_lst)
    iRow                    = trl_lst(i);
    js1_dir                 = tbl1.trl_js_dir{iRow};
    js2_dir                 = tbl2.trl_js_dir{iRow};
    rdp_coh                 = tbl1.trl_rdp_coh{iRow};
    ts                      = tbl1.trl_frme_ts{iRow};
    js_dev                  = rad2deg(circ_dist(deg2rad(js1_dir),deg2rad(js2_dir)));  % Minimum Joystick-Joystick difference
    
    for iSmple = 1:size(ts,2)-1
        js1_dff(iSmple)    	= rad2deg(circ_dist(deg2rad(js1_dir(iSmple)),deg2rad(js1_dir(iSmple+1))));
        js2_dff(iSmple)    	= rad2deg(circ_dist(deg2rad(js2_dir(iSmple)),deg2rad(js2_dir(iSmple+1))));
    end
    
    ex                      = isnan(js1_dff) | isnan(js2_dff);
    js1_dff(ex)             = [];
    js2_dff(ex)             = [];
    
    ex                      = isnan(js1_dir) | isnan(js_dev) | isnan(js2_dir) | isnan(rdp_coh);
    js1_dir(ex)             = [];
    js2_dir(ex)             = [];
    js_dev(ex)              = [];
    rdp_coh(ex)             = [];
    ts(ex)                  = [];
    
    % Extract coherence chunks
    cindx                	= diff(rdp_coh) ~=0;
    cid                   	= rdp_coh([true cindx]);
    cts                     = ts([true cindx]);
    
    for iCoh = 1:size(cts,2)
        
        % Build coherence index
        if iCoh < size(cts,2)
            cIdx                = ts >= cts(iCoh) & ts < cts(iCoh+1);
        else
            cIdx                = ts >= cts(iCoh) & ts <= ts(end);
        end
        
        if sum(cIdx) == 0
            continue
        end
        
        % Correlation analysis
        count                   = count + 1;
        coh(count)              = cid(iCoh);
        [xc(count,:), lags]    	= xcorr(abs(js1_dff(cIdx(1:end-1))),abs(js2_dff(cIdx(1:end-1))),nLag,'coeff');
        sxc(count,:)            = smoothdata(xc(count,:) ,'gaussian',20);         % Smooth correlation output with gaussian kernel
        maxR(count)             = max(sxc(count,:));                        % Max cross-correlation coefficient
        posPk(count)            = find(sxc(count,:) == maxR(count));
        cc(count)               = circ_corrcc(deg2rad(js1_dir(cIdx)),deg2rad(js2_dir(cIdx)));	% Circular correlation between stimulus and joystick direction
        acc(count)              = abs(1 - median(abs(js_dev(cIdx))) / 180);
    end
end

snr = unique(coh);
for iCoh = 1:size(snr,2)
    ind             = coh == snr(iCoh)';
    macc{iCoh}      = acc(ind);
    msxc(iCoh,:)    = mean(squeeze(sxc(ind,:)));
    mpk{iCoh}       = posPk(ind);
end

%% Accuracy

figure
dat =[];
cat = [];
for i = 1:length(macc)
    dat = [dat; macc{i}'];
    cat = [cat; repmat(i,[1,length(macc{i})])'];
end

v = violinplot(dat,cat);
for i = 1:length(v)
    v(i).ViolinColor = [.5 .5 .5];
    v(i).ViolinAlpha = .3;
end

title('Similarity between joysticks')
xlim([0 length(mpk)+1])
xlabel('Coherence')
ylabel('Joystick aligment')
set(gca,'FontSize',16)
set(gca,'XTickLabels',snr)

%% XC Peak distribution

figure; hold on
dat =[];
cat = [];
for i = 1:length(mpk)
    dat = [dat; mpk{i}'];
    cat = [cat; repmat(i,[1,length(mpk{i})])'];
end

line([0 length(mpk)+1],[0 0], 'Color','k')
v = violinplot(dat-nLag,cat);
for i = 1:length(v)
    v(i).ViolinColor = [.5 .5 .5];
    v(i).ViolinAlpha = .3;
end

title('xCorr [Js-Js] peak distribution')
xlim([0 length(mpk)+1])
xlabel('Coherence')
ylabel('xCorr peak position')
set(gca,'FontSize',16)
set(gca,'XTickLabels',snr)

%% XC curves

figure;hold on
alp = linspace(0,1,17);
for iCoh = 1:size(msxc,1)
    plot(lags*1e3*1/120,msxc(iCoh,:),'Color',[alp(iCoh) 0 0], 'LineWidth',1)
%     plot(msxc(iCoh,:)./max(msxc(iCoh,:),[],2),'Color',[alp(iCoh) 0 0], 'LineWidth',1)
end

title('xCorr between JS1 and JS2')
set(gca,'FontSize', 16)
xlim([-1250 1250])
xlabel('Lag [ms]')
ylabel('xCorr coeff')

%% Raw XC

figure
imagesc(msxc./max(msxc,[],2))
ylabel('Coherence')
xlabel('Lag')
title('xCorr between JS1 and JS2 [norm]')
set(gca,'FontSize', 16)
set(gca,'YTick',[1:17])
set(gca,'YTickLabels',snr)
set(gca,'XTick',[1 151 301])
set(gca,'XTickLabels',{[lags(1) lags(151) lags(301)]})
colorbar

%%
% Determine similarity between joystick signals:
% Accuracy
% Displacement
% Circular correlation
% Assessment of psychometric function

% Determine if Leader-Follower present:
% Crosscorrelation -> Lag distribution
% Granger causality
% Shifted, Kalman-filtered signal to predict timecourse of other joystick

