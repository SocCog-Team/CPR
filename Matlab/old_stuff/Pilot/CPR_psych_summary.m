clear all 
close all

addpath /Users/fschneider/Documents/GitHub/CPR/Matlab
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a

cd('/Users/fschneider/ownCloud/Shared/SFB_1528_A1_Project - @SFB1528 [ownCloud, FS, ST]/Summary_plots')
load('pop_summary.mat')
load('pop_tbl.mat')

%% ANALYSIS OF TIME WINDOW

nSamples                = 29;                                                   	% Number of samples before direction changes
cohPool                 = unique(t.trl_coh);                                        % Tested coherence levels
cohPool(cohPool == 0)   = [];

for iCoh = 1:size(cohPool,1)                                                        % For each coherence level
    str_arr           	= [];                                                       % Reset temporary variables
    dir_arr             = [];
    trg_shown           = [];
    trg_hit             = [];
    
    cohIdx              = t.trl_coh == cohPool(iCoh);                               % Coherence index
    str_raw{iCoh}       = t.js_str(cohIdx);                                         % Joystick strength
    dir_raw{iCoh}       = t.js_dir(cohIdx);                                         % Joystick direction
    rdp_dir{iCoh}       = t.rdp_dir(cohIdx);                                        % Stimulus direction
    
    for iSS = 1:size(str_raw{iCoh},1)                                               % For each steady state...
        str_arr(iSS,:)  	= str_raw{iCoh}{iSS}(end-nSamples:end);                 % Take last nSamples+1 data points [Fs = 100Hz -> 10ms] of each steady state
        dir_arr(iSS,:)   	= dir_raw{iCoh}{iSS}(end-nSamples:end);
    end
      
    js_dev{iCoh}            = rad2deg(circ_dist(deg2rad(dir_arr),deg2rad(t.rdp_dir(cohIdx))));  % Raw RDP-Joystick difference
    js_acc{iCoh}            = abs(1 - abs(js_dev{iCoh}) / 180);                                 % Joystick accuracy
    trg_hit                	= t.trg_hit(cohIdx);
    trg_shown            	= t.trg_shown(cohIdx);
    
    sm.str_mean_dist{iCoh}	= mean(str_arr,2);
    sm.str_std_dist{iCoh} 	= std(str_arr,[],2);
    sm.acc_dist{iCoh}    	= mean(js_acc{iCoh},2);
    
    sm.str_mean(iCoh)       = mean(mean(str_arr,2));                               	 % Average across steady states
    sm.str_std(iCoh)        = std(std(str_arr,[],2));                                % Standard deviation across steady states
    sm.HIr(iCoh)            = sum(trg_hit) / sum(trg_shown);                         % Hit rate
    sm.acc(iCoh)            = mean(mean(js_acc{iCoh},2));                            % Response accuracy
end


%%% SUBJECT DATA %%%
clear mlg slg pk perf
for iSub = 1:length(out)
    
    avg_str(iSub,:)         = out{iSub}.str_mean;
    avg_acc(iSub,:)         = out{iSub}.acc;
    perf(iSub,:)            = (out{iSub}.HIr);

    clear cpool cval cidx
    cpool                   = out{iSub}.coh;
    cval                    = unique(cpool);
    
    for iCoh = 1:size(cval,1)
        cidx                = cpool == cval(iCoh);
        pk(iSub,iCoh)       = nanmedian(out{iSub}.auPk(cidx));
        slg(iSub,iCoh)      = nanstd(out{iSub}.posPk(cidx));
        mlg(iSub,iCoh)      = nanmedian(out{iSub}.posPk(cidx));
    end
end

%%% POOLED STEADY STATE DATA %%%
bx                          = [];
by                          = [];
bbx                         = [];
bby                         = [];

snr                         = unique(cohPool);
rsnr                        = round(snr,2);
arr_str                  	= cell(1,length(snr));
arr_acc                     = cell(1,length(snr));

for iCoh = 1:size(snr,1)
    arr_str{iCoh}         	= [arr_str{iCoh} sm.str_mean_dist{iCoh}'];
    arr_acc{iCoh}         	= [arr_acc{iCoh} sm.acc_dist{iCoh}'];
    
    by                      = [by; arr_str{iCoh}'];
    bx                      = [bx; repmat(iCoh,length(arr_str{iCoh}),1)];
    
    bby                 	= [bby; arr_acc{iCoh}'];
    bbx                  	= [bbx; repmat(iCoh,length(arr_acc{iCoh}),1)];
end

%%% PLOT %%%
alph                        = .15;
yoff                        = .12;
cfs                         = 12;
fs                          = 14;
cl                        	= linspace(1,.4,size(arr_str,2));
% cl                        	= ones(1,size(arr_str,2));
yofs                        = .03;
clm                         = [.1 .6];
row                         = [.76 .47 .12];
dim                         = [0.3 0.2];

f                          	= figure('Units', 'normalized', 'Position', [0 0 .5 1]);
% ax                          = subplot(3,2,3); hold on
ax                       	= axes('Position',[clm(1) row(2) dim]); hold on
im                          = imagesc(avg_acc./mean(avg_acc,2));
ax.YLabel.String            = 'Subject';
cb                          = colorbar;
cb.Location                 = 'southoutside';
cb.Position(3:4)            = [.1 .01];
cb.Position(1)              = ax.Position(1) + ax.Position(3)/3;
cb.Position(2)              = cb.Position(2) - yoff;
cb.Title.String             = 'Avg tracking accuracy [norm]';
cb.FontSize                 = cfs;
cb.Ticks                    = [.8 1.2];
caxis([0.8 1.2])
colormap(gray)

yyaxis right
vs                        	= violinplot(bby,bbx);
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor            = [cl(iSub) 0 0];
    vs(iSub).ViolinPlot.LineStyle   = 'none';
    vs(iSub).ViolinPlot.Marker   	= 'none';
    vs(iSub).ViolinAlpha            = alph;
    vs(iSub).BoxPlot.Marker       	= 'none';
    vs(iSub).BoxPlot.FaceColor      = [1 1 1];
    vs(iSub).WhiskerPlot.LineStyle  = 'none';
    vs(iSub).WhiskerPlot.Marker   	= 'none';
    vs(iSub).ScatterPlot.MarkerFaceAlpha = alph;
end
ax.YLabel.String            = 'Tracking accuracy';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = {'Motion tracking accuracy'};
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = fs;
ax.YAxis(2).Color           = [1 0 0];
box off
axis tight

ax                       	= axes('Position',[clm(1) row(3) dim]); hold on
im                          = imagesc(avg_str./mean(avg_str,2));
ax.YLabel.String            = 'Subject';
cb                          = colorbar;
cb.Location                 = 'southoutside';
cb.Position(3:4)            = [.1 .01];
cb.Position(1)              = ax.Position(1) + ax.Position(3)/3;
cb.Position(2)              = cb.Position(2) - yoff;
cb.Title.String             = 'Avg joystick displacement [norm]';
cb.FontSize                 = cfs;
cb.Ticks                    = [.6 1.4];
caxis([0.6 1.4])
colormap(gray)

yyaxis right
vs                        	= violinplot(by./max(by),bx);
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor            = [cl(iSub) 0 0];
    vs(iSub).ViolinPlot.LineStyle   = 'none';
    vs(iSub).ViolinPlot.Marker   	= 'none';
    vs(iSub).ViolinAlpha            = alph;
    vs(iSub).BoxPlot.Marker       	= 'none';
    vs(iSub).BoxPlot.FaceColor      = [1 1 1];
    vs(iSub).WhiskerPlot.LineStyle  = 'none';
    vs(iSub).WhiskerPlot.Marker   	= 'none';
    vs(iSub).ScatterPlot.MarkerFaceAlpha = alph;
end
ax.YLabel.String            = 'Joystick displacement [norm]';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = {'Radial joystick displacement'};
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = fs;
ax.YAxis(2).Color           = [1 0 0];
box off
axis tight


%%% Condition-wise performance %%%
HIidx               = t.trg_hit(logical(t.trg_shown));
HIr                 = sum(HIidx) / length(HIidx);
MIr                 = sum(~HIidx) / length(HIidx);
tcoh                = t.trl_coh(logical(t.trg_shown));
clvl                = unique(tcoh);
for iCoh = 1:length(clvl)
    cindx           = tcoh == clvl(iCoh);
    hir(iCoh)    	= sum(HIidx(cindx)) / length(HIidx(cindx));
    mir(iCoh)     	= sum(~HIidx(cindx)) / length(HIidx(cindx));
end

% s                   = subplot(3,2,1);
ax                  = axes('Position',[clm(1) row(1) dim]); hold on
bp                	= bar([hir',mir'], 'stacked');
col                 = [1 0];

for iOutc = 1:2
    bp(iOutc).FaceColor    	= [col(iOutc) 0 0];
    bp(iOutc).EdgeColor     = [col(iOutc) 0 0];
end

ax.XTick         	= [1:length(clvl)];
ax.XLim          	= [0.5 5.5];
ax.Title.String  	= 'Performance';
ax.YLabel.String   	= 'Rate';
ax.XLabel.String 	= 'Coherence level';
ax.XTickLabel    	= {round(clvl,2)};
ax.FontSize      	= fs;
axis tight
box off
lg               	= legend({'Hit','Miss'});
lg.Location         = 'northwest';


ax                  = axes('Position',[clm(2) row(1) dim]); hold on
nLag                = 150;                                                      	% XCorr Lag
yofs                = .015;
tmp                 = cellfun(@size, t.trl_rdp_dir, 'UniformOutput', false);
for i = 1:size(tmp,1)
    indx(i,:) = tmp{i}(2) > 1;
end

rdp_dir             = t.trl_rdp_dir(indx);
rdp_dir_ts          = t.trl_rdp_dir_ts(indx);
js_dir              = t.trl_js_dir(indx);
js_dir_ts           = t.trl_js_dir_ts(indx);
js_str              = t.trl_js_str(indx);
js_str_ts           = t.trl_js_str_ts(indx);
sub                 = t.ID(indx);

for i = 1:size(rdp_dir_ts,1)
    
    % Adjust vector length
    vec             = [];
    for ii = 1:size(rdp_dir_ts{i},2)-1
        ssIdx       = [];
        ssIdx       = js_dir_ts{i} >= rdp_dir_ts{i}(ii) & js_dir_ts{i} < rdp_dir_ts{i}(ii+1);
        vec         = [vec repmat(rdp_dir{i}(ii),1,sum(ssIdx))];
    end
    
    ssIdx           = js_dir_ts{i} >= rdp_dir_ts{i}(end);
    vec             = [vec repmat(rdp_dir{i}(end),1,sum(ssIdx))];
    
    if size(vec,2) ~= size(js_dir{i},2)
        vec         = [vec nan(1,size(js_dir{i},2) - size(vec,2))];
    end
    
    % Correct for circular space
    clear js_dff js_corr rdp_dff rdp_corr
    
    js_dff          = mod(diff(js_dir{i}) + 180,360) - 180;
    js_corr         = js_dff(1);
    rdp_dff         = mod((diff(vec)+180),360)-180;
    rdp_corr        = rdp_dff(1);
    
    for iSample = 1:length(js_dir{i})-1
        js_corr(iSample+1)  = js_dff(iSample) + js_corr(iSample);
        rdp_corr(iSample+1) = rdp_dff(iSample) + rdp_corr(iSample);
    end
    
    xc(i,:)               	= xcorr(double(abs(rdp_dff)),abs(js_dff),nLag);               	 % Crosscorrelation between stimulus and joystick direction
    maxR(i)                 = max(xc(i,1:nLag+1));
    posPk(i,:)              = find(xc(i,1:nLag+1) == max(xc(i,1:nLag+1)));
    cc(i,:)               	= circ_corrcc(deg2rad(js_dir{i}),deg2rad(vec));               	 % Circular correlation between stimulus and joystick direction

    try
        auPk(i,:)           = trapz(xc(i,posPk(i,:)-10:posPk(i,:)+10));
    catch
        auPk(i,:)         	= nan;
    end
    
    ps                      = plot(xc(i,1:nLag),'Color',[.5 .5 .5 alph], 'Marker','none', 'LineWidth',2);    % Plot trial-wise cross-correlation
    
end

ax.XTick                    = [1 nLag/2 nLag];
ax.XLim                     = [1 nLag];
ax.XTickLabel               = {['-' num2str(nLag*10)],['-' num2str((nLag/2)*10)],'0'};
ax.XLabel.String            = 'Lag [ms]';
ax.YLabel.String            = 'XCorr coeff';
ax.FontSize                 = fs;
ax.Title.String             = {'Cross-correlation'};
ax.Title.Interpreter        = 'none';
yyaxis right
pm                          = plot(mean(xc),'Color',[1 0 0], 'Marker','none', 'LineWidth',3);
ax.YAxis(2).Color           = [1 0 0];
lg                          = legend([ps, pm], {'Trial','Mean'});
lg.Location                 = 'northwest';
axis tight

%%% Coherence %%%
bx                          = [];
by                          = [];
bbx                         = [];
bby                         = [];

coh                         = t.trl_coh(indx);
snr                         = unique(coh);
% aucPeak                     = cell(1,length(snr));
ccorr                       = cell(1,length(snr));
xcLag                       = cell(1,length(snr));
sbjct                       = cell(1,length(snr));

for iCoh = 1:size(snr,1)
    cIdx                 	= coh == snr(iCoh);
%     aucPeak{iCoh}           = [aucPeak{iCoh} auPk(cIdx)'];
    ccorr{iCoh}             = [ccorr{iCoh} cc(cIdx)'];
    xcLag{iCoh}             = [xcLag{iCoh} posPk(cIdx)' - nLag];
    sbjct{iCoh}             = [sbjct{iCoh} sub(cIdx)];
    
%     by                      = [by; aucPeak{iCoh}'];
%     bx                      = [bx; repmat(iCoh,length(aucPeak{iCoh}),1)];

    by                      = [by; ccorr{iCoh}'];
    bx                      = [bx; repmat(iCoh,length(ccorr{iCoh}),1)];
    
    bby                     = [bby; xcLag{iCoh}'];
    bbx                     = [bbx; repmat(iCoh,length(xcLag{iCoh}),1)];
end

usub                        = unique(sub);
for iSub = 1:length(usub)
    for iCoh = 1:length(snr)
        sIdx                = strcmp(sbjct{iCoh},usub(iSub));
        mcc(iSub,iCoh)      = nanmean(ccorr{iCoh}(sIdx));
    end
end

alph = .4;
ax                       	= axes('Position',[clm(2) row(3) dim]); hold on
% im                          = imagesc(pk./mean(pk,2));
im                          = imagesc(mcc./mean(mcc,2));
ax.YLabel.String            = 'Subject';
cb                          = colorbar;
cb.Location                 = 'southoutside';
cb.Position(3:4)            = [.1 .01];
cb.Position(1)              = ax.Position(1) + ax.Position(3)/3;
cb.Position(2)              = cb.Position(2) - yoff;
cb.Title.String             = 'Corr coeff [norm]';
cb.FontSize                 = cfs;
cb.Ticks                    = [-.5 3];
caxis([-.5 3])
colormap(gray)

yyaxis right
vs                        	= violinplot(by,bx);
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor            = [cl(iSub) 0 0];
    vs(iSub).ViolinPlot.LineStyle   = 'none';
    vs(iSub).ViolinPlot.Marker   	= 'none';
    vs(iSub).ViolinAlpha            = alph;
    vs(iSub).BoxPlot.Marker       	= 'none';
    vs(iSub).BoxPlot.FaceColor      = [1 1 1];
    vs(iSub).WhiskerPlot.LineStyle  = 'none';
    vs(iSub).WhiskerPlot.Marker   	= 'none';
    vs(iSub).ScatterPlot.MarkerFaceAlpha = alph;
end
% ax.YLabel.String            = 'Area under peak';
ax.YLabel.String            = 'Correlation coeff';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
% ax.Title.String             = 'Area under XCorr peak';
ax.Title.String             = 'Circular correlation';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = fs;
ax.YAxis(2).Color           = [1 0 0];
box off
axis tight

ax                       	= axes('Position',[clm(2) row(2) dim]); hold on
im                          = imagesc(slg./mean(slg,2));
ax.YLabel.String            = 'Subject';
cb                          = colorbar;
cb.Location                 = 'southoutside';
cb.Position(3:4)            = [.1 .01];
cb.Position(1)              = ax.Position(1) + ax.Position(3)/3;
cb.Position(2)              = cb.Position(2) - yoff;
cb.Title.String             = 'Avg peak variability [norm]';
cb.FontSize                 = cfs;
cb.Ticks                    = [0 1.5];
caxis([0 1.5])
colormap(gray)

yyaxis right
vs                        	= violinplot(bby*10,bbx);
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor            = [cl(iSub) 0 0];
    vs(iSub).ViolinPlot.LineStyle   = 'none';
    vs(iSub).ViolinPlot.Marker   	= 'none';
    vs(iSub).ViolinAlpha            = alph;
    vs(iSub).BoxPlot.Marker       	= 'none';
    vs(iSub).BoxPlot.FaceColor      = [1 1 1];
    vs(iSub).WhiskerPlot.LineStyle  = 'none';
    vs(iSub).WhiskerPlot.Marker   	= 'none';
    vs(iSub).ScatterPlot.MarkerFaceAlpha = alph;
end
ax.YLabel.String            = 'Lag [ms]';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = 'XCorr peak lag variability';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = fs;
ax.YAxis(2).Color           = [1 0 0];
box off
axis tight

ax0 = axes('Position',[0 0 1 1],'Visible','off');
text(.03,.98, 'A', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.03,.7, 'B', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.03,.35, 'C', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.51,.98, 'D', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.51,.7, 'E', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')
text(.51,.35, 'F', 'Parent', ax0, 'FontSize', 30, 'Color', 'k', 'FontWeight', 'bold')

dest_dir = '/Users/fschneider/Desktop/';
print(f, [dest_dir 'summary_psych'], '-r300', '-dpng');