%% Concatenate subject tables
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a

clear all 
close all

pth     = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/Pilot_free_viewing/';
cd(pth)

% fnames  = {
%     'fxs_cpr_20210204_mac.mwk2';
%     'fxs_cpr_20210204_bew.mwk2';
%     'fxs_cpr_20210204_kan.mwk2';
%     'fxs_cpr_20210204_fxs.mwk2';
%     'fxs_cpr_20210205_pas.mwk2';
%     'fxs_cpr_20210205_sem.mwk2';
%     'fxs_cpr_20210205_stm.mwk2';
%     'fxs_cpr_20210205_lac.mwk2';
%     'fxs_cpr_20210210_ilv.mwk2';
%     'fxs_cpr_20210210_nes.mwk2';
%     'fxs_cpr_20210224_piy.mwk2';
%     'fxs_cpr_20210318_sut.mwk2';
%     };
% 
% tpop = [];
% 
% for iSubj = 1:size(fnames,1)
%     clear t
%     load([fnames{iSubj} '_tbl.mat'])    
%     tpop = [tpop; t];
%     
%     load([fnames{iSubj} '.mat'])
%     [~,~,out{iSubj}] = CPR_psych_import(fnames{iSubj}(end-7:end-5),pth,fnames{iSubj},d);
%     close all
% end
% 
% clear t
% t = tpop;
% save('pop_tbl.mat', 't', '-v7.3')                              % Save as .mat file
% save('pop_summary.mat', 'out', '-v7.3')          

%% PERFORMANCE

load('/Users/fschneider/ownCloud/CPR_data/Pilot_free_viewing/CPR_pop_tbl.mat')
load('/Users/fschneider/ownCloud/CPR_data/Pilot_free_viewing/CPR_pop_summary.mat')

f                   = figure('Units', 'normalized', 'Position', [0 0 .65 1]); set(gcf,'color', [1 1 1]);

HIidx               = t.trg_hit(logical(t.trg_shown));
HIr                 = sum(HIidx) / length(HIidx);
MIr                 = sum(~HIidx) / length(HIidx);

%%% Performance pie chart
s                   = subplot(3,3,1);
p                   = pie([HIr,MIr], [1,1]);
pPatch              = findobj(p,'Type','patch');
pText               = findobj(p,'Type','text');
percentValues       = get(pText,'String');
txt                 = {'HIr: ';'MIr: '};
combinedtxt         = strcat(txt,percentValues);
s.Title.String      = {'Overall', 'Performance'};
s.Title.FontSize    = 14;
s.Title.Position(1) = -1;
% s.Position(1)       = 0;
col                 = {[1 0 0],[0 0 0]};

for i = 1:size(pText,1)
    pText(i).String         = combinedtxt(i);
    pText(i).FontSize       = 14;
    pPatch(i).FaceColor     = col{i};
end

% Condition-wise performance
tcoh                = t.ss_coh(logical(t.trg_shown));
clvl                = unique(tcoh);
for iCoh = 1:length(clvl)
    cindx           = tcoh == clvl(iCoh);
    hir(iCoh)    	= sum(HIidx(cindx)) / length(HIidx(cindx));
    mir(iCoh)     	= sum(~HIidx(cindx)) / length(HIidx(cindx));
end

s                   = subplot(3,3,2);
bp                	= bar([hir',mir'], 'stacked');
cl                  = [1 0];

for iOutc = 1:2
    bp(iOutc).FaceColor    	= [cl(iOutc) 0 0];
    bp(iOutc).EdgeColor     = [cl(iOutc) 0 0];
end

s.XTick            = [1:length(clvl)];
s.XLim             = [0.5 5.5];
s.YLabel.String   	= 'Rate';
s.XLabel.String 	= 'Coherence level';
s.XTickLabel    	= {round(clvl,2)};
s.FontSize      	= 12;
axis tight

%% ANALYSIS OF TIME WINDOW

nSamples                = 29;                                                   	% Number of samples before direction changes
cohPool                 = unique(t.ss_coh);                                        % Tested coherence levels
cohPool(cohPool == 0)   = [];

for iCoh = 1:size(cohPool,1)                                                        % For each coherence level
    str_arr           	= [];                                                       % Reset temporary variables
    dir_arr             = [];
    trg_shown           = [];
    trg_hit             = [];
    
    cohIdx              = t.ss_coh == cohPool(iCoh);                               % Coherence index
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
    
    sm.str_mean(iCoh)      = mean(mean(str_arr,2));                               	% Average across steady states
    sm.str_std(iCoh)       = std(std(str_arr,[],2));                                % Standard deviation across steady states
    sm.HIr(iCoh)           = sum(trg_hit) / sum(trg_shown);                         % Hit rate
    sm.acc(iCoh)           = mean(mean(js_acc{iCoh},2));                            % Response accuracy
end

%%% STEADY STATE DATA %%%
bx                          = [];
by                          = [];
bbx                         = [];
bby                         = [];

snr                         = unique(cohPool);
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

% sct = scatterhist(by,bby, ...
%     'Kernel', 'on', ...
%     'Direction', 'out',...
%     'Location','SouthEast',...
%     'Color','k',...
%     'Marker','.',...
%     'MarkerSize',6);
% xlim([0 1.25])
% ylabel('Accuracy')
% xlabel('Confidence')
% set(gca, 'fontsize',16);

rsnr = round(snr,2);

ax                          = subplot(3,3,3); hold on
vs                        	= violinplot(bby,bbx);
cl                        	= linspace(0,1,size(vs,2));
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor    = [cl(iSub) 0 0];
end
ax.YLabel.String            = 'Tracking accuracy';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = {'Avg motion tracking accuracy', 'Final 300ms/steady state'};
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
box off
axis tight


ax                          = subplot(3,3,6); hold on
vs                        	= violinplot(by,bx);
cl                        	= linspace(0,1,size(vs,2));
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor    = [cl(iSub) 0 0];
end
ax.YLabel.String            = 'Joystick strength [norm]';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = {'Avg radial joystick displacement', 'Final 300ms/steady state'};
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
box off
axis tight


%%% TRIAL DATA %%%
tmp = cellfun(@size, t.trl_rdp_dir, 'UniformOutput', false);
for i = 1:size(tmp,1)
    indx(i,:) = tmp{i}(2) > 1;
end
coh                         = t.ss_coh(indx);
trl_str                     = t.trl_js_str(indx);

for iTrl = 1:size(trl_str,1)
    mStr(iTrl)              = mean(trl_str{iTrl});
    sdStr(iTrl)             = std(trl_str{iTrl});
end

by                          = [];
bx                          = [];
bby                         = [];
bbx                         = [];

for iCoh = 1:length(clvl)
    cidx                    = [];
    cidx                    = coh == clvl(iCoh);
    
    by                      = [by; sdStr(cidx)'];
    bx                      = [bx; repmat(iCoh,length(sdStr(cidx)),1)];
   
    bby                     = [bby; mStr(cidx)'];
    bbx                     = [bbx; repmat(iCoh,length(mStr(cidx)),1)];
end

ax                        	= subplot(3,3,5); hold on
vs                        	= violinplot(by,bx);
cl                        	= linspace(0,1,size(vs,2));
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor    = [cl(iSub) 0 0];
end
ax.YLabel.String            = 'JS strength variability';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr};
ax.Title.String             = {'Radial joystick displacement','Trial-wise variability'};
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
box off
axis tight

ax                        	= subplot(3,3,4); hold on
vs                        	= violinplot(bby,bbx);
cl                        	= linspace(0,1,size(vs,2));
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor    = [cl(iSub) 0 0];
end
ax.YLabel.String            = 'Avg JS strength';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr};
ax.Title.String             = {'Radial joystick displacement','Trial-wise average'}';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
box off
axis tight

%% XCORR ANALYSIS

ax                  = subplot(3,3,7); hold on
nLag                = 150;                                                      	% XCorr Lag

rdp_dir             = t.trl_rdp_dir(indx);
rdp_dir_ts          = t.trl_rdp_dir_ts(indx);
js_dir              = t.trl_js_dir(indx);
js_dir_ts           = t.trl_js_dir_ts(indx);
js_str              = t.trl_js_str(indx);
js_str_ts           = t.trl_js_str_ts(indx);

for i = 1:size(rdp_dir_ts,1)
    
    % Adjust vector length
    vec         = [];
    for ii = 1:size(rdp_dir_ts{i},2)-1
        ssIdx   = [];
        ssIdx   = js_dir_ts{i} >= rdp_dir_ts{i}(ii) & js_dir_ts{i} < rdp_dir_ts{i}(ii+1);
        vec     = [vec repmat(rdp_dir{i}(ii),1,sum(ssIdx))];
    end
    
    ssIdx       = js_dir_ts{i} >= rdp_dir_ts{i}(end);
    vec         = [vec repmat(rdp_dir{i}(end),1,sum(ssIdx))];
    
    if size(vec,2) ~= size(js_dir{i},2)
        vec     = [vec nan(1,size(js_dir{i},2) - size(vec,2))];
    end
    
    % Correct for circular space
    clear js_dff js_corr rdp_dff rdp_corr
    
    js_dff      = mod(diff(js_dir{i}) + 180,360) - 180;
    js_corr     = js_dff(1);
    rdp_dff     = mod((diff(vec)+180),360)-180;
    rdp_corr    = rdp_dff(1);
    
    for iSample = 1:length(js_dir{i})-1
        js_corr(iSample+1)  = js_dff(iSample) + js_corr(iSample);
        rdp_corr(iSample+1) = rdp_dff(iSample) + rdp_corr(iSample);
    end
    
    xc(i,:)               	= xcorr(double(abs(rdp_dff)),abs(js_dff),nLag);               	 % Crosscorrelation between stimulus and joystick direction
    maxR(i)                 = max(xc(i,1:nLag+1));
    posPk(i,:)              = find(xc(i,1:nLag+1) == max(xc(i,1:nLag+1)));
    
    try
        auPk(i,:)           = trapz(xc(i,posPk(i,:)-10:posPk(i,:)+10));
    catch
        auPk(i,:)         	= nan;
    end
    
    ps                      = plot(xc(i,1:nLag),'Color',[.5 .5 .5 .1], 'Marker','none', 'LineWidth',1.5);    % Plot trial-wise cross-correlation
    
end

ax.XTick                    = [1 nLag/2 nLag];
ax.XLim                     = [1 nLag];
ax.XTickLabel               = {['-' num2str(nLag*10)],['-' num2str((nLag/2)*10)],'0'};
ax.XLabel.String            = 'Lag [ms]';
ax.YLabel.String            = 'XCorr coeff';
ax.FontSize                 = 14;
ax.Title.String             = 'XCorr(RDP, JS)';
ax.Title.Interpreter        = 'none';
yyaxis right
pm                          = plot(mean(xc),'Color',[1 0 0], 'Marker','none', 'LineWidth',3);
ax.YAxis(2).Color           = [1 0 0];
lg                          = legend([ps, pm], {'Trials','Mean'});
lg.Position(1)              = 0.03;
axis tight

%%% Coherence %%%
bx                          = [];
by                          = [];
bbx                         = [];
bby                         = [];

snr                         = unique(coh);
aucPeak                     = cell(1,length(snr));
xcLag                       = cell(1,length(snr));

for iCoh = 1:size(snr,1)
    cIdx                 	= coh == snr(iCoh);
    aucPeak{iCoh}           = [aucPeak{iCoh} auPk(cIdx)'];
    xcLag{iCoh}             = [xcLag{iCoh} posPk(cIdx)' - nLag];
    
    by                      = [by; aucPeak{iCoh}'];
    bx                      = [bx; repmat(iCoh,length(aucPeak{iCoh}),1)];
    
    bby                     = [bby; xcLag{iCoh}'];
    bbx                     = [bbx; repmat(iCoh,length(xcLag{iCoh}),1)];
end

ax                        	= subplot(3,3,8); hold on
vs                        	= violinplot(by,bx);
cl                        	= linspace(0,1,size(vs,2));
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor    = [cl(iSub) 0 0];
end
ax.YLabel.String            = 'Area under peak';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = 'Area under XCorr peak';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
box off
axis tight

ax                        	= subplot(3,3,9); hold on
vs                        	= violinplot(bby*10,bbx);
cl                        	= linspace(0,1,size(vs,2));
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor    = [cl(iSub) 0 0];
end
ax.YLabel.String            = 'Lag [ms]';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = 'XCorr peak lag';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
box off
axis tight

dest_dir = '/Users/fschneider/ownCloud/CPR_data/Pilot_free_viewing/';
print(f, [dest_dir 'summary_pop'], '-r300', '-dpng');

%% Subject wise

clear lg pk

for iSub = 1:length(out)
    
    avg_str(iSub,:) = out{iSub}.str_mean;
    avg_acc(iSub,:) = out{iSub}.acc;
    
    clear cpool cval cidx
    cpool           = out{iSub}.coh;
    cval            = unique(cpool);
    
    for iCoh = 1:size(cval,1)
        cidx = cpool == cval(iCoh);
        pk(iSub,iCoh) = nanmedian(out{iSub}.auPk(cidx));
        lg(iSub,iCoh) = nanstd(out{iSub}.posPk(cidx));
        mlg(iSub,iCoh) = nanmedian(out{iSub}.posPk(cidx));
    end
end
    
f                           = figure('Units', 'normalized', 'Position', [0 0 1 1]); set(gcf,'color', [1 1 1]);

ax = subplot(2,3,1);
imagesc(avg_acc./mean(avg_acc,2))
% imagesc(avg_acc)
ax.YLabel.String            = 'Subject';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = 'Avg motion tracking accuracy';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
colorbar
cb                          = colorbar;
cb.Label.String             = 'Tracking accuracy [norm]';
cb.Label.FontSize           = 14;

ax = subplot(2,3,2);
imagesc(avg_str./mean(avg_str,2))
% imagesc(avg_str)
ax.YLabel.String            = 'Subject';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = 'Avg joystick displacement';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
cb                          = colorbar;
cb.Label.String             = 'Joystick displacement [norm]';
cb.Label.FontSize           = 14;

% ax = subplot(2,3,3); 
% imagesc(avg_str./mean(avg_str,2))
% ax.YLabel.String            = 'Subject';
% ax.XLabel.String            = 'Coherence level';
% ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
% ax.Title.String             = 'Avg Strength [norm]';
% ax.XColor                   = [0 0 0];
% ax.YColor                   = [0 0 0];
% ax.FontSize                 = 14;
% cb                          = colorbar;
% cb.Label.String             = 'Avg strength [norm]';

ax = subplot(2,3,4);
imagesc(pk./mean(pk,2))
ax.YLabel.String            = 'Subject';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = 'Avg area under XC peak';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
cb                          = colorbar;
cb.Label.String             = 'Area under xcorr peak [norm]';
cb.Label.FontSize           = 14;

ax = subplot(2,3,5);
imagesc(mlg./mean(mlg,2))
% imagesc(mlg)
ax.YLabel.String            = 'Subject';
ax.XLabel.String            = 'Coherence ';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = 'Avg XC peak lag';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
cb                          = colorbar;
cb.Label.String             = 'Xcorr peak lag [norm]';
cb.Label.FontSize           = 14;

ax = subplot(2,3,6);
imagesc(lg./mean(lg,2))
% imagesc(lg)
ax.YLabel.String            = 'Subject';
ax.XLabel.String            = 'Coherence ';
ax.XTickLabel               = {rsnr(1),rsnr(2),rsnr(3),rsnr(4),rsnr(5)};
ax.Title.String             = 'Avg XC peak lag variability';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
cb                          = colorbar;
cb.Label.String             = 'Xcorr peak lag variability [norm]';
cb.Label.FontSize           = 14;

cm = [linspace(0,1,256)' zeros(256,2)];
colormap(cm)

dest_dir = '/Users/fschneider/ownCloud/CPR_data/Pilot_free_viewing/';
print(f, [dest_dir 'summary_subj'], '-r300', '-dpng');