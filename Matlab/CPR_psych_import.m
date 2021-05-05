function [t,d, out] = CPR_psych_import(subj, pth, fname, d)

% Add relevant directories
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

cd(pth)

%% Import data

% Define variables for import
var_import = {
    'ML_', ...
    'CTRL_', ...
    'RDP_', ...
    'INFO_', ...
    'TRIAL_', ...
    'IO_joystickDirection', ...
    'IO_joystickStrength' };

if nargin < 4                                                             
    d                	= CPR_import_mwk2(fname, var_import, false);
end

dte                     = regexp(fname,'\d{8}');
if strcmp(subj, 'cla') || strcmp(subj, 'nil')
    fid              	= fname(dte(2):dte(2)+7);
else
    fid               	= fname(dte:dte+7);
end

%% Extract data and construct table

% Initialise variables
idx                         = [];
trl                         = [];

% Create variable-specific indices
idx.tOn                     = d.event == 'ML_trialStart';
idx.tEnd                    = d.event == 'ML_trialEnd';
idx.steady_duration         = d.event == 'CTRL_SteadyStateDuration_ms';
idx.RDP_dir                 = d.event == 'RDP_direction';
idx.RDP_coh                 = d.event == 'RDP_coherence';
idx.JS_dir                  = d.event == 'IO_joystickDirection';
idx.JS_str                  = d.event == 'IO_joystickStrength';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.trg                     = d.event == 'TRIAL_reactionTrigger';
idx.stateOn                 = d.event == 'INFO_SteadyStateCounter';

if strcmp(subj, 'cla') || strcmp(subj, 'nil')
    idx.reward              = d.event == 'INFO_Juice_ml';
    rew_str                 = 'ml';
else
    if str2num(fid) < 20210401
        idx.reward      	= d.event == 'INFO_ExtraCash';
        rew_str         	= 'EUR';
    else
        idx.reward      	= d.event == 'INFO_Cash';
        rew_str          	= 'EUR';
    end
end

if sum(idx.tOn) == 0
    idx.tOn                 = d.event == 'TRIAL_start';
    idx.tEnd                = d.event == 'TRIAL_end';
end

% Trial timestamps
trl.tOn                     = d.time(idx.tOn);
trl.tEnd                    = d.time(idx.tEnd);

% % Remove trials with duration <30us
% excl                        = trl.tEnd - trl.tOn < 30;
% trl.tOn(excl)               = [];
% trl.tEnd(excl)              = [];

t                           = CPR_construct_table(subj, fid, d, trl, idx, true);
            
%% PLOT PERFORMANCE

f                           = figure('Units', 'normalized', 'Position', [0 0 .8 1]); set(gcf,'color', [1 1 1]);

HIidx                       = t.trg_hit(t.trg_shown);
HIr                         = sum(HIidx) / length(HIidx);
MIr                         = sum(~HIidx) / length(HIidx);
trg_ts                      = t.trg_ts(~isnan(t.trg_ts));
trg_ts                      = (trg_ts - trl.tOn(1)) ./ 1e6;
rew                         = d.value(idx.reward);
exIdx                       = strcmp(cellfun(@class, rew, 'UniformOutput', false),'double');
rew                         = cell2mat(rew(exIdx));
rew_ts                      = d.time(idx.reward);
rew_ts                      = rew_ts(exIdx);
rew_ts                      = (rew_ts - trl.tOn(1)) ./ 1e6;

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
tcoh                        = t.ss_coh(t.trg_shown);
clvl                        = unique(tcoh);
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
s.YLabel.String             = ['Cumulative reward' rew_str];

x                           = [p.XData(1),repelem(p.XData(2:end),2)];       % Fill area
y                           = [repelem(p.YData(1:end-1),2),p.YData(end)];
fl                          = fill([x,fliplr(x)],[y,0*ones(size(y))], col{3});
fl.FaceAlpha                = .5;
fl.EdgeAlpha                = .5;
fl.FaceColor                = col{3};
fl.EdgeColor                = col{3};

yyaxis right
p2                          = plot(trg_ts,movmean(HIidx,5));
p2.LineWidth                = 2;
p2.Color                    = col{1};
p3                          = plot(trg_ts,movmean(~HIidx,5));
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
l                           = legend([p2,p3,fl],{'HIr', 'MIr', 'EUR'});
l.Location                  = 'west';

%% ANALYSIS OF TIME WINDOW

out                         = CPR_time_window(t,29);
cohPool                     = unique(t.ss_coh);                           	% Tested coherence levels
cohPool(cohPool == 0)       = [];                                         	% Zero not tested, throw out

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

[P,ANOVATAB,STATS]          = kruskalwallis(bby,bbx,'off');

if P > .05
    ax.Title.String     	= {'Avg motion tracking accuracy [states]','n.s.'};
end

ax                          = subplot(3,3,6); hold on
vs                        	= violinplot(by,bx);
cl                        	= linspace(0,1,size(vs,2));
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor    = [cl(iSub) 0 0];
end
ax.YLabel.String            = 'Joystick strength [norm]';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr};
ax.Title.String             = 'Avg radial joystick displacement [states]';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
box off


%%% PLOT TRIAL DATA %%%
tmp = cellfun(@size, t.trl_rdp_dir, 'UniformOutput', false);
for i = 1:size(tmp,1)
    indx(i,:) = tmp{i}(2) > 1;
end

% Get trial data
trl_str                     = t.trl_js_str(indx);                           % Joystick strength
trl_str_ts                  = t.trl_js_str_ts(indx);
tmp_coh                     = t.trl_rdp_coh(indx);                          % RDP coherence
tmp_coh_ts                  = t.trl_rdp_coh_ts(indx);
tmp_coh{1}(1)               = [];                                           % Remove first entry of first trial
tmp_coh_ts{1}(1)            = [];

for iTrl = 1:size(trl_str,1)
    clear trl_data trl_data_ts trl_coh trl_coh_ts
    trl_data = trl_str{iTrl};
    trl_data_ts = trl_str_ts{iTrl};
    trl_coh = tmp_coh{iTrl};
    trl_coh_ts = tmp_coh_ts{iTrl};
    
    for iCoh = 1:size(trl_coh,2)
        % Coherence index
        if iCoh < size(trl_coh,2)
            cohIdx       	= trl_data_ts >= trl_coh_ts(iCoh) & trl_data_ts < trl_coh_ts(iCoh+1);
        else
            cohIdx          = trl_data_ts >= trl_coh_ts(iCoh) & trl_data_ts <= trl_data_ts(end);
        end
        
        cohID(iTrl,iCoh)  	= trl_coh(iCoh);                                % Coherence ID
        mStr(iTrl,iCoh) 	= mean(trl_data(cohIdx));                       % Average strength for given 
        sdStr(iTrl,iCoh)  	= std(trl_data(cohIdx));                        % Standard deviation
    end
end

cohID(cohID == 0)           = nan;
by                          = [];
bx                          = [];
bby                         = [];
bbx                         = [];

for iCoh = 1:length(clvl)
    cidx                    = [];
    cidx                    = cohID == clvl(iCoh);
    
    by                      = [by; sdStr(cidx)];
    bx                      = [bx; repmat(iCoh,length(sdStr(cidx)),1)];
   
    bby                     = [bby; mStr(cidx)];
    bbx                     = [bbx; repmat(iCoh,length(mStr(cidx)),1)];
end

% ax                        	= subplot(3,3,5); hold on
% vs                        	= violinplot(by,bx);
% cl                        	= linspace(0,1,size(vs,2));
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
ax.Title.String             = 'Avg radial joystick displacement [trial]';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
box off

%% CORRELATION ANALYSIS

ax                              = subplot(3,3,7); hold on
nLag                            = 150;
[cr,ps]                         = CPR_correlation(t,indx, nLag, true);

cohID(cohID == 0)               = nan;              
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
end

ax.YAxis(2).Color               = [cl(iCoh) 0 0];
lg                              = legend([ps, pm], {sprintf('All Trials\n[smoothed]'),sprintf('Mean\n[smoothed]')});
lg.Position(1)                  = .13;

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

if strcmp(subj, 'cla') || strcmp(subj, 'nil')
    dest_dir = '/Users/fschneider/Documents/MWorks/XBI/Plots/';
    print(f, [dest_dir 'summary_' subj '_' fid], '-r300', '-dpng');
else
    print(f, [pth '/summary_' subj '_' fid], '-r300', '-dpng');
end
end
