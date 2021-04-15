function [t,d] = XBI_cpr_train_import(monk, fname, d)

% Add relevant directories
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a

cd /Users/fschneider/Documents/MWorks/XBI
% cd /Users/fschneider/Documents/MWorks/Data

%% Check input arguments & import data

% Define variables for import
var_import = {
    'ML_', ...
    'CTRL_', ...
    'RDP_', ...
    'INFO_', ...
    'TRIAL_', ...
    'IO_joystickDirection', ...
    'IO_joystickStrength' };

if nargin < 3                                                               % If no data structure provided...
    if iscell(fname)                                                        % If multiple files ...
        d               = mergeFiles(fname, var_import);                    % ... merge to single structure
    else
        if isfile([fname '.mat'])                                           % If .mat file available...
            tmp         = load([fname '.mat']);                            	% ...load .mat file
            d           = tmp.d;
        else
            d           = MW_readFile(fname, 'include', var_import);       	% Import .mwk2 sesion file
          
            disp('Save struct...')
            save([fname '.mat'], 'd', '-v7.3')                              % Save as .mat file
            disp('Done!')
        end
    end
end

if iscell(fname)                                                            % Get session ID
    fid                	= fname{1}(31:38);
else
    fid               	= fname(31:38);
end

%% Initialise variables

% Counter
cc                      = 0;
idx                     = [];
trl                     = [];

% Table specs
var     = {'ID', 'string'; ...
    'trl_no', 'double'; ...
    'trl_dur', 'double'; ...
    'trl_coh', 'double'; ...
    'ss_no', 'double'; ...
    'ss_dur', 'double'; ...
    'rdp_dir', 'double';...
    'trg_shown', 'logical'; ...
    'trg_hit', 'logical'; ...
    'trg_ts', 'double'; ...
    'js_dir', 'cell'; ...
    'js_dir_ts', 'cell'; ...
    'js_str', 'cell'; ...
    'js_str_ts', 'cell'; ...
    'trl_rdp_dir', 'cell'; ...
    'trl_rdp_dir_ts', 'cell'; ...
    'trl_js_dir', 'cell'; ...
    'trl_js_dir_ts', 'cell'; ...
    'trl_js_str', 'cell'; ...
    'trl_js_str_ts', 'cell'};

% Initialise table
t       = table('Size',[5000, size(var,1)],...
    'VariableTypes',var(:,2),...
    'VariableNames',var(:,1));

%% Extract data

% Create variable-specific indices
idx.tOn                 = d.event == 'TRIAL_start';
idx.tEnd                = d.event == 'TRIAL_end';
idx.steady_duration     = d.event == 'CTRL_SteadyStateDuration_ms';
idx.RDP_dir             = d.event == 'RDP_direction';
idx.RDP_coh             = d.event == 'RDP_coherence';
idx.JS_dir              = d.event == 'IO_joystickDirection';
idx.JS_str              = d.event == 'IO_joystickStrength';
idx.outcome             = d.event == 'TRIAL_outcome';
idx.trg                 = d.event == 'TRIAL_reactionTrigger';
idx.reward              = d.event == 'INFO_Juice_ml';

% Trial timestamps
trl.tOn                 = d.time(idx.tOn);
trl.tEnd                = d.time(idx.tEnd);

% Remove trials with duration <30us
excl                    = trl.tEnd - trl.tOn < 30;
trl.tOn(excl)           = [];
trl.tEnd(excl)          = [];

% Trial loop
for iTrl = 1:length(trl.tEnd)
    
    trlIdx                  = [];
    trlIdx                  = d.time >= trl.tOn(iTrl) & d.time <= trl.tEnd(iTrl);   	% Trial index
    trl.dir_ts{iTrl}        = getTrialData(d.time, trlIdx, idx.RDP_dir);             	% Timestamps of RDP direction changes [Steady state start, Ignore first value -> randomly drawn direction]
    trl.SSno(iTrl)          = length(getTrialData(d.value, trlIdx, idx.RDP_dir));   	% Number of steady states
    trl.ssdur{iTrl}         = getTrialData(d.value, trlIdx, idx.steady_duration);
    
    % Steady state loop
%     for iSS = 1:trl.SSno(iTrl)
    for iSS = 1:trl.SSno(iTrl)-1
        
        % Steady state index
        if iSS < trl.SSno(iTrl)-1
            ssIdx           = d.time >= trl.dir_ts{iTrl}(iSS+1) & d.time < trl.dir_ts{iTrl}(iSS+2); % 1st timestamp irrelevant
        else
            ssIdx           = d.time >= trl.dir_ts{iTrl}(iSS+1) & d.time <= trl.tEnd(iTrl);
        end
        
        % Fill in table: Trial/State/Stimulus parameter
        cc                  = cc+1;                                                     % Steady state counter
        t.ID{cc}            = [monk '_' fid];                                           % Subject ID
        t.trl_no(cc)        = iTrl;                                                     % Trial number
        t.trl_dur(cc)       = (trl.tEnd(iTrl) - trl.tOn(iTrl)) ./ 1e6;                  % Trial duration [s]
        tmp_coh             = getTrialData(d.value, trlIdx, idx.RDP_coh);
        t.trl_coh(cc)       = tmp_coh(end);                                             % Trial coherence level
        t.ss_no(cc)         = cc;                                                       % Steady state counter
        t.ss_dur(cc)        = trl.ssdur{iTrl}(iSS);                                     % Steady state duration
        t.rdp_dir(cc)       = mod(getTrialData(d.value, ssIdx, idx.RDP_dir),360);       % Stimulus direction
        t.trg_shown(cc)     = iscell(getTrialData(d.value, ssIdx, idx.outcome));        % Target shown?
        t.trg_hit(cc)       = strcmp(getTrialData(d.value, ssIdx, idx.outcome), 'hit'); % Target collected?
%         t.trg_ts(cc)        = getTrialData(d.time, ssIdx, idx.outcome);                 % Target timestamp
        
        trg_val             = getTrialData(d.value, ssIdx, idx.trg); % Target timestamp
        trg_ts              = getTrialData(d.time, ssIdx, idx.trg); % Target timestamp
        if sum(trg_val) == 1
            t.trg_ts(cc) 	= trg_ts(logical(trg_val));
        else
            t.trg_ts(cc) 	= nan;
        end
        clear trg_val trg_ts
        
        % Fill in table: Behavioural response
        % (1) For each steady state...
        t.js_dir{cc}        = getTrialData(d.value, ssIdx, idx.JS_dir);                 % Joystick direction
        t.js_dir_ts{cc}     = getTrialData(d.time, ssIdx, idx.JS_dir);                  % Timestamps: Joystick direction
        t.js_str{cc}        = getTrialData(d.value, ssIdx, idx.JS_str);                 % Joystick strength
        t.js_str_ts{cc}     = getTrialData(d.time, ssIdx, idx.JS_str);                  % Timestamps: Joystick strength
        
        % (2) For entire trial...
        if iSS == 1
            t.trl_rdp_dir{cc}       = mod(getTrialData(d.value, trlIdx, idx.RDP_dir),360);
            t.trl_rdp_dir_ts{cc}    = getTrialData(d.time, trlIdx, idx.RDP_dir);
            t.trl_js_dir{cc}        = mod(getTrialData(d.value, trlIdx, idx.JS_dir),360);
            t.trl_js_dir_ts{cc}     = getTrialData(d.time, trlIdx, idx.JS_dir);
            t.trl_js_str{cc}        = getTrialData(d.value, trlIdx, idx.JS_str);
            t.trl_js_str_ts{cc}     = getTrialData(d.time, trlIdx, idx.JS_str);
        else
            t.trl_rdp_dir{cc}       = nan;
            t.trl_rdp_dir_ts{cc}    = nan;
            t.trl_js_dir{cc}        = nan;
            t.trl_js_dir_ts{cc}     = nan;
            t.trl_js_str{cc}        = nan;
            t.trl_js_str_ts{cc}     = nan;
        end
    end
end

t(ismissing(t.ID),:)                = [];
 
%% PERFORMANCE

f                   = figure('Units', 'normalized', 'Position', [0 0 .8 1]); set(gcf,'color', [1 1 1]);

HIidx               = t.trg_hit(t.trg_shown);
HIr                 = sum(HIidx) / length(HIidx);
MIr                 = sum(~HIidx) / length(HIidx);
trg_idx             = cell2mat(d.value(idx.trg));
trg_ts              = d.time(idx.trg);
trg_ts              = trg_ts(logical(trg_idx));
trg_ts              = (trg_ts - trg_ts(1)) ./ 1e6;
rew                 = d.value(idx.reward);
rew                 = cell2mat(rew(4:end));
rew_ts              = d.time(idx.reward);
rew_ts              = rew_ts(4:end);
rew_ts              = (rew_ts - rew_ts(1)) ./ 1e6;

%%% Performance pie chart
s                   = subplot(3,2,1);
p                   = pie([HIr,MIr], [1,1]);
pPatch              = findobj(p,'Type','patch');
pText               = findobj(p,'Type','text');
percentValues       = get(pText,'String');
txt                 = {'HIr: ';'MIr: '};
combinedtxt         = strcat(txt,percentValues);
s.Title.String      = {'Overall', 'Performance'};
s.Title.FontSize    = 14;
s.Title.Position(1) = -1;
s.Position(1)       = 0;
col                 = {[1 0 0],[0 0 0]};

for i = 1:size(pText,1)
    pText(i).String         = combinedtxt(i);
    pText(i).FontSize       = 14;
    pPatch(i).FaceColor     = col{i};
end

% Condition-wise performance
tcoh                = t.trl_coh(t.trg_shown);
clvl                = unique(tcoh);
for iCoh = 1:length(clvl)
    cindx           = tcoh == clvl(iCoh);
    hir(iCoh)    	= sum(HIidx(cindx)) / length(HIidx(cindx));
    mir(iCoh)     	= sum(~HIidx(cindx)) / length(HIidx(cindx));
end

ax                  = axes('Position',[.345 .75 .15 .15]); hold on
bp                	= bar([hir',mir'], 'stacked');
cl                  = [1 0];

for iOutc = 1:2
    bp(iOutc).FaceColor    	= [cl(iOutc) 0 0];
    bp(iOutc).EdgeColor     = [cl(iOutc) 0 0];
end

ax.XTick            = [1:length(clvl)];
% ax.XLim             = [0.5 5.5];
ax.XLim             = [0.5 length(clvl)+.5];
ax.YLabel.String   	= 'Rate';
ax.XLabel.String 	= 'Coherence level';
ax.XTickLabel    	= {round(clvl,2)};
ax.FontSize      	= 12;

%%% Performance over time
col                 = {[1 0 0], [0 0 0], [.5 .5 .5]};                   % Color specs

s                   = subplot(3,2,2); hold on
p                   = stairs(rew_ts,rew*1.2);
p.LineWidth         = 2;
p.LineStyle         = '-';
p.Color             = [1 1 1];
s.YLim           	= [0 rew(end)];
s.YLabel.String     = 'Cumulative reward [ml]';

x                  	= [p.XData(1),repelem(p.XData(2:end),2)];       % Fill area
y                   = [repelem(p.YData(1:end-1),2),p.YData(end)];
fl                  = fill([x,fliplr(x)],[y,0*ones(size(y))], col{3});
fl.FaceAlpha        = .5;
fl.EdgeAlpha        = .5;
fl.FaceColor        = col{3};
fl.EdgeColor        = col{3};

yyaxis right
p2                  = plot(trg_ts,movmean(HIidx,5));
p2.LineWidth        = 2;
p2.Color            = col{1};
p3                  = plot(trg_ts,movmean(~HIidx,5));
p3.LineWidth        = 2;
p3.LineStyle        = ':';
p3.Color            = col{2};
s.YLim           	= [0 1];
s.XLim           	= [1 rew_ts(end)];
s.Title.String      = 'Performance [movmean, 5 targets]';
s.XLabel.String     = 'Time [s]';
s.YLabel.String     = 'Rate';
s.FontSize          = 14;
s.Box           	= 'off';
s.YAxis(1).Color    = col{3};
s.YAxis(2).Color    = col{1};
l                   = legend([p2,p3,fl],{'HIr', 'MIr', 'ml'});
l.Location          = 'west';

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
    
    out.str_mean_dist{iCoh}	= mean(str_arr,2);
    out.str_std_dist{iCoh} 	= std(str_arr,[],2);
    out.acc_dist{iCoh}    	= mean(js_acc{iCoh},2);
    
    out.str_mean(iCoh)      = mean(mean(str_arr,2));                               	% Average across steady states
    out.str_std(iCoh)       = std(std(str_arr,[],2));                                % Standard deviation across steady states
    out.HIr(iCoh)           = sum(trg_hit) / sum(trg_shown);                         % Hit rate
    out.acc(iCoh)           = mean(mean(js_acc{iCoh},2));                            % Response accuracy
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


%%% TRIAL DATA %%%
tmp = cellfun(@size, t.trl_rdp_dir, 'UniformOutput', false);
for i = 1:size(tmp,1)
    indx(i,:) = tmp{i}(2) > 1;
end
coh                         = t.trl_coh(indx);
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

out.xc                      = xc;
out.maxR                    = maxR;
out.posPk                   = posPk;
out.auPk                    = auPk;
out.coh                     = coh;
                
ax.XTick                    = [1 nLag/2 nLag];
ax.XLim                     = [1 nLag];
ax.XTickLabel               = {['-' num2str(nLag)],['-' num2str(nLag/2)],'0'};
ax.XLabel.String            = 'Lag';
ax.YLabel.String            = 'XCorr coeff';
ax.FontSize                 = 14;
ax.Title.String             = 'XCorr(RDP, JS)';
ax.Title.Interpreter        = 'none';
yyaxis right
pm                          = plot(mean(xc),'Color',[1 0 0], 'Marker','none', 'LineWidth',3);
ax.YAxis(2).Color           = [1 0 0];
lg                          = legend([ps, pm], {'Trials','Mean'});
lg.Position(1)              = .15;

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
    xcLag{iCoh}             = [xcLag{iCoh} posPk(cIdx)'-nLag];
    
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
ax.XTickLabel               = {rsnr};
ax.Title.String             = 'Area under XCorr peak';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
box off

ax                        	= subplot(3,3,9); hold on
vs                        	= violinplot(bby,bbx);
cl                        	= linspace(0,1,size(vs,2));
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor    = [cl(iSub) 0 0];
end
ax.YLabel.String            = 'Lag';
ax.XLabel.String            = 'Coherence level';
ax.XTickLabel               = {rsnr};
ax.Title.String             = 'XCorr peak lag';
ax.XColor                   = [0 0 0];
ax.YColor                   = [0 0 0];
ax.FontSize                 = 14;
box off

dest_dir = '/Users/fschneider/Documents/MWorks/XBI/Plots/';
print(f, [dest_dir 'summary_' monk '_' fid], '-r300', '-dpng');

end
