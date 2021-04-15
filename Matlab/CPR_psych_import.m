function [t,d, out] = CPR_psych_import(subj, pth, fname, d)

% Add relevant directories
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

cd(pth)

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


if nargin < 4                                                               % If no data structure provided...
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

fid                     = fname(9:16);

%% Initialise variables

% Counter
cc                      = 0;
idx                     = [];
trl                     = [];

% Table specs
var     = {'ID', 'string'; ...
    'trl_no', 'double'; ...
    'trl_dur', 'double'; ...
    'ss_no', 'double'; ...
    'ss_coh', 'double'; ...
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
    'trl_rdp_coh', 'cell'; ...
    'trl_rdp_coh_ts', 'cell'; ...
    'trl_js_dir', 'cell'; ...
    'trl_js_dir_ts', 'cell'; ...
    'trl_js_str', 'cell'; ...
    'trl_js_str_ts', 'cell'};

% Initialise table
t       = table('Size',[2000, size(var,1)],...
    'VariableTypes',var(:,2),...
    'VariableNames',var(:,1));

%% Extract data

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
idx.reward                  = d.event == 'INFO_Cash';
idx.stateOn                 = d.event == 'INFO_SteadyStateCounter';

if sum(idx.tOn) == 0
    idx.tOn                 = d.event == 'TRIAL_start';
    idx.tEnd                = d.event == 'TRIAL_end';
end

% Trial timestamps
trl.tOn                     = d.time(idx.tOn);
trl.tEnd                    = d.time(idx.tEnd);

% Remove trials with duration <30us
excl                        = trl.tEnd - trl.tOn < 30;
trl.tOn(excl)               = [];
trl.tEnd(excl)              = [];

% Trial loop
for iTrl = 1:length(trl.tEnd)
    
    trlIdx                  = [];
    trlIdx                  = d.time >= trl.tOn(iTrl) & d.time <= trl.tEnd(iTrl);   	% Trial index
    trl.dir_ts{iTrl}        = getTrialData(d.time, trlIdx, idx.stateOn);             	% Timestamps of RDP direction changes [Steady state start, Ignore first value -> randomly drawn direction]
    trl.SSno(iTrl)          = length(getTrialData(d.value, trlIdx, idx.RDP_dir))-1;     % Number of steady states
    trl.ssdur{iTrl}         = getTrialData(d.value, trlIdx, idx.steady_duration);
    trl.coh{iTrl}           = getTrialData(d.value, trlIdx, idx.RDP_coh);               % Trial coherence level

    % Steady state loop
    for iSS = 1:trl.SSno(iTrl)
        
        % Steady state index
        if iSS < trl.SSno(iTrl)
            ssIdx           = d.time >= trl.dir_ts{iTrl}(iSS+1) & d.time < trl.dir_ts{iTrl}(iSS+2);
        else
            ssIdx           = d.time >= trl.dir_ts{iTrl}(iSS+1) & d.time <= trl.tEnd(iTrl);
        end
        
        % Fill in table: Trial/State/Stimulus parameter
        cc                  = cc+1;                                                     % Steady state counter
        t.ID{cc}            = [subj '_' fid];                                           % Subject ID
        t.trl_no(cc)        = iTrl;                                                     % Trial number
        t.trl_dur(cc)       = (trl.tEnd(iTrl) - trl.tOn(iTrl)) ./ 1e6;                  % Trial duration [s]
        t.ss_no(cc)         = cc;                                                       % Steady state counter
        t.ss_coh(cc)        = getTrialData(d.value, ssIdx, idx.RDP_coh);                % Steady state coherence
        t.ss_dur(cc)        = trl.ssdur{iTrl}(iSS);                                     % Steady state duration
        t.rdp_dir(cc)       = mod(getTrialData(d.value, ssIdx, idx.RDP_dir),360);       % Stimulus direction
        t.trg_shown(cc)     = iscell(getTrialData(d.value, ssIdx, idx.outcome));        % Target shown?
        t.trg_hit(cc)       = strcmp(getTrialData(d.value, ssIdx, idx.outcome), 'hit'); % Target collected?        
        trg_val             = getTrialData(d.value, ssIdx, idx.trg);                    % Target trigger
        trg_ts              = getTrialData(d.time, ssIdx, idx.trg);                     % Target timestamp
        if sum(trg_val) == 1
            t.trg_ts(cc)  	= trg_ts(logical(trg_val));
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
            t.trl_rdp_coh{cc}       = getTrialData(d.value, trlIdx, idx.RDP_coh);
            t.trl_rdp_coh_ts{cc}    = getTrialData(d.time, trlIdx, idx.RDP_coh);
            t.trl_js_dir{cc}        = mod(getTrialData(d.value, trlIdx, idx.JS_dir),360);
            t.trl_js_dir_ts{cc}     = getTrialData(d.time, trlIdx, idx.JS_dir);
            t.trl_js_str{cc}        = getTrialData(d.value, trlIdx, idx.JS_str);
            t.trl_js_str_ts{cc}     = getTrialData(d.time, trlIdx, idx.JS_str);
        else
            t.trl_rdp_dir{cc}       = nan;
            t.trl_rdp_dir_ts{cc}    = nan;
            t.trl_rdp_coh{cc}       = nan;
            t.trl_rdp_coh_ts{cc}    = nan;
            t.trl_js_dir{cc}        = nan;
            t.trl_js_dir_ts{cc}     = nan;
            t.trl_js_str{cc}        = nan;
            t.trl_js_str_ts{cc}     = nan;
        end   
    end
    
    % Temporary fix: coherence per steady state
    b = repmat(trl.coh{iTrl},10,1);
    if iTrl == 1
        c = reshape(b(:,2:end),size(b,1)*(size(b,2)-1),1);
    else
        c = reshape(b,size(b,1)*(size(b,2)),1);
    end
    strt = (iTrl*100)-99;
    t.ss_coh(strt:strt+length(c)-1) = c;
end

t(ismissing(t.ID),:)                = [];
disp('Save table...')
save([fname '_tbl.mat'], 't', '-v7.3')                                                   % Save as .mat file
disp('Done!')
            
%% PERFORMANCE

f                   = figure('Units', 'normalized', 'Position', [0 0 .8 1]); set(gcf,'color', [1 1 1]);

HIidx               = t.trg_hit(t.trg_shown);
HIr                 = sum(HIidx) / length(HIidx);
MIr                 = sum(~HIidx) / length(HIidx);
trg_ts              = t.trg_ts(~isnan(t.trg_ts));
trg_ts              = (trg_ts - trl.tOn(1)) ./ 1e6;
rew                 = d.value(idx.reward);
exIdx               = strcmp(cellfun(@class, rew, 'UniformOutput', false),'double');
rew                 = cell2mat(rew(exIdx));
rew_ts              = d.time(idx.reward);
rew_ts              = rew_ts(exIdx);
rew_ts              = (rew_ts - trl.tOn(1)) ./ 1e6;

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
tcoh                = t.ss_coh(t.trg_shown);
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
ax.XLim             = [0.5 length(clvl)+.5];
ax.YLabel.String   	= 'Rate';
ax.XLabel.String 	= 'Coherence level';
ax.XTickLabel    	= {round(clvl,2)};
ax.FontSize      	= 12;

%%% Performance over time
col                 = {[1 0 0], [0 0 0], [.5 .5 .5]};                   % Color specs

s                   = subplot(3,2,2); hold on
p                   = stairs(rew_ts,rew);
p.LineWidth         = 2;
p.LineStyle         = '-';
p.Color             = [1 1 1];
s.YLim           	= [0 rew(end)];
s.YLabel.String     = 'Cumulative reward [EUR]';

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
l                   = legend([p2,p3,fl],{'HIr', 'MIr', 'EUR'});
l.Location          = 'west';

%% ANALYSIS OF TIME WINDOW

nSamples                = 29;                                                   	% Number of samples before direction changes
cohPool                 = unique(t.ss_coh);                                         % Tested coherence levels
cohPool(cohPool == 0) = [];                                                         % Tested coherence levels

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


%%% TRIAL DATA %%%
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

%% XCORR ANALYSIS

ax                  = subplot(3,3,7); hold on
nLag                = 150;                                                  % XCorr Lag
rdp_dir             = t.trl_rdp_dir(indx);                                  % RDP direction
rdp_dir_ts          = t.trl_rdp_dir_ts(indx);
js_dir              = t.trl_js_dir(indx);                                   % Joystick direction
js_dir_ts           = t.trl_js_dir_ts(indx);
% js_str              = t.trl_js_str(indx);                                 % Joystick strength
% js_str_ts           = t.trl_js_str_ts(indx);
rdp_coh            	= t.trl_rdp_coh(indx);                                  % RDP coherence
rdp_coh_ts        	= t.trl_rdp_coh_ts(indx);
rdp_coh{1}(1)   	= [];                                                   % Remove first entry of first trial
rdp_coh_ts{1}(1)  	= [];
count               = 0;

for iTrl = 1:size(rdp_dir_ts,1)
    
    % Adjust vector length
    dir_vec                     = [];   
    for iSS = 1:size(rdp_dir_ts{iTrl},2)-1
        ssIdx                   = [];
        ssIdx                   = js_dir_ts{iTrl} >= rdp_dir_ts{iTrl}(iSS) & js_dir_ts{iTrl} < rdp_dir_ts{iTrl}(iSS+1);
        dir_vec               	= [dir_vec repmat(rdp_dir{iTrl}(iSS),1,sum(ssIdx))];  
    end
    
    ssIdx                       = js_dir_ts{iTrl} >= rdp_dir_ts{iTrl}(end);
    dir_vec                     = [dir_vec repmat(rdp_dir{iTrl}(end),1,sum(ssIdx))];
    
    if size(dir_vec,2) ~= size(js_dir{iTrl},2)
        dir_vec                 = [dir_vec nan(1,size(js_dir{iTrl},2) - size(dir_vec,2))];
    end
    
    % Correct for circular space
    clear js_dff js_corr rdp_dff rdp_corr
    
    js_dff                      = [0 mod(diff(js_dir{iTrl}) + 180, 360) - 180];
    js_corr                     = js_dff(1);
    rdp_dff                     = [0 mod((diff(dir_vec) + 180), 360) - 180];
    rdp_corr                    = rdp_dff(1);
    
    for iSample = 1:length(js_dir{iTrl})-1
        js_corr(iSample+1)      = js_dff(iSample) + js_corr(iSample);
        rdp_corr(iSample+1)     = rdp_dff(iSample) + rdp_corr(iSample);
    end
     
    % Extract coherence chunks
    for iCoh  = 1:size(rdp_coh{iTrl},2)
        
        % Build coherence index
        if iCoh < size(rdp_coh{iTrl},2)
            cIdx                = js_dir_ts{iTrl} >= rdp_coh_ts{iTrl}(iCoh) & js_dir_ts{iTrl} < rdp_coh_ts{iTrl}(iCoh+1);
        else
            cIdx                = js_dir_ts{iTrl} >= rdp_coh_ts{iTrl}(iCoh) & js_dir_ts{iTrl} <= js_dir_ts{iTrl}(end);
        end
        
        % Cross-correlation analysis
        count                   = count + 1;
        coh(count)              = rdp_coh{iTrl}(iCoh);                                          % Coherence ID
        xc(count,:)             = xcorr(double(abs(rdp_dff(cIdx))),abs(js_dff(cIdx)),nLag);    	% Crosscorrelation between stimulus and joystick direction
        maxR(count)           	= max(xc(iTrl,1:nLag+1));                                      	% Max correlation coefficient
        posPk(count)        	= find(xc(count,1:nLag+1) == max(xc(count,1:nLag+1)));           	% Peak position of cross-correlation
    
        try
            auPk(count,:)     	= trapz(xc(count,posPk(count)-10:posPk(count)+10));           	% Area under cross-correlation peak
        catch
            auPk(count,:)     	= nan;
        end
        
        % Plot trial-wise cross-correlation 
        ps                      = plot(xc(count,1:nLag),'Color',[.5 .5 .5 .1], 'Marker','none', 'LineWidth',2);    
        
    end
end

cohID(cohID == 0)               = nan;
out.xc                          = xc;
out.maxR                        = maxR;
out.posPk                       = posPk;
out.auPk                        = auPk;
out.coh                         = coh;
                
ax.XTick                        = [1 nLag/2 nLag];
ax.XLim                         = [1 nLag];
ax.XTickLabel                   = {['-' num2str(nLag)],['-' num2str(nLag/2)],'0'};
ax.XLabel.String                = 'Lag';
ax.YLabel.String                = 'XCorr coeff';
ax.FontSize                     = 14;
ax.Title.String                 = 'XCorr(RDP, JS)';
ax.Title.Interpreter            = 'none';
yyaxis right

pm                              = plot(nanmedian(xc),'Color',[1 0 0], 'Marker','none', 'LineWidth',4);
ax.YAxis(2).Color               = [1 0 0];
lg                              = legend([ps, pm], {'Trials','Median'});
lg.Position(1)                  = .15;

% if max(median(xc(:,1:nLag))) < 2 * median(std(xc(:,1:nLag)))
%     ax.Title.String     	= 'XCorr(RDP, JS) -> n.s.';
% end

%%% Coherence %%%
bx                              = [];
by                              = [];
bbx                             = [];
bby                             = [];

snr                             = unique(coh);
aucPeak                         = cell(1,length(snr));
xcLag                           = cell(1,length(snr));

for iCoh = 1:size(snr,2)
    cIdx                        = coh == snr(iCoh);
    aucPeak{iCoh}               = [aucPeak{iCoh} auPk(cIdx)];
    xcLag{iCoh}                 = [xcLag{iCoh} posPk(cIdx)'-nLag];
        
    by                          = [by; aucPeak{iCoh}];
    bx                          = [bx; repmat(iCoh,length(aucPeak{iCoh}),1)];
    
    bby                         = [bby; xcLag{iCoh}];
    bbx                         = [bbx; repmat(iCoh,length(xcLag{iCoh}),1)];
end

ax                              = subplot(3,3,8); hold on
vs                              = violinplot(by,bx);
cl                              = linspace(0,1,size(vs,2));
for iSub = 1:size(vs,2)
    vs(iSub).ViolinColor        = [cl(iSub) 0 0];
end
ax.YLabel.String                = 'Area under peak';
ax.XLabel.String                = 'Coherence level';
ax.XTickLabel                   = {rsnr};
ax.Title.String                 = 'Area under XCorr peak';
ax.XColor                       = [0 0 0];
ax.YColor                       = [0 0 0];
ax.FontSize                     = 14;
box off

ax                              = subplot(3,3,9); hold on
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

print(f, [pth '/summary_' subj '_' fid], '-r300', '-dpng');

end
