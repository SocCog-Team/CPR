function d = XBI_arc_target_import(monk, fname, d)

% Add relevant directories
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab

cd('/Users/fschneider/Documents/MWorks/XBI')
% cd('/Users/fschneider/Documents/MWorks/Data')

close all

%% Check input arguments & import data

% Define variables for import
var_import = {'ML_', 'CTRL_', 'TRIAL_', 'RDP_','IO_joystickDirection'};

if nargin < 3                                                              
    d                	= CPR_import_mwk2(fname, var_import, false);
end

if iscell(fname)                                                            % Get session ID
    fid                	= fname{1}(31:38);
else
    fid               	= fname(31:38);
end

%% PLOT

f                       = figure('Units', 'normalized', 'Position', [0 0 1 1]);

%%% RT distribution
ax                      = subplot(2,8,1:4); hold on
[RT, ~, ~, ~]        	= MW_getRT(d, [], 30, [], [1 0], ax);                      	% Get reaction time distributions

%%% Conditionwise RT
ax                      = subplot(2,8,5:8); hold on
[RT, ~, ~, ~]        	= MW_getRT(d, [], 30, [], [0 1], ax);                     	% Get reaction time distributions
ax.Position(2)          = 0.5838; 
ax.Position(4)          = 0.3412; 

%%% Performance pie chart
ax                    	= subplot(2,8,9); hold on
[perf, pIdx, ~]      	= MW_getPerformance(d,50,[1 0 0 0 0],ax);

%%% Moving hit rate
ax                   	= subplot(2,8,11:12); hold on
[perf, pIdx, ~]        	= MW_getPerformance(d,50,[0 1 0 0 0],ax);
ax.Position             = [0.285    0.1100    0.2    0.3412]; 

%%% Condition-wise performance
ax                     	= subplot(2,8,13:16); hold on
[perf, pIdx, ~]         = MW_getPerformance(d,50,[0 0 1 0 0],ax);

dest_dir = '/Users/fschneider/Documents/MWorks/XBI/Plots/';
print(f, [dest_dir 'summary_ ' monk '_' fid], '-r300', '-dpng');
 

%% Directional response

% Create indices
idx.tOn         	= d.event == 'TRIAL_start';
idx.tEnd        	= d.event == 'TRIAL_end';
idx.js_dir          = d.event == 'IO_joystickDirection';
idx.rdp_dir         = d.event == 'RDP_direction';
idx.outcome         = d.event == 'TRIAL_outcome';
idx.trg             = d.event == 'TRIAL_reactionTrigger';

if sum(idx.tOn) == 0
    idx.tOn     	= d.event == 'ML_trialStart';
    idx.tEnd        = d.event == 'ML_trialEnd';
end

% Trial timestamps
trl.tOn             = d.time(idx.tOn);
trl.tEnd            = d.time(idx.tEnd);

% Remove last trial onset if unequal number
if size(trl.tOn,2)	~= size(trl.tEnd,2)
    trl.tOn(end)	= [];
end

tmp             	= trl.tEnd - trl.tOn < 5;
trl.tOn(tmp)     	= [];
trl.tEnd(tmp)   	= [];

for iTrl = 1:length(trl.tEnd)
    
    % Extract trial data
    trlIdx              = [];
    trlIdx              = d.time > trl.tOn(iTrl) & d.time < trl.tEnd(iTrl);
    trl.js_dir{iTrl}	= getTrialData(d.value, trlIdx, idx.js_dir);
    trl.rdp_dir{iTrl}	= getTrialData(d.value, trlIdx, idx.rdp_dir);
    trl.js_dir_ts{iTrl}	= getTrialData(d.time, trlIdx, idx.js_dir);
    trl.rdp_dir_ts{iTrl}= getTrialData(d.time, trlIdx, idx.rdp_dir);
    trl.outcome(iTrl)   = getTrialData(d.value, trlIdx, idx.outcome);
    trl.trg_ts{iTrl}    = getTrialData(d.time, trlIdx, idx.trg);
    
    dur = 29;
    if size(trl.rdp_dir{iTrl},2) == 1
        tmp = []; js_dev = [];js_acc = []; acc_tmp = [];
        indx                = trl.js_dir_ts{iTrl} <= trl.trg_ts{iTrl}(1);
        tmp                 = trl.js_dir{iTrl}(indx);
        if size(tmp,2) > dur
            js_dev              = tmp(end-dur:end) - trl.rdp_dir{iTrl};
            js_acc              = abs(1 - abs(mod(js_dev,360)) / 180);
            trl.avg_acc(iTrl)   = mean(js_acc);
        else
            trl.avg_acc(iTrl)   = nan;
        end
    else
        for iDir = 1:size(trl.rdp_dir{iTrl},2)
            tmp = []; js_dev = []; js_acc = []; acc_tmp = [];
            if iDir == size(trl.rdp_dir{iTrl},2)
                indx        = trl.js_dir_ts{iTrl} <= trl.trg_ts{iTrl}(1);
            else
                indx        = trl.js_dir_ts{iTrl} <= trl.rdp_dir_ts{iTrl}(iDir+1);
            end
            
            tmp             = trl.js_dir{iTrl}(indx);
            js_dev          = tmp(end-dur:end) - trl.rdp_dir{iTrl}(iDir);
            js_acc          = abs(1 - abs(mod(js_dev,360)) / 180);
            acc_tmp(iDir)   = mean(js_acc);
        end
        trl.avg_acc(iTrl)   = mean(acc_tmp);
    end
end

if ~isnan(nanmean(trl.avg_acc(iTrl)))
    %%% PLOT ACCURACY
    addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
    
    steCnt              = cellfun(@length,trl.rdp_dir);
    outc                = cell2mat(cellfun(@(c)strcmp(c,'hit'),trl.outcome, 'UniformOutput', false));
    
    f                   = figure;
    vs                  = violinplot(trl.avg_acc(outc),steCnt(outc));
    ax                  = gca;
    ax.XLabel.String    = 'No. States before Target';
    ax.YLabel.String    = 'Accuracy [norm]';
    ax.Title.String     = {'Average motion tracking accuracy', '[300ms window before target presentation]'};
    ax.FontSize         = 14;
    
    print(f, [dest_dir 'accuracy_' monk '_' fid], '-r300', '-dpng');
end
end