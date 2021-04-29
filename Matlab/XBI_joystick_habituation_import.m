function [d, idx, trl] = XBI_joystick_habituation_import(fname, d)

% Add relevant directories
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/ownCloud/Documents/Matlab/Helper_functions

cd('/Users/fschneider/Documents/MWorks/XBI')

close all

%% Check input arguments & import data

% Define variables for import
var_import            	= {'ML_', 'CTRL_', 'TRIAL_', 'RDP_'};

if nargin < 2                                                               % If no data structure provided...
    if iscell(fname)                                                        % If multiple files ...
        d               = mergeFiles(fname, var_import);                    % ... merge to single structure
    else
        if isfile([fname '.mat'])                                           % If .mat file available...
            tmp         = load([fname '.mat']);                            	% ...load .mat file
            d           = tmp.d;
        else
            d           = MW_readFile(fname, 'include', var_import);       	% Import .mwk2 session file
            
%             disp('Save struct...')
%             save([fname '.mat'], 'd', '-v7.3')                              % Save as .mat file
%             disp('Done!')
        end
    end
end

if iscell(fname)                                                            % Get session ID
    fid                	= fname{1}(38:45);
else
    fid               	= fname(38:45);
end

%% % Get data

lst                     = categories(d.event);                              % List of variables, just fyi
idx                     = [];                                            	% Initialise
trl                     = [];

% Create indices
idx.tOn                 = d.event == 'ML_trialStart';
idx.tEnd                = d.event == 'ML_trialEnd';
idx.response            = d.event == 'CTRL_time_response';
idx.previous            = d.event == 'CTRL_time_previous_response';
idx.RT                  = d.event == 'TRIAL_RT';
idx.RespDir             = d.event == 'TRIAL_RespDir';
idx.outcome             = d.event == 'TRIAL_outcome';

% Trial timestamps
trl.tOn                 = d.time(idx.tOn);
trl.tEnd                = d.time(idx.tEnd);

% Remove last trial onset if unequal number
if size(trl.tOn,2)      ~= size(trl.tEnd,2)
    trl.tOn(end)        = [];
end

tmp                     = trl.tEnd - trl.tOn < 5;
trl.tOn(tmp)            = [];
trl.tEnd(tmp)           = [];

%Initialise
trl.previous            = nan(1,length(trl.tEnd));
trl.response          	= nan(1,length(trl.tEnd));
trl.RT                  = nan(1,length(trl.tEnd));
trl.RespDir             = nan(1,length(trl.tEnd));
trl.outcome             = nan(1,length(trl.tEnd));

for iTrl = 1:length(trl.tEnd)
    trlIdx              = [];
    trlIdx              = d.time > trl.tOn(iTrl) & d.time < trl.tEnd(iTrl);
    trl.previous(iTrl)	= getTrialData(d.value, trlIdx, idx.previous);
    trl.response(iTrl)	= getTrialData(d.value, trlIdx, idx.response);
    trl.RT(iTrl)        = getTrialData(d.value, trlIdx, idx.RT);
    trl.RespDir(iTrl)   = getTrialData(d.value, trlIdx, idx.RespDir);
    tmp_outc            = getTrialData(d.value, trlIdx, idx.outcome);
    
    if iscell(tmp_outc)
        if strcmp(tmp_outc, 'hit') || strcmp(tmp_outc, 'success')
            trl.outcome(iTrl)	= true;
        else
            trl.outcome(iTrl)   = false;
        end
    else
        trl.outcome(iTrl)       = tmp_outc;
    end
end

trl.RT(1)               = nan;

f                       = figure;
left_color              = [0 0 0];
right_color             = [1 0 0];
set(f,'defaultAxesColorOrder',[left_color; right_color]);
bwidth                  = .05;
foSize                  = 12;

% Distribution of response directions
s1                      = subplot(2,2,3);
h1                      = histogram(trl.RespDir);
s1.XLim                 = [0 360];
s1.Title.String         = 'Direction Joystick Response';
s1.XLabel.String        = 'Angle [deg]';
s1.YLabel.String        = 'No. Trials';
s1.FontSize             = foSize;
s1.Box                  = 'off';
h1.BinWidth             = s1.XLim(2) * bwidth;
h1.FaceColor            = 'k';
h1.EdgeColor            = 'k';
h1.FaceAlpha            = 1;

% Inter-response interval
s2                      = subplot(2,2,1:2); hold on
tmpRT                   = trl.RT;
tmpRT(tmpRT > 10000)    = 10000;
sc21                    = scatter(1:size(tmpRT(trl.outcome == 0),2),tmpRT(trl.outcome == 0));
sc21.MarkerFaceColor    = [.75 .75 .75];
sc21.MarkerEdgeColor    = [.75 .75 .75];
sc22                    = scatter(1:size(tmpRT(trl.outcome == 1),2),tmpRT(trl.outcome == 1));
sc22.MarkerFaceColor    = [0 0 0];
sc22.MarkerEdgeColor    = [0 0 0];
ff                      = fit((2:length(tmpRT))', tmpRT(2:end)', 'poly1');
ls2                     = plot(ff);
ls2.LineWidth           = 2;
ls2.Color               = [1 0 0];
lg                      = legend([ls2,sc21,sc22],{'lsline','error','correct'});
lg.Box                  = 'on';
lg.FontSize             = 12;
s2.Title.String         = 'Inter-response interval';
s2.XLabel.String        = 'No. Trial';
s2.YLabel.String        = 'IRI [ms]';
s2.FontSize             = foSize;
s2.Box                  = 'off';
s2.YLim                 = [0 max(tmpRT)];
s2.XLim                 = [1 length(trl.RT)];

s3                      = subplot(2,2,4); hold on
h3                      = histogram(tmpRT(trl.outcome == 1),20);
h3.FaceColor            = 'k';
h3.EdgeColor            = 'k';
h3.FaceAlpha            = 1;
s3.Title.String         = 'Inter-response interval [correct]';
s3.YLabel.String        = 'No. Trial';
s3.XLabel.String        = 'IRI [ms]';
s3.FontSize             = foSize;
s3.XLim                 = [1000 ceil(h3.BinEdges(end))];

print(['./Plots/summary_nil_' fname(38:45)], '-dpng', '-r300')
end