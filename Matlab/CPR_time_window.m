function out = CPR_time_window(t,nSamples)

% This function extracts stimulus and joystick data in a specified time
% window at the end of a steady state
%
% Input:  	.t              Table, Contains steady state information
%          	.nSamples       Integer, Number of samples before direction
%                           changes (End of state)
%
% Output:   .out            Structure, Contains averages and performance
%                           measures for each coherence level
%
% Example:	out = CPR_time_window(tbl,29)
%
% Known bugs:
%
% Feature wish list:
%
% Version history
%   1.0     (fxs 2020-05-04) Initial version.

%% Check input

addpath /Users/fschneider/Documents/MATLAB/CircStat2012a

if nargin < 2
    nSamples            = 25;
end

if nargin < 1 || ~istable(t)
    error('Input must be a table')
end

%% Time window analysis

cohPool                 = unique(t.ss_coh);                                         % Tested coherence levels
cohPool(cohPool == 0)   = [];                                                   	% Zero not tested, throw out

for iCoh = 1:size(cohPool,1)                                                        % For each coherence level
    str_arr           	= [];                                                       % Reset temporary variables
    dir_arr             = [];
    trg_shown           = [];
    trg_hit             = [];
    
    cohIdx              = t.ss_coh == cohPool(iCoh);                                % Coherence index
    str_raw{iCoh}       = t.js_str(cohIdx);                                         % Joystick strength
    dir_raw{iCoh}       = t.js_dir(cohIdx);                                         % Joystick direction
    rdp_dir{iCoh}       = t.rdp_dir(cohIdx);                                        % Stimulus direction
    
    for iSS = 1:size(str_raw{iCoh},1)                                               % For each steady state...
        str_arr(iSS,:)  	= str_raw{iCoh}{iSS}(end-nSamples:end);                 % Take last nSamples+1 data points [Fs = 100Hz -> 10ms] of each steady state
        dir_arr(iSS,:)   	= dir_raw{iCoh}{iSS}(end-nSamples:end);
    end
    
    js_dev{iCoh}            = rad2deg(circ_dist(deg2rad(dir_arr),deg2rad(t.rdp_dir(cohIdx))));  % Minimum RDP-Joystick difference
    js_acc{iCoh}            = abs(1 - abs(js_dev{iCoh}) / 180);                                 % Joystick accuracy
    trg_hit                	= t.trg_hit(cohIdx);
    trg_shown            	= t.trg_shown(cohIdx);
    
    out.str_mean_dist{iCoh}	= nanmedian(str_arr,2);
    out.str_std_dist{iCoh} 	= nanstd(str_arr,[],2);
    out.acc_dist{iCoh}    	= nanmedian(js_acc{iCoh},2);
    
    out.str_mean(iCoh)      = nanmedian(nanmedian(str_arr,2));                  	% Average across steady states
    out.str_std(iCoh)       = nanstd(nanstd(str_arr,[],2));                       	% Standard deviation across steady states
    out.HIr(iCoh)           = sum(trg_hit) / sum(trg_shown);                     	% Hit rate
    out.acc(iCoh)           = nanmedian(nanmedian(js_acc{iCoh},2));               	% Response accuracy
end
end