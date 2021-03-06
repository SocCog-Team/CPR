function [out, ps] = CPR_correlation(t, trl, tindx, nLag, plotFlag)

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
%   1.0     (fxs 2021-05-04) Initial version.
%   1.1     (fxs 2021-05-25) Improved vectorisation procedure of data

%% Check input

addpath /Users/fschneider/Documents/MATLAB/CircStat2012a

if nargin < 4 || isempty(plotFlag)
    plotFlag     	= true;
end

if nargin < 3 || isempty(nLag)
    nLag            = 150;
end

if nargin < 2 || isempty(tindx)
    error('Index pointing to trial/block data is missing')
end

if nargin < 1 || ~istable(t)
    error('Input must be a table')
end

%% Assign data from table

rdp_dir             = t.trl_rdp_dir(tindx);                              	% RDP direction
rdp_dir_ts          = t.trl_rdp_dir_ts(tindx);
js_dir              = t.trl_js_dir(tindx);                               	% Joystick direction
js_ts             	= t.trl_js_ts(tindx);
% js_str              = t.trl_js_str(tindx);                                % Joystick strength
rdp_coh            	= t.trl_rdp_coh(tindx);                               	% RDP coherence
rdp_coh_ts        	= t.trl_rdp_coh_ts(tindx);
% rdp_coh{1}(1)   	= [];                                                   % Remove first entry of first trial
% rdp_coh_ts{1}(1)  = [];
count               = 0;

% excl                = cell2mat(cellfun(@(x) sum(x == 0),rdp_coh,'UniformOutput', false)) > 0;
% rdp_coh(excl)       = [];
% rdp_coh_ts(excl)    = [];

%% Loop through blocks and correlate stimulus and behaviour

for iTrl = 1:size(rdp_dir_ts,1)
    
    % Adjust vector length
    rdp_vec                     = []; 
    js_vec                      = []; 
    ts                          = []; 
    
    for iDir = 1:size(rdp_dir{iTrl},2)-1
        % JS index
        jsIdx                   = [];
        jsIdx               	= js_ts{iTrl} >= rdp_dir_ts{iTrl}(iDir) & js_ts{iTrl} < rdp_dir_ts{iTrl}(iDir+1);
        js_vec               	= [js_vec js_dir{iTrl}(jsIdx)];
        ts                      = [ts js_ts{iTrl}(jsIdx)];
        rdp_vec              	= [rdp_vec repmat(rdp_dir{iTrl}(iDir),1,sum(jsIdx))];
    end
    
    jsIdx                       = [];
    jsIdx                       = js_ts{iTrl} >= rdp_dir_ts{iTrl}(end) & js_ts{iTrl} <= trl.tEnd(iTrl);
    js_vec                      = [js_vec js_dir{iTrl}(jsIdx)];
    ts                          = [ts js_ts{iTrl}(jsIdx)];
    rdp_vec                     = [rdp_vec repmat(rdp_dir{iTrl}(end),1,sum(jsIdx))];
    
    % Sample-by-sample difference
    clear js_dff js_corr rdp_dff rdp_corr
    
    for iSample = 1:size(js_vec,2)-1
        js_dff(iSample)     	= rad2deg(circ_dist(deg2rad(js_vec(iSample)),deg2rad(js_vec(iSample+1))));
        rdp_dff(iSample)     	= rad2deg(circ_dist(deg2rad(rdp_vec(iSample)),deg2rad(rdp_vec(iSample+1))));
    end
 
    js_dff                      = [0 js_dff];
    rdp_dff                     = [0 rdp_dff];

%     js_dff                      = [0 mod(diff(js_dir{iTrl}) + 180, 360) - 180];
%     rdp_dff                     = [0 mod((diff(rdp_vec) + 180), 360) - 180];
%
%     % For plotting: Unwrap data [cumulative sample-by-sample difference]
%     js_corr                     = js_dff(1);
%     rdp_corr                    = rdp_dff(1);
%     
%     for iSample = 1:length(js_dir{iTrl})-1
%         js_corr(iSample+1)      = js_dff(iSample) + js_corr(iSample);
%         rdp_corr(iSample+1)     = rdp_dff(iSample) + rdp_corr(iSample);
%     end
     
    % Extract coherence chunks
    cindx                       = [true diff(rdp_coh{iTrl}) ~=0];
    cts                         = rdp_coh_ts{iTrl}(cindx);
    
    for iCoh = 1:size(cts,2)
        
        % Build coherence index
        if iCoh < size(cts,2)
            cIdx                = ts >= cts(iCoh) & ts < cts(iCoh+1);
        else
            cIdx                = ts >= cts(iCoh) & ts <= trl.tEnd(iTrl);
        end
        
        if sum(cIdx) == 0
            continue
        end
        
        % Correlation analysis
        count                   = count + 1;
        coh(count)              = rdp_coh{iTrl}(iCoh);                                          % Coherence ID
        cc(count)               = circ_corrcc(deg2rad(js_vec(cIdx)),deg2rad(rdp_vec(cIdx)));  	% Circular correlation between stimulus and joystick direction
        xc(count,:)             = xcorr(double(abs(rdp_dff(cIdx))),abs(js_dff(cIdx)),nLag);    	% Cross-correlation between stimulus and joystick direction
        sxc(count,:)            = smoothdata(xc(count,:),'gaussian',20);                        % Smooth correlation output with gaussian kernel
        maxR(count)           	= max(sxc(count,1:nLag+1));                                     % Max cross-correlation coefficient
        posPk(count)        	= find(sxc(count,1:nLag+1) == maxR(count));                     % Peak position of cross-correlation

        try
            auPk(count,:)     	= trapz(sxc(count,posPk(count)-10:posPk(count)+10));           	% Area under cross-correlation peak
        catch
            auPk(count,:)     	= nan;
        end
        
        if plotFlag
            % Plot trial-wise cross-correlation
            ps              	= plot(sxc(count,1:nLag),'Color',[.5 .5 .5 .1], 'Marker','none', 'LineWidth',2);
        else 
            ps                  = [];
        end
    end
end

out.coh                         = coh;
out.cc                          = cc;
out.xc                          = xc;
out.sxc                         = sxc;
out.maxR                        = maxR;
out.posPk                       = posPk;
out.auPk                        = auPk;

end