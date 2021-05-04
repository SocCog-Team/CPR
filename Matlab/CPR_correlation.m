function [out, ps] = CPR_correlation(t,tindx, nLag, plotFlag)

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
js_dir_ts           = t.trl_js_dir_ts(tindx);
% js_str              = t.trl_js_str(tindx);                                % Joystick strength
% js_str_ts           = t.trl_js_str_ts(tindx);
rdp_coh            	= t.trl_rdp_coh(tindx);                               	% RDP coherence
rdp_coh_ts        	= t.trl_rdp_coh_ts(tindx);
rdp_coh{1}(1)   	= [];                                                   % Remove first entry of first trial
rdp_coh_ts{1}(1)  	= [];
count               = 0;

%% Loop through blocks and correlate stimulus and behaviour

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
    
    % Correct for circular space: Unwrap data
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
        
        % Correlation analysis
        count                   = count + 1;
        coh(count)              = rdp_coh{iTrl}(iCoh);                                          % Coherence ID
        cc(count)               = circ_corrcc(deg2rad(js_dir{iTrl}),deg2rad(dir_vec));          % Circular correlation between stimulus and joystick direction
        xc(count,:)             = xcorr(double(abs(rdp_dff(cIdx))),abs(js_dff(cIdx)),nLag);    	% Cross-correlation between stimulus and joystick direction
        maxR(count)           	= max(xc(iTrl,1:nLag+1));                                      	% Max cross-correlation coefficient
        posPk(count)        	= find(xc(count,1:nLag+1) == max(xc(count,1:nLag+1)));          % Peak position of cross-correlation

        try
            auPk(count,:)     	= trapz(xc(count,posPk(count)-10:posPk(count)+10));           	% Area under cross-correlation peak
        catch
            auPk(count,:)     	= nan;
        end
        
        if plotFlag
            % Plot trial-wise cross-correlation
            ps              	= plot(xc(count,1:nLag),'Color',[.5 .5 .5 .1], 'Marker','none', 'LineWidth',2);
        else 
            ps                  = [];
        end
    end
end

out.coh                         = coh;
out.cc                          = cc;
out.xc                          = xc;
out.maxR                        = maxR;
out.posPk                       = posPk;
out.auPk                        = auPk;

end