function [out, ps] = CPR_correlation_analysis(t, nLag, plotFlag)

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
%   1.2     (fxs 2021-11-25) Simplify code with frame-wise analysis

%% Check input

addpath /Users/fschneider/Documents/MATLAB/CircStat2012a

if nargin < 3 || isempty(plotFlag)
    plotFlag     	= true;
end

if nargin < 2 || isempty(nLag)
    nLag            = 150;
end

if nargin < 1 || ~istable(t)
    error('Input must be a table')
end

%% Assign data from table

rdp_dir             = t.trl_rdp_dir;                              	% RDP direction
js_dir              = t.trl_js_dir;                               	% Joystick direction
% js_str              = t.trl_js_str;                                % Joystick strength
rdp_coh            	= t.trl_rdp_coh;                               	% RDP coherence
ts                  = t.trl_frme_ts;
count               = 0;

%% Loop through blocks and correlate stimulus and behaviour

for iTrl = 1:size(ts,1)
    
    % Sample-by-sample difference
    clear js_dff rdp_dff
    
    excl                        = isnan(rdp_coh{iTrl}) | isnan(rdp_dir{iTrl}) | isnan(js_dir{iTrl}) | isnan(ts{iTrl});
    rdp_coh{iTrl}(excl)         = [];
    rdp_dir{iTrl}(excl)        	= [];
    js_dir{iTrl}(excl)         	= [];
    ts{iTrl}(excl)              = [];
    
    if isempty(ts{iTrl})
        continue
    end
    
    for iSample = 1:size(ts{iTrl},2)-1
        js_dff(iSample)     	= rad2deg(circ_dist(deg2rad(js_dir{iTrl}(iSample)),deg2rad(js_dir{iTrl}(iSample+1))));
        rdp_dff(iSample)     	= rad2deg(circ_dist(deg2rad(rdp_dir{iTrl}(iSample)),deg2rad(rdp_dir{iTrl}(iSample+1))));
    end
 
    js_dff                      = [0 js_dff];
    rdp_dff                     = [0 rdp_dff];
    
    % NANs coming from table-building or MWorks?!
    ex                          = isnan(js_dff) | isnan(rdp_dff);
    js_dff(ex)                  = [];
    rdp_dff(ex)                 = [];
    
    % Extract coherence chunks
    cindx                       = diff(rdp_coh{iTrl}) ~=0;
    cid                         = rdp_coh{iTrl}([true cindx]);
    cts                         = ts{iTrl}([true cindx]);
            
    if sum(cid) > 3 
        continue
    end
        
    for iCoh = 1:size(cts,2)
                
        % Build coherence index
        if iCoh < size(cts,2)
            cIdx                = ts{iTrl} >= cts(iCoh) & ts{iTrl} < cts(iCoh+1);
        else
            cIdx                = ts{iTrl} >= cts(iCoh) & ts{iTrl} <= ts{iTrl}(end);
        end
        
        if sum(cIdx) < 10 %|| isempty(find(diff(rdp_dir{iTrl}(cIdx))))
            continue
        end
        
        % Correlation analysis
        count                   = count + 1;
        coh(count)              = cid(iCoh);                                                                % Coherence ID
        cc(count)               = circ_corrcc(deg2rad(js_dir{iTrl}(cIdx)),deg2rad(rdp_dir{iTrl}(cIdx)));	% Circular correlation between stimulus and joystick direction
        xc(count,:)             = xcorr(double(abs(rdp_dff(cIdx))),abs(js_dff(cIdx)),nLag);                 % Cross-correlation between stimulus and joystick direction
        sxc(count,:)            = smoothdata(xc(count,:),'gaussian',20);                                    % Smooth correlation output with gaussian kernel
        maxR(count)           	= max(sxc(count,1:nLag+1));                                                 % Max cross-correlation coefficient
        posPk(count)        	= find(sxc(count,1:nLag+1) == maxR(count));                                 % Peak position of cross-correlation

        try
            auPk(count,:)     	= trapz(sxc(count,posPk(count)-10:posPk(count)+10));                        % Area under cross-correlation peak
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