function [out, ps] = CPR_correlation_analysis(in, nLag, plotFlag)

% This function extracts stimulus and joystick data in a specified time
% window at the end of a steady state
%
% Input:  	.in             Table, Contains stimulus state information
%          	.nLag           Integer, Number of samples before direction
%                           changes (End of state)
%           .plotFlag       Boolean, Plot: yes - no
%
% Output:   .out            Structure, Contains averages and performance
%                           measures for each coherence level
%
% Example:	out = CPR_time_window(tbl,29,true)
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

if nargin < 1
    error('Input missing')
end

%% Loop through blocks and correlate stimulus and behaviour

count                           = 0;

for iTrl = 1:size(in.frme_ts,2)
    
    % Sample-by-sample difference
    clear js_dff rdp_dff
    
    excl                        = isnan(in.rdp_coh{iTrl}) | isnan(in.rdp_dir{iTrl}) | isnan(in.js_dir{iTrl}) | isnan(in.frme_ts{iTrl});
    if ~isempty(excl)
        in.rdp_coh{iTrl}(excl)	= [];
        in.rdp_dir{iTrl}(excl) 	= [];
        in.js_dir{iTrl}(excl)  	= [];
        in.frme_ts{iTrl}(excl) 	= [];
    end
    
    rfr{iTrl}                   = median(in.refresh{iTrl});
    
    if isempty(in.frme_ts{iTrl})
        continue
    end
    
    for iSample = 1:size(in.frme_ts{iTrl},2)-1
        js_dff(iSample)     	= rad2deg(circ_dist(deg2rad(in.js_dir{iTrl}(iSample)),deg2rad(in.js_dir{iTrl}(iSample+1))));
        rdp_dff(iSample)     	= rad2deg(circ_dist(deg2rad(in.rdp_dir{iTrl}(iSample)),deg2rad(in.rdp_dir{iTrl}(iSample+1))));
    end
 
    js_dff                      = [0 js_dff];
    rdp_dff                     = [0 rdp_dff];
    
    % NANs coming from table-building or MWorks?!
    ex                          = isnan(js_dff) | isnan(rdp_dff);
    js_dff(ex)                  = [];
    rdp_dff(ex)                 = [];
    
    % Extract coherence chunks
    cindx                       = diff(in.rdp_coh{iTrl}) ~=0;
    cid                         = in.rdp_coh{iTrl}([true cindx]);
    cts                         = in.frme_ts{iTrl}([true cindx]);
            
    if length(cid) > 3 || sum(rdp_dff) == 0
        continue
    end
        
    for iCoh = 1:size(cts,2)
                
        % Build coherence index
        if iCoh < size(cts,2)
            cIdx                = in.frme_ts{iTrl} >= cts(iCoh) & in.frme_ts{iTrl} < cts(iCoh+1);
        else
            cIdx                = in.frme_ts{iTrl} >= cts(iCoh) & in.frme_ts{iTrl} <= in.frme_ts{iTrl}(end);
        end
        
        if sum(cIdx) < 10 %|| isempty(find(diff(rdp_dir{iTrl}(cIdx))))
            continue
        end
        
        % Correlation analysis
        count                   = count + 1;
        coh(count)              = cid(iCoh);                                                                % Coherence ID
        cc(count)               = circ_corrcc(deg2rad(in.js_dir{iTrl}(cIdx)),deg2rad(in.rdp_dir{iTrl}(cIdx)));	% Circular correlation between stimulus and joystick direction
        xc(count,:)             = xcorr(abs(js_dff(cIdx)),double(abs(rdp_dff(cIdx))),nLag,'normalized');   	% Cross-correlation between stimulus and joystick direction
        sxc(count,:)            = smoothdata(xc(count,:),'gaussian',20);                                    % Smooth correlation output with gaussian kernel
        maxR(count)           	= max(sxc(count,nLag:end));                                                 % Max cross-correlation coefficient
        posPk(count)        	= find(sxc(count,nLag:end) == maxR(count));                                 % Peak position of cross-correlation
        
        fs(count)               = rfr{iTrl};
        lg(count)             	= (find(sxc(count,nLag:end) == max(sxc(count,nLag:end))) * fs(count)) / 1e3; % Calculate lag using sample rate

        try
            auPk(count,:)     	= trapz(sxc(count,posPk(count)-10:posPk(count)+10));                        % Area under cross-correlation peak
        catch
            auPk(count,:)     	= nan;
        end
        
        if plotFlag
            % Plot trial-wise cross-correlation
            ps              	= plot(sxc(count,nLag:end),'Color',[.5 .5 .5 .1], 'Marker','none', 'LineWidth',2);
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
out.lag                       	= lg;

end