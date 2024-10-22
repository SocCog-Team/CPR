addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/mat_to_summary/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/preprocessing

file_pth                    = '/Users/fschneider/Desktop/social_context_study/';
dest_pth                    = '/Users/fschneider/Desktop/social_context_study/stim/';
lb_fs                       = 14;
lw                          = 1.5;

for iDyad = 10
    % Extract all files
    h5Files = dir(fullfile([file_pth 'Dyad' num2str(iDyad) '/h5/'], '*.h5'));
    
    for iFile = 1%:length(h5Files)
        d                       = MW_readH5([file_pth '/Dyad' num2str(iDyad) '/h5/' h5Files(iFile).name]); % ...load .h5 file
        
        idx.cOn               	= d.event == 'TRIAL_start';
        idx.cEnd             	= d.event == 'TRIAL_end';
        idx.frame             	= d.event == 'STIM_displayUpdate';
        idx.RDP_dir           	= d.event == 'STIM_RDP_direction';
        idx.RDP_coh          	= d.event == 'STIM_RDP_coherence';
        idx.RDP_dot       		= d.event == 'STIM_RDP_dotPositions';
        idx.JS_dir              = d.event == 'IO_joystickDirection';
        idx.JS_str              = d.event == 'IO_joystickStrength';
        idx.JS2_dir          	= d.event == 'IO_joystickDirection2';
        idx.JS2_str         	= d.event == 'IO_joystickStrength2';
        idx.trg_on           	= d.event == 'STIM_target_onset';
        
        % Get cycle timestamps
        cyc.cOn              	= d.time(idx.cOn);
        cyc.cEnd              	= d.time(idx.cEnd);
        
        for iCyc = 1:length(cyc.cEnd)-1
            % Trial index
            cycIdx          	= [];
            cycIdx          	= d.time >= cyc.cOn(iCyc) & d.time <= cyc.cEnd(iCyc);
            c                   = 0;
            
            % Extract frame times
            dp              	= d.value(idx.RDP_dot & cycIdx);
            dp_ts              	= d.time(idx.RDP_dot & cycIdx);
            rdp_dir             = d.value(idx.RDP_dir);
            rdp_ts              = d.time(idx.RDP_dir);
            rdp_coh             = d.value(idx.RDP_coh);
            rdp_coh_ts          = d.time(idx.RDP_coh);
            dir_p1              = cell2mat(d.value(idx.JS_dir & cycIdx));
            dir_p2              = cell2mat(d.value(idx.JS2_dir & cycIdx));
            js_ts               = double(d.time(idx.JS2_dir & cycIdx))./1e6;

            n                   = length(dp);
            good_frme           = 1;
            clear dts scoh sdir

            for iFrme = 2:n
                disp(['Frame: ' num2str(iFrme)])
                
                if sum(isnan(dp{iFrme}))
                    continue
                end
                                
                % Dot position
                c                   = c+1;
                xIdx                = logical(mod([1:size(dp{iFrme},1)],2));
                xpos                = dp{iFrme}(xIdx);
                ypos                = dp{iFrme}(~xIdx);
                xpos_last           = dp{good_frme}(xIdx);      	% Last frame: x-position
                ypos_last         	= dp{good_frme}(~xIdx);      	% Last frame: y-position
                dts(c)              = dp_ts(iFrme);
                scoh(c)             = cell2mat(rdp_coh(find(dp_ts(iFrme) > rdp_coh_ts, 1,'last')));
                sdir(c)             = mod(cell2mat(rdp_dir(find(dp_ts(iFrme) > rdp_ts, 1,'last'))),360);
                
                for iDot = 2:length(ypos)
                    vs              = [xpos_last(iDot), ypos_last(iDot)];    	% dot position last frame
                    ve              = [xpos(iDot), ypos(iDot)];                 % dot position this frame
                    dist(iDot)      = pdist([vs;ve],'euclidean');               % distance/vector length between points
                    
                    delta               = ve - vs;                              % vector that points from x1 to x2
                    dot_dir{c}(iDot)	= mod(atan2d(delta(1),delta(2)),360);   % Alternative: rad2deg(cart2pol(delta(2), delta(1)))
                    dt{c}(iDot,:)       = delta;
                end
                
                dot_idx{c}              = dist >= median(dist)-.0001 & dist <= median(dist)+.0001; % Exclude new dots appearances
                good_frme               = iFrme;
            end
            
            ofs = 1;
            clear resultant res_length res_ts rcoh
            for i = 2:n-1
                % Avg frame direction
                mdf                     = mean(dt{i}(dot_idx{i},:));     	% Mean x/y of all dots that didn't jump
                resultant(i)            = mod(atan2d(mdf(1),mdf(2)),360); 	% Resultant vector [deg]
                res_length(i)           = pdist([[0 0];mdf],'euclidean');
                res_ts(i)               = double(dts(i)-dts(1))/1e6;
                
                % Check coherence level based on dot direction
                dots = dot_dir{i}(dot_idx{i});
                signal = abs(rad2deg(circ_dist(deg2rad(dots),deg2rad(sdir(i))))) < ofs; % Number of dots moving in nominal direction
                rcoh(i) = sum(signal)/length(signal);
            end
            
            %%% Polar histogram of individual cycle %%%
            f                           = figure('units','centimeters','position',[0 0 10 10]);
            ax                          = polaraxes;
            pp                          = polarhistogram(deg2rad(resultant), 180);
            pp.FaceColor                = [.2 .2 .2];
            pp.EdgeColor                = 'none';
            ax.ThetaZeroLocation        = 'top';
            ax.ThetaDir                 = 'clockwise';
            ax.FontSize                 = lb_fs;
            print(f, [dest_pth '/cyc' num2str(iCyc) '_rdp_dir_rvec_hist'], '-r500', '-dpng');
            
            %%% Derivative of individual cycle %%%
            f                           = figure('units','centimeters','position',[0 0 10 10]);
            pp                          = histogram(diff(resultant), 150);
            pp.FaceColor                = [.2 .2 .2];
            pp.EdgeColor                = 'none';
            ax                          = gca;
            ax.XTick                    = [-50 -25 0 25 50];
            ax.XLim                     = [-50 50];
            ax.YLabel.String            = 'Samples [#]';
            ax.XLabel.String            = 'delta [deg]';
            ax.FontSize                 = lb_fs;
            print(f, [dest_pth '/cyc' num2str(iCyc) '_deriv_hist'], '-r500', '-dpng');
            
            %%% Line plot of actual direction in individual cycle for stimulus & human response %%%
            f                           = figure('units','centimeters','position',[0 0 20 10]);
            ax1                         = subplot(2,1,1); hold on
            ps                          = plot(res_ts, resultant, 'Color', [.2 .2 .2], 'LineWidth', lw);
            p1                          = plot(js_ts-js_ts(1), dir_p1, 'Color', [0.3711 0.1328 0.8555], 'LineWidth', lw);
            p1                          = plot(js_ts-js_ts(1), dir_p2, 'Color', [0 .7 0], 'LineWidth', lw);
            ax1.YTick                   = [0 180 360];
            ax1.YLim                    = [0 360];
            ax1.YLabel.String           = 'Direction [deg]';
            ax1.XLabel.String           = 'Time [s]';
            ax1.FontSize                = lb_fs;
            ax1.Box                     = 'off';
%             lg                          = legend('RDP','P1','P2');
%             lg.FontSize                 = 10;
%             lg.Location                 = 'eastoutside';
            
            %%% Nominal and actual coherence %%%
            ax2                         = subplot(2,1,2); hold on   
            ps                          = plot(double(dts-dts(1))./1e6, scoh, 'Color', [.2 .2 .2], 'LineWidth', lw);
            sc                          = scatter(res_ts, rcoh, 'filled');
            sc.MarkerFaceAlpha          = .1;
            sc.MarkerFaceColor          = [.6 .2 .2];
            sc.MarkerEdgeColor          = 'none';
            ax2.YLim                    = [0 1];
            ax2.YLabel.String           = 'Coherence [%]';
            ax2.XLabel.String           = 'Time [s]';
            ax2.FontSize                = lb_fs;
            ax2.Box                     = 'off';
%             lg                          = legend('Nominal','Actual');
%             lg.FontSize                 = 10;
%             lg.Location                 = 'eastoutside';
            print(f, [dest_pth '/cyc' num2str(iCyc) '_timecourse'], '-r500', '-dpng');
        end
    end
end

