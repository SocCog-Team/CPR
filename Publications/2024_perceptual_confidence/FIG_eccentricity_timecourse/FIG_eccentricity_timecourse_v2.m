close all
clear all

% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/preprocessing/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/mat_to_summary/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

% Import subject summary spreadsheet
pth                         = '/Volumes/T7_Shield/CPR_psychophysics/';      % Local hard drive
x                           = readtable([pth 'Subjects_summary.xlsx']);     % Spreadsheet
sbj_lst                     = x.Abbreviation;                               % Subject ID list
sbj_lst(cellfun(@isempty,sbj_lst)) = [];

cc                          = 0;

% For all subjects
for iSubj = 1:length(sbj_lst)
    
    disp(['Processing subject: ' sbj_lst{iSubj}])
    data_pth                                = [pth sbj_lst{iSubj} '/raw/'];     % Data path
    
    if isfolder(data_pth)
        cd(data_pth)
        h5_files                            = dir('*.h5');                      % Get all .mat files in directory
        
        % For all files in directory
        for iFile = 1:length(h5_files)
            
            % Solo experiments
            if contains(h5_files(iFile).name,'CPRsolo')
                
                % Load file
                d                        	= MW_readH5(h5_files(iFile).name); % ...load .h5 file
                exp_info                    = strsplit(h5_files(iFile).name,'_');
                
                % Initialise variables
                idx                         = [];
                cyc                         = [];
                
                % Create variable-specific indices
                idx.cOn                     = d.event == 'TRIAL_start';
                idx.cEnd                    = d.event == 'TRIAL_end';
                idx.frame                   = d.event == 'STIM_displayUpdate';
                idx.RDP_onset               = d.event == 'STIM_RDP_onset';
                idx.RDP_dir                 = d.event == 'STIM_RDP_direction';
                idx.RDP_coh                 = d.event == 'STIM_RDP_coherence';
                idx.RDP_dot                	= d.event == 'STIM_RDP_dotPositions';
                idx.trg_on                  = d.event == 'STIM_target_onset';
                idx.JS_dir                  = d.event == 'IO_joystickDirection';
                idx.JS_str                  = d.event == 'IO_joystickStrength';
                idx.fixation              	= d.event == 'IO_fixation_flag';
                idx.outcome                 = d.event == 'TRIAL_outcome';
                idx.trg                     = d.event == 'TRIAL_reactionTrigger';
                idx.eye_x_dva              	= d.event == 'EYE_x_dva';
                idx.eye_y_dva             	= d.event == 'EYE_y_dva';
                
                % Get trial timestamps
                cyc.cOn                     = d.time(idx.cOn);
                cyc.cEnd                    = d.time(idx.cEnd);
                
                % Stimulus cycle (i.e. trial) loop
                for iCyc = 1:length(cyc.cEnd)
                    
                    % Cycle index
                    cycIdx                  = [];
                    cycIdx                  = d.time >= cyc.cOn(iCyc) & d.time <= cyc.cEnd(iCyc);
                    
                    % Extract nominal stimulus direction
                    tmp_dir_ts              = getTrialData(d.time, cycIdx, idx.RDP_dir);                    % RDP_direction timestamps
                    tmp_dir                 = getTrialData(d.value, cycIdx, idx.RDP_dir);                   % RDP_direction
                    tmp_coh                 = getTrialData(d.value, cycIdx, idx.RDP_coh);               	% RDP coherence
                    tmp_coh_ts              = getTrialData(d.time, cycIdx, idx.RDP_coh);                 	% RDP coherence timestamps
                    cyc.state_ts{iCyc}      = tmp_dir_ts(2:end);                                            % State ON duration timestamps
                    
                    % Stimulus state loop
                    for iSS = 1:length(cyc.state_ts{iCyc})
                        
                        disp(['Cycle: ' num2str(iCyc) ', State: ' num2str(iSS)])
                        
                        % Stimulus state index
                        if iSS < length(cyc.state_ts{iCyc})
                            ssIdx           = d.time >= cyc.state_ts{iCyc}(iSS) & d.time < cyc.state_ts{iCyc}(iSS+1);
                        else
                            ssIdx           = d.time >= cyc.state_ts{iCyc}(iSS) & d.time <= cyc.cEnd(iCyc);
                        end
                        
                        % Trial/State/Stimulus parameter
                        cc                  = cc+1;                                                         % Stimulus state counter
                        rdp_coh(cc)         = tmp_coh(find(tmp_coh_ts <= cyc.state_ts{iCyc}(iSS),1,'last'));% Stimulus state coherence
                        rdp_dir(cc)         = tmp_dir(find(tmp_dir_ts <= cyc.state_ts{iCyc}(iSS),1,'last'));% Stimulus direction
                        
                        %%%% TIME WINDOW ANALYSIS %%%
                        % Extract dot position
                        dp                  = d.value(ssIdx & idx.RDP_dot);
                        dp_ts               = d.time(ssIdx & idx.RDP_dot); % Frame timestamp
                        
                        % Calculate motion energy
                        frme_vec{cc}      	= extract_resultant_vector(dp,dp_ts,rdp_dir(cc),rdp_coh(cc));
                        
                        % Extract joystick samples
                        tmp_js_ts{cc}       = getTrialData(d.time, ssIdx, idx.JS_str);                      % Timestamps: Joystick strength
                        tmp_js_dir{cc}      = getTrialData(d.value, ssIdx, idx.JS_dir);                     % Joystick direction
                        tmp_js_ecc{cc}      = getTrialData(d.value, ssIdx, idx.JS_str);                     % Joystick strength
                        
                        % Build frame-wise vector for joystick data
                        for iFrme = 1:length(dp_ts)
                            fIdx                  	= [];
                            fIdx                 	= find(tmp_js_ts{cc} < dp_ts(iFrme),1,'last');  % Extract last entry before frame onset
                            if sum(fIdx) == 0 || isempty(fIdx) || fIdx > length(tmp_js_dir{cc}) || fIdx > length(tmp_js_ecc{cc})
                                js_dir{cc}(iFrme)	= nan;                                                  % Write to vector
                                js_ecc{cc}(iFrme)	= nan;
                            else
                                js_dir{cc}(iFrme)	= tmp_js_dir{cc}(fIdx);
                                js_ecc{cc}(iFrme)	= tmp_js_ecc{cc}(fIdx);
                            end
                        end
                    end
                end
            end
        end
    end
    
    save(['/Users/fschneider/Desktop/vector_subk_'  sbj_lst{iSubj}], 'frme_vec', '-v7.3')
    save(['/Users/fschneider/Desktop/joystick_subk_' sbj_lst{iSubj}], 'js_dir','js_ecc', '-v7.3')
end
                        
%%
clear avg_ecc avg_acc

win = 60;

for iSubj = 1:length(sbj_lst)
    
    mat_acc = nan(10000,300);
    mat_ecc = nan(10000,300);
    c = 0;
    
    load(['/Users/fschneider/Desktop/resultant_vec/vector_subj_' sbj_lst{iSubj} '.mat'])
    load(['/Users/fschneider/Desktop/resultant_vec/joystick_subj_' sbj_lst{iSubj} '.mat'])
    
    clear ncoh slen ndir
    
    for iState = 1:size(frme_vec,2)
        
        if length(js_ecc{iState}) >= win && length(js_ecc{iState}) == length(frme_vec{iState}.resultant_ang_error)+1
            c           = c+1;
            slen(c)     = length(js_ecc{iState}(win:end));
            ncoh(c)    	= frme_vec{iState}.nominal_coh;
            ndir(c)    	= frme_vec{iState}.nominal_dir_deg;
            
            mat_acc(c,end-slen(c)+1:end)  	= abs(1 - abs(frme_vec{iState}.resultant_ang_error(win-1:end)) / 180);
            mat_ecc(c,end-slen(c)+1:end) 	= js_ecc{iState}(win:end);
        end
    end
    
    snr_lst = unique(ncoh);
    if length(snr_lst)>7
       snr_lst = snr_lst(1:2:end);
    end
    
    mat_acc(sum(isnan(mat_acc),2) == size(mat_acc,2),:) = [];
    mat_ecc(sum(isnan(mat_ecc),2) == size(mat_ecc,2),:) = [];
    
    for iCoh = 1:length(snr_lst)
        cidx = ncoh == snr_lst(iCoh);
        avg_acc(iSubj,iCoh,:) = nanmean(mat_acc(cidx,:));
        avg_ecc(iSubj,iCoh,:) = nanmean(mat_ecc(cidx,:));
    end
end

%%% PLOT
f = figure;
col = cool(length(snr_lst));
ax1 = subplot(1,2,1); hold on
ax2 = subplot(1,2,2); hold on

for iCoh = 1:length(snr_lst)
   
    [ax1, pl(iCoh)] = plot_averages(ax1,squeeze(avg_acc(:,iCoh,:)),col(iCoh, :));
    [axw, pl(iCoh)] = plot_averages(ax2,squeeze(avg_ecc(:,iCoh,:)),col(iCoh, :));

end

ax1 = adjust_axes(ax1,'accuracy [norm]');
ax2 = adjust_axes(ax2,'tilt [norm]');

print(f, '/Users/fschneider/Desktop/acg_timeline_pop', '-r500', '-dpng');


%% Crosscorrelation between vector length and joystick displacement

addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/

nLag        = 150;
smooth_win  = 20;
n_ofs       = 5;
nshift      = 35;

clear avg_xc_acc_vec avg_xc_ecc_vec avg_xc_ecc_acc avg_xc_control

for iSubj = 1:10%length(sbj_lst) 
    
    clear frme_vec js_ecc js_dir xc_acc_vec xc_ecc_vec xc_ecc_acc state_coh xc_control 

    load(['/Users/fschneider/Desktop/resultant_vec/vector_subj_' sbj_lst{iSubj} '.mat'])
    load(['/Users/fschneider/Desktop/resultant_vec/joystick_subj_' sbj_lst{iSubj} '.mat'])
    cnt = 0;

    for iState = 1:size(frme_vec,2)
        clear v_smooth v_smooth_detrend ecc_detrend js_dir_error joystick_acc acc_detrend
        
        if length(js_ecc{iState}) == length(frme_vec{iState}.resultant_length)+1
            cnt = cnt+1;
            v_smooth = smoothdata(frme_vec{iState}.resultant_length,'gaussian',smooth_win);
            v_smooth_detrend = [0 v_smooth - nanmean(v_smooth)];
            
            ecc_detrend = js_ecc{iState} - nanmean(js_ecc{iState});
            
            js_dir_error = rad2deg(circ_dist(deg2rad(js_dir{iState}), deg2rad(frme_vec{iState}.nominal_dir_deg)));
            joystick_acc = abs(1 - abs(js_dir_error) / 180);
            acc_detrend = joystick_acc - nanmean(joystick_acc);
            
            
            state_coh(cnt) = frme_vec{iState}.nominal_coh;
            
            [xc_acc_vec(cnt,:) lags] = xcorr(ecc_detrend(n_ofs:end),v_smooth_detrend(n_ofs:end),nLag,'normalized');
            [xc_ecc_vec(cnt,:) lags] = xcorr(acc_detrend(n_ofs:end),v_smooth_detrend(n_ofs:end),nLag,'normalized');
            [xc_ecc_acc(cnt,:) lags] = xcorr(ecc_detrend(n_ofs:end),acc_detrend(n_ofs:end),nLag,'normalized');
            [xc_control(cnt,:) lags] = xcorr(ecc_detrend(n_ofs:end-nshift+n_ofs),ecc_detrend(nshift:end),nLag,'normalized');
        end
    end
    
%     snr_lst = unique(state_coh);
%     if length(snr_lst)>7
%         snr_lst = snr_lst(1:2:end);
%     end
    
    for iCoh = 1:length(snr_lst)
        cidx = state_coh == snr_lst(iCoh);
        str{iCoh} = num2str(round(snr_lst(iCoh),2));
        avg_xc_acc_vec(iSubj,iCoh,:) = nanmean(xc_acc_vec(cidx,:));
        avg_xc_ecc_vec(iSubj,iCoh,:) = nanmean(xc_ecc_vec(cidx,:));
        avg_xc_ecc_acc(iSubj,iCoh,:) = nanmean(xc_ecc_acc(cidx,:));
        avg_xc_control(iSubj,iCoh,:) = nanmean(xc_control(cidx,:));
    end
end

%% PLOT
f = figure('units','centimeters','position',[0 0 50 10]);

for iPlot = 1:4
    clear dat
    if iPlot == 1
        dat = avg_xc_control;
        ylim = [-.5 1];
        tstr = 'control - shifted by 250ms';
    elseif iPlot == 2
        dat = avg_xc_ecc_acc;
        tstr = 'accuracy - ecc';
        ylim = [-.1 .35];
    elseif iPlot == 3
        dat = avg_xc_ecc_vec;
        tstr = 'ecc - res.vec';
        ylim = [-.1 .1];
    elseif iPlot == 4
        dat = avg_xc_acc_vec;
        tstr = 'acc - res.vec';
        ylim = [-.1 .1];
    end
    
    ax = subplot(1,4,iPlot); hold on
    ln = line([nLag nLag],[-1 1], 'Color', 'k','LineStyle', ':', 'LineWidth', 2);    

    for iCoh = 1:length(snr_lst)
%         plt(iCoh) = plot(lags./120,mean(squeeze(dat(:,iCoh,:))),'Color', col(iCoh,:), 'LineWidth', 2);
         [ax, pl] = xc_plot_averages(ax,squeeze(dat(:,iCoh,:)),col(iCoh,:));
    end
    
    ax.XLabel.String = 'Shift [ms]';
    ax.YLabel.String = 'XCorr coef';
    ax.Title.String = tstr;
    ax.FontSize = 14;
    ax.YLim = ylim;
    ax.XTick = [0 75 150 225 300];
    
    for iLab = 1:length(ax.XTickLabel)
        ax.XTickLabel{iLab} = num2str( round((str2num(ax.XTickLabel{iLab})-nLag) * (1000/120)));
    end
    
    axis tight
end

lg = legend(plt,str,'Location','West');

% print(f, '/Users/fschneider/Desktop/xcorr', '-r500', '-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = extract_resultant_vector(dp, dp_ts, nominal_dir, nominal_coh)

% Initialize variables
good_frame = 1;
c = 0;

% Loop through each frame of the dot positions
for iFrame = 2:length(dp)
    % Skip frame if it contains NaN values
    if any(isnan(dp{iFrame}))
        continue
    end
    
    % Increment frame count
    c = c + 1;
    
    % Separate x and y positions
    xIndices = logical(mod(1:size(dp{iFrame}, 1), 2));
    xPos = dp{iFrame}(xIndices);
    yPos = dp{iFrame}(~xIndices);
    
    % Get positions from the last valid frame
    xPosLast = dp{good_frame}(xIndices);
    yPosLast = dp{good_frame}(~xIndices);
    
    % Record timestamp of the current frame
    dts(c) = dp_ts(iFrame);
    
    % Calculate vector properties for each dot
    for iDot = 1:length(yPos)
        % Vector start and end positions
        start_pos = [xPosLast(iDot), yPosLast(iDot)];
        end_pos = [xPos(iDot), yPos(iDot)];
        
        % Distance between start and end positions
        dist(iDot) = pdist([start_pos; end_pos], 'euclidean');
        
        % Calculate the direction and delta between points
        delta = end_pos - start_pos;
        dot_dir{c}(iDot) = mod(atan2d(delta(1), delta(2)), 360);
        dt{c}(iDot, :) = delta;
    end
    
    % Identify indices of dots with median distance (filter for stable points)
    dot_idx{c} = dist >= median(dist) - 0.0001 & dist <= median(dist) + 0.0001;
    
    % Update the last valid frame
    good_frame = iFrame;
end

% Initialize output arrays
ofs = 1;
clear resultant res_length res_ts rcoh

% Loop through each calculated vector set to find the resultant
for iVec = 1:length(dt)
    % Calculate mean difference vector for stable dots
    meanDelta = mean(dt{iVec}(dot_idx{iVec}, :));
    
    % Compute resultant vector direction and length
    res_deg(iVec) = mod(atan2d(meanDelta(1), meanDelta(2)), 360);
    res_length(iVec) = pdist([[0, 0]; meanDelta], 'euclidean') / median(dist(dot_idx{iVec}));
    res_ts(iVec) = double(dts(iVec) - dts(1)) / 1e6;
    
    % Calculate angular error with nominal direction
    ang_error(iVec) = rad2deg(circ_dist(deg2rad(res_deg(iVec)), deg2rad(nominal_dir)));
    
    % Check coherence of dot directions relative to nominal direction
    dot_angles = dot_dir{iVec}(dot_idx{iVec});
    signal = abs(rad2deg(circ_dist(deg2rad(dot_angles), deg2rad(nominal_dir)))) < ofs;
    rcoh(iVec) = sum(signal) / length(signal);
end

% Compile results into output structure
out.nominal_dir_deg = nominal_dir;
out.nominal_coh = nominal_coh;
out.resultant_ang_error = ang_error;
out.resultant_deg = res_deg;
out.resultant_length = res_length;
out.resultant_ts = res_ts;
out.actual_coherence = rcoh;

end

function [ax, pl] = plot_averages(ax,dat,col)
    axes(ax)
    
    % Boostrap confidence intervals
    nRep        	= 5000;
    ofs             = 100;
    [CI,~]        	= bootci(nRep,{@mean,dat(:,100:end)},'Alpha',0.05);
    
    % Prepare filled area
    vec           	= 1:length(CI);
    x_spacing     	= [vec fliplr(vec)]+ofs;
    ci            	= [CI(1,:) fliplr(CI(2,:))];

    % Overlay confidence intervals
    fl           	= fill(x_spacing,ci,col,'EdgeColor','none', 'FaceAlpha', .3);
    
    % Plot mean curve
    pl              = plot(nanmean(dat), 'LineWidth', 2, 'Color', col); 
    
end

function ax = adjust_axes(ax,ystr)
    ax.YLabel.String 	= ystr;
    ax.XLabel.String 	= 'time [ms]';
    ax.YLim             = [0 1];
    ax.XLim             = [100 300];
    ax.XTick            = 100:100:300;
    tmp_lab             = cellfun(@str2num, ax.XTickLabel, 'UniformOutput', false);
    ax.FontSize         = 16;
    
    for iLab = 1:length(tmp_lab)
        ax.XTickLabel{iLab} = num2str(round((tmp_lab{iLab} - ax.XTick(end)) * 8.333));
    end
end

function [ax, pl] = xc_plot_averages(ax,dat,col)
    axes(ax); hold on
    
    % Boostrap confidence intervals
    nRep        	= 5000;
    [CI,~]        	= bootci(nRep,{@mean,dat},'Alpha',0.05);
    
    % Prepare filled area
    vec           	= 1:length(CI);
    x_spacing     	= [vec fliplr(vec)];
    ci            	= [CI(1,:) fliplr(CI(2,:))];

    % Overlay confidence intervals
    fl           	= fill(x_spacing,ci,col,'EdgeColor','none', 'FaceAlpha', .3);
    
    % Plot mean curve
    pl              = plot(nanmean(dat), 'LineWidth', 2, 'Color', col); 
    
end