close all
clear all

% Add relevant directories
% addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
% addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/

% Load file
d                        	= MW_readH5('D:\Temp/20221206_aaa_CPRsolo_block2_psycho3_fxs.h5'); % ...load .h5 file

% Initialise variables
idx                         = [];
cyc                         = [];
cc                          = 0;

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
for iCyc = 1%:length(cyc.cEnd)
    
    % Cycle index
    cycIdx                  = [];
    cycIdx                  = d.time >= cyc.cOn(iCyc) & d.time <= cyc.cEnd(iCyc);
    
    % Extract nominal stimulus direction
    tmp_dir_ts              = getTrialData(d.time, cycIdx, idx.RDP_dir);                    % RDP_direction timestamps
    tmp_dir                 = getTrialData(d.value, cycIdx, idx.RDP_dir);                   % RDP_direction
    tmp_coh                 = getTrialData(d.value, cycIdx, idx.RDP_coh);               	% RDP coherence
    tmp_coh_ts              = getTrialData(d.time, cycIdx, idx.RDP_coh);                 	% RDP coherence timestamps
    cyc.state_ts{iCyc}      = tmp_dir_ts(2:end);                                            % State ON duration timestamps
    cyc.dir{iCyc}           = tmp_dir(2:end);                                               % Cycle RDP direction
    
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
        cc                  = cc+1;                                                         % Steady state counter
        rdp_coh(cc)         = tmp_coh(find(tmp_coh_ts <= cyc.state_ts{iCyc}(iSS),1,'last'));% Stimulus state coherence
        rdp_dir(cc)         = mod(cyc.dir{iCyc}(iSS),360);                                  % Stimulus direction
        
        % Extract dot position
        dp                  = d.value(ssIdx & idx.RDP_dot);
        dp_ts               = d.time(ssIdx & idx.RDP_dot); % Frame timestamp
        
        % Calculate motion energy
        frme_vec{cc}      	= extract_resultant_vector(dp,dp_ts,rdp_dir(cc));
        vec_len{cc}      	= frme_vec{cc}.resultant_length;
        act_coh{cc}      	= frme_vec{cc}.actual_coherence;
        
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

%% scatter: avg ecc vs avg vector length

nSample         = 99;                                           % Time window size [samples]

indx        	= ~(cell2mat(cellfun(@length,js_ecc', 'UniformOutput', false)) <= nSample+1);
roi_ecc         = cell2mat(cellfun(@(x) x(end-nSample:end),js_ecc(indx)', 'UniformOutput', false));
roi_vec         = cell2mat(cellfun(@(x) x(end-nSample:end),vec_len(indx)', 'UniformOutput', false));
roi_coh         = cell2mat(cellfun(@(x) x(end-nSample:end),act_coh(indx)', 'UniformOutput', false));

e               = mean(roi_ecc,2);
v              	= mean(roi_vec,2);
c               = mean(roi_coh,2);

figure;hold on
for i = 1:length(e)
    sc          = scatter(e(i),v(i) ,'filled', 'MarkerFaceAlpha',.3, 'MarkerFaceColor', [c(i) 0 0]);
end

xlabel('avg eccentricity')
ylabel('avg vector length')
set(gca, 'fontsize', 14)

%% Cycle plot


figure;
plot(roi_vec');

idx = 1:iSS;
f = figure('units','normalized','position',[0 0 1 1]); hold on
subplot(2,3,1)
imagesc(roi_ecc(idx,:))
xlabel('# frames before direction change')
ylabel('# stimulus state')
title('joystick eccentricity')
set(gca, 'fontsize', 14)
set(gca, 'xtick', [1 50 100])
set(gca, 'xticklabels', [-100 -50 0])
colormap('gray')
c = colorbar;
c.Label.String = 'joystick eccentricity';

subplot(2,3,2)
imagesc(roi_vec(idx,:))
xlabel('# frames before direction change')
ylabel('# stimulus state')
title('resultant vector length')
set(gca, 'fontsize', 14)
set(gca, 'xtick', [1 50 100])
set(gca, 'xticklabels', [-100 -50 0])
colormap('gray')
c = colorbar;
c.Label.String = 'vector length';

subplot(2,3,3)
imagesc(roi_coh(idx,:))
xlabel('# frames before direction change')
ylabel('# stimulus state')
title('physical coherence')
set(gca, 'xtick', [1 50 100])
set(gca, 'xticklabels', [-100 -50 0])
set(gca, 'fontsize', 14)
colormap('gray')
c = colorbar;
c.Label.String = 'Dot coherence';

subplot(2,3,4)
imagesc(diff([nan(size(roi_ecc(idx,:),1),1) roi_ecc(idx,:)],1,2))
xlabel('# frames before direction change')
ylabel('# stimulus state')
set(gca, 'xtick', [1 50 100])
set(gca, 'xticklabels', [-100 -50 0])
set(gca, 'fontsize', 14)
colormap('gray')
c = colorbar;
c.Label.String = 'Derivative';

subplot(2,3,5)
imagesc(diff([nan(size(roi_vec(idx,:),1),1) roi_vec(idx,:)],1,2))
xlabel('# frames before direction change')
ylabel('# stimulus state')
set(gca, 'xtick', [1 50 100])
set(gca, 'xticklabels', [-100 -50 0])
set(gca, 'fontsize', 14)
colormap('gray')
c = colorbar;
c.Label.String = 'Derivative';

subplot(2,3,6)
imagesc(diff([nan(size(roi_coh(idx,:),1),1) roi_coh(idx,:)],1,2))
xlabel('# frames before direction change')
ylabel('# stimulus state')
set(gca, 'xtick', [1 50 100])
set(gca, 'xticklabels', [-100 -50 0])
set(gca, 'fontsize', 16)
colormap('gray')
c = colorbar;
c.Label.String = 'Derivative';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = extract_resultant_vector(dp, dp_ts, nominal_dir)

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
clear res_deg res_length res_ts rcoh ang_error

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
out.resultant_ang_error = ang_error;
out.resultant_deg = res_deg;
out.resultant_length = res_length;
out.resultant_ts = res_ts;
out.actual_coherence = rcoh;

end