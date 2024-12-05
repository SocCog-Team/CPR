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
for iSubj = 23%:length(sbj_lst)
    
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
            end
        end
    end
end

save('/Users/fschneider/Desktop/vec_stuff', 'frme_vec', '-v7.3')
                        
%%
win = 120;
nSample = 29;
mat_acc = nan(10000,300);
mat_ecc = nan(10000,300);
c = 0;
snr = unique(ncoh);
col = cool(length(snr));

for iState = 1:size(frme_vec,2)

    if length(js_ecc{iState}) >= win
        c           = c+1;
        slen(c)     = length(js_ecc{iState}(win:end));
        ncoh(c)      = frme_vec{iState}.nominal_coh;
        ndir(c)      = frme_vec{iState}.nominal_dir_deg;
         
        mat_acc(c,end-slen(c)+1:end)  	= abs(1 - abs(frme_vec{iState}.resultant_ang_error(win-1:end)) / 180);
        mat_ecc(c,end-slen(c)+1:end) 	= js_ecc{iState}(win:end);
        
    end
end

mat_acc(sum(isnan(mat_acc),2) == size(mat_acc,2),:) = [];
mat_ecc(sum(isnan(mat_ecc),2) == size(mat_ecc,2),:) = [];

f = figure;
[v,i] = sort(slen);
subplot(2,2,1);hold on
imagesc(mat_acc(i,:))
cb = colorbar;
cb.Label.String = 'accuracy';
ylabel('State [#]')
xlabel('State duration [remaining samples]')
title('Accuracy')
set(gca,'fontsize',16)

subplot(2,2,2);hold on
imagesc(mat_ecc(i,:))
colormap('turbo')
cb = colorbar;
cb.Label.String = 'displacement';
title('Displacement')
set(gca,'fontsize',16)

subplot(2,2,3);hold on
for iCoh = 1:length(snr)
    comb = [i find(ncoh == snr(iCoh))];
    [uniqueNumbers, ~, indices] = unique(comb);
    counts = accumarray(indices, 1);
    indx = uniqueNumbers(counts > 1);

    scatter(slen(indx), nanmean(mat_acc(indx,end-nSample:end),2),'filled', 'MarkerFaceColor', col(iCoh,:), 'MarkerFaceAlpha', .5)
end
ylabel('JS accuracy')
ylim([0 1])
xlabel('State duration [remaining samples]')
set(gca,'fontsize',16)

subplot(2,2,4);hold on
for iCoh = 1:length(snr)
    comb = [i find(ncoh == snr(iCoh))];
    [uniqueNumbers, ~, indices] = unique(comb);
    counts = accumarray(indices, 1);
    indx = uniqueNumbers(counts > 1);
    scatter(slen(indx), nanmean(mat_ecc(indx,end-nSample:end),2),'filled', 'MarkerFaceColor', col(iCoh,:), 'MarkerFaceAlpha', .5)
end
ylim([0 1])
set(gca,'fontsize',16)

print(f, '/Users/fschneider/Desktop/duration_acc_vs_ecc', '-r500', '-dpng');

% figure; hold on
% for iL = 1:size(mat_acc,1)
% plot(mat_acc(iL,:),'Color',col(snr == ncoh(iL),:))
% end
%% Crosscorrelation between vector length and joystick dispalcement
nLag = 150;
smooth_win = 20;
clear xc
for iState = 1:size(frme_vec,2)
    clear v_smooth c_smooth_detrend ecc_detrend
    v_smooth = smoothdata(frme_vec{iState}.resultant_length,'gaussian',smooth_win);
    v_smooth_detrend = [0 v_smooth] - nanmean(v_smooth);
    v_smooth_detrend(isnan(v_smooth_detrend)) = 0;

    ecc_detrend = js_ecc{iState} - nanmean(js_ecc{iState});
    ecc_detrend(isnan(ecc_detrend)) = 0;


%     plot(xcorr(abs(v_smooth_detrend),abs(ecc_detrend),nLag));
%     plot(xcorr(v_smooth,js_ecc{iState},nLag));
    xc(iState,:) = xcorr(v_smooth_detrend,ecc_detrend,nLag,'normalized');
end
% print(f, '/Users/fschneider/Desktop/xcorr', '-r500', '-dpng');

f=figure
ax = subplot(2,1,1);hold on
imagesc(xc)
ln = line([0 0],[0 2000], 'Color', 'k','LineStyle', ':', 'LineWidth', 2);


ax = subplot(2,1,2);hold on
plot(mean(xc))
ln = line([0 0],[-1 1], 'Color', 'k','LineStyle', ':', 'LineWidth', 2);

% figure; hold on
% plot([0 v_smooth_detrend])
% plot(ecc_detrend)


%% Correlation between avg joystick displacement and avg vector length
nSample         = 29; % Time window size [samples]
indx        	= ~(cell2mat(cellfun(@length,js_ecc', 'UniformOutput', false)) <= nSample+1);
roi_ecc         = cell2mat(cellfun(@(x) x(end-nSample:end),js_ecc(indx)', 'UniformOutput', false));
roi_vec         = cell2mat(cellfun(@(x) x(end-nSample:end),vec_len(indx)', 'UniformOutput', false));
roi_coh         = cell2mat(cellfun(@(x) x(end-nSample:end),act_coh(indx)', 'UniformOutput', false));

e = mean(roi_ecc,2);
v = mean(roi_vec,2);
c = mean(roi_coh,2);
for i = 1:size(frme_vec,2); nc(i) = frme_vec{i}.nominal_coh;end

figure;hold on
for i = 1:length(e)
    scatter(e(i),v(i) ,'filled', 'MarkerFaceAlpha',.5, 'MarkerFaceColor', [c(i) 0 0])
end
xlabel('avg eccentricity')
ylabel('avg vector length')
set(gca, 'fontsize', 14)

%%% Correlation within coherence condition
f = figure('units','centimeters','position',[0 0 50 10]);
snr = unique(rdp_coh);
col = cool(7);
for iCoh = 1:7
    subplot(1,7,iCoh); hold on
    cidx = nc == snr(iCoh);
    scatter(v(cidx),e(cidx),'filled', 'MarkerFaceAlpha',.75, 'MarkerFaceColor', col(iCoh,:))
    lsline
    [r,pv]               	= corrcoef(v(cidx),e(cidx));
    ylim([0 1])
    ylabel('vector length')
    xlabel('displacement')
    title({['r: ' num2str(round(r(2),3))]; ['p: ' num2str(round(pv(2),3))]})
    set(gca, 'fontsize', 14)
end

%% Cycle plot

idx = 31:60;
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


%% angular error vs vector length

cmap = cool(256);
close all
figure
nSamples = 99;
time = 1:100;
for iState = 1:length(frme_vec)
    angular_error(iState,:) = frme_vec{iState}.resultant_ang_error(end-nSamples:end);
    vec_length(iState,:) = frme_vec{iState}.resultant_length(end-nSamples:end);
    ncoh(iState) = mean(frme_vec{iState}.actual_coherence);
    
    p3 = plot3(time, vec_length(iState,:),angular_error(iState,:)); hold on
    p3.Color = cmap(ceil(256*ncoh(iState)),:);
    
end

xlabel('time [frames]')
ylabel('length of resultant vector [coherence]')
zlabel('angular error between nominal and actual direction')
set(gca, 'fontsize', 20)
azimuth = 245; % Horizontal rotation (in degrees)
elevation = 25; % Vertical elevation (in degrees)
view(azimuth, elevation);

%% joystick error vs eccentricity

cmap = cool(256);
close all
figure
nSamples = 99;
time = 1:100;
pcoh = rdp_coh;
pcoh(pcoh==0) = 0.0001;
for iState = 1:length(frme_vec)
    if length(js_ecc{iState}) < nSamples+1
        continue
    end
    %     js_angular_error(iState,:) = rad2deg(circ_dist(deg2rad(js_dir{iState}(end-nSamples:end)), deg2rad(rdp_dir(iState))));
    %     p3 = plot3(time, js_ecc{iState}(end-nSamples:end),js_angular_error(iState,:)); hold on
    p3 = plot3(time, js_ecc{iState}(end-nSamples:end),frme_vec{iState}.resultant_ang_error(end-nSamples:end)); hold on
    p3.Color = [cmap(ceil(256*pcoh(iState)),:) .5];
    
end

xlabel('time [frames]')
ylabel('joystick eccentricity')
% zlabel('angular error between nominal and joystick direction')
zlabel('angular error between nominal and actual direction')
set(gca, 'fontsize', 20)
azimuth = 245; % Horizontal rotation (in degrees)
elevation = 25; % Vertical elevation (in degrees)
view(azimuth, elevation)

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
