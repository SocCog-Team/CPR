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
                        
                        % Extract timestamp of first target
                        tmp_trg_on          = getTrialData(d.value, ssIdx, idx.trg_on); 
                        tmp_trg_ts          = getTrialData(d.time, ssIdx, idx.trg_on);
                        try
                            onset_all_ts{cc}  	= tmp_trg_ts(logical(tmp_trg_on));
                            onset_all_ms        = double(tmp_trg_ts(logical(tmp_trg_on)) - dp_ts(1)) ./1e3;
                            trg1_onset_ms(cc)   = onset_all_ms(1);
                        catch
                            onset_all_ts{cc}  	= nan;
                            trg1_onset_ms(cc)   = nan;
                        end
                        
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
                            
                            % Frame before forst target?
                            if dp_ts(iFrme) < onset_all_ts{cc}(1)
                                trg1_index{cc}(iFrme)	= true;             % prior to first target
                            else
                                trg1_index{cc}(iFrme)	= false;            % after first target
                            end
                        end
                    end
                end
            end
        end
    end
    
    save(['/Users/fschneider/Desktop/resultant_vec/vector_subj_'  sbj_lst{iSubj}], 'frme_vec', '-v7.3')
    save(['/Users/fschneider/Desktop/resultant_vec/trg_subj_'  sbj_lst{iSubj}],  'trg1_index','trg1_onset_ms', '-v7.3')
    save(['/Users/fschneider/Desktop/resultant_vec/joystick_subj_' sbj_lst{iSubj}], 'js_dir','js_ecc', '-v7.3')
end
        
load('/Users/fschneider/Desktop/subj_lst.mat', 'sbj_lst', '-v7.3')

%% Plot timeline of joystick response; Extract physical coherence level
%%% Remove all sample after first target appearance!
% Show only 100 samples before direction change

clear avg_ecc avg_acc avg_pcoh

load('/Users/fschneider/Desktop/resultant_vec/subj_lst.mat')

win                         = 60;
nSample                     = 29;

for iSubj = 1:length(sbj_lst)
    clear ncoh slen ndir
    mat_acc                 = nan(10000,300);
    mat_ecc                 = nan(10000,300);
    mat_pcoh                = nan(10000,300);
    c = 0;
    
    load(['/Users/fschneider/Desktop/resultant_vec/vector_subj_' sbj_lst{iSubj} '.mat'])
    load(['/Users/fschneider/Desktop/resultant_vec/joystick_subj_' sbj_lst{iSubj} '.mat'])
        
    for iState = 1:size(frme_vec,2)
        
        if length(js_ecc{iState}) >= win && length(js_ecc{iState}) == length(frme_vec{iState}.resultant_ang_error)+1
            % Nominal parameters
            c           = c+1;
            slen(c)     = length(js_ecc{iState}(win:end));
            ncoh(c)    	= frme_vec{iState}.nominal_coh;
            ndir(c)    	= frme_vec{iState}.nominal_dir_deg;
            
            % Joystick response
            js_dir_error = rad2deg(circ_dist(deg2rad(js_dir{iState}), deg2rad(frme_vec{iState}.nominal_dir_deg)));
            joystick_acc = abs(1 - abs(js_dir_error) / 180);
            mat_acc(c,end-slen(c)+1:end)    = joystick_acc(win:end);
            mat_ecc(c,end-slen(c)+1:end) 	= js_ecc{iState}(win:end);
            
            % Physical coherence level
            mat_pcoh(c,end-slen(c)+1:end)  	= frme_vec{iState}.actual_coherence(win-1:end); % shorter bcs of pairwise frame computation

        end
    end
    
    snr_lst = unique(ncoh);
    if length(snr_lst)>7
       snr_lst = snr_lst(1:2:end);
    end
    
    mat_acc(sum(isnan(mat_acc),2) == size(mat_acc,2),:) = [];
    mat_ecc(sum(isnan(mat_ecc),2) == size(mat_ecc,2),:) = [];
    mat_pcoh(sum(isnan(mat_pcoh),2) == size(mat_pcoh,2),:) = [];
    
    for iCoh = 1:length(snr_lst)
        cidx                            = ncoh == snr_lst(iCoh);
        avg_acc(iSubj,iCoh,:)           = nanmean(mat_acc(cidx,:));
        avg_ecc(iSubj,iCoh,:)           = nanmean(mat_ecc(cidx,:));
        
        %%% Median split %%%
        clear pcoh_data avg_ecc_time_win idx_larger_than_median
        ecc_data                            = nanmean(mat_ecc(cidx,end-nSample:end),2);
        pcoh_data                           = nanmean(mat_pcoh(cidx,end-nSample:end),2);
        acc_data                            = nanmean(mat_acc(cidx,end-nSample:end),2);
        idx_larger_than_median              = ecc_data > median(ecc_data);
 
        avg_pcoh_larger_median{iSubj,iCoh}  = pcoh_data(idx_larger_than_median);
        avg_pcoh_smaller_median{iSubj,iCoh} = pcoh_data(~idx_larger_than_median);
        
        avg_acc_larger_median{iSubj,iCoh}   = acc_data(idx_larger_than_median);
        avg_acc_smaller_median{iSubj,iCoh}  = acc_data(~idx_larger_than_median);
    end
end

%%% Plot timeline %%%
f = figure;
col = cool(length(snr_lst));
ax1 = subplot(1,2,1); hold on
ax2 = subplot(1,2,2); hold on

for iCoh = 1:length(snr_lst)
   
    [ax1, pl(iCoh)] = plot_averages(ax1,squeeze(avg_acc(:,iCoh,:)),col(iCoh, :));
    [ax2, pl(iCoh)] = plot_averages(ax2,squeeze(avg_ecc(:,iCoh,:)),col(iCoh, :));

end

ax1 = adjust_axes(ax1,'accuracy [norm]');
ax2 = adjust_axes(ax2,'tilt [norm]');

print(f, '/Users/fschneider/Desktop/avg_timeline_pop', '-r500', '-dpng');

%% Physical coherence for high vs low joystick tilt

iplot   = 0;

for iSubj = 1:length(sbj_lst)
    for iCoh = 1:length(snr_lst)
        [p(iSubj,iCoh),h(iSubj,iCoh)] = ranksum(avg_pcoh_larger_median{iSubj,iCoh},avg_pcoh_smaller_median{iSubj,iCoh});
        pcoh_group_average(iSubj,iCoh,:) = [median(avg_pcoh_larger_median{iSubj,iCoh}) median(avg_pcoh_smaller_median{iSubj,iCoh})];
    end
end

group_dff = pcoh_group_average(:,:,1) - pcoh_group_average(:,:,2);

%%% Correction %%%
sum(sum(p < .05/(size(p,1)*size(p,2)))) % Number of significant results after Bonferroni correction

%%% PLOT %%%
f = figure;hold on
ln = line([0 0],[0 20], 'Color', 'k','LineStyle', ':', 'LineWidth', 2);
edges = -.002 :.00025: .002;
for iCoh = 1:length(snr_lst)
    [pt(iCoh), ht(iCoh)] = ttest(group_dff(:,iCoh));
    
    histogram(group_dff(:,iCoh),edges,'FaceColor',col(iCoh,:),'FaceAlpha',.2,'EdgeAlpha',0)
end

ax = gca;
ax.YLabel.String = '#';
ax.XLabel.String = 'Median coherence difference';
ax.FontSize = 14;
ax.XTickLabelRotation = 0;

print(f, '/Users/fschneider/Desktop/physical_coherence_dff', '-r500', '-dpng');

%%% Check significant example
mean(group_dff(:,logical(pt)))
figure;histogram(group_dff(:,logical(pt)),20)

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

function ax = plot_coh_histograms(ax, iplot, col, in_larger, in_smaller)
ofs = .012;
step = .002;
edges = median(in_larger)-ofs:step:median(in_larger)+ofs;

h1 = histogram(in_larger,edges);
h1.FaceAlpha = .35;
h1.FaceColor = col(iplot,:);
h1.EdgeAlpha = 0;

h2 = histogram(in_smaller,edges);
h2.FaceAlpha = .35;
h2.FaceColor = [0 0 0];
h2.EdgeAlpha = 0;

if iplot == 1
    ax.YLabel.String = '#';
end

if iplot == 4
    ax.XLabel.String = 'Physical Coherence';
end

if iplot == 7
    lg = legend('high tilt','low tilt', 'Location', 'northwest');
    lg.FontSize = 10;
    lg.Position(2) = .46;
end

ax.XTick = round(median(in_larger),3);
ax.FontSize = 14;
ax.XTickLabelRotation = 0;
end

function ax = plot_acc_histograms(ax, iplot, col, in_larger, in_smaller)
edges                   = 0:.05:1;

h1 = histogram(in_larger,edges);
h1.FaceAlpha = .35;
h1.FaceColor = col(iplot,:);
h1.EdgeAlpha = 0;

h2 = histogram(in_smaller,edges);
h2.FaceAlpha = .35;
h2.FaceColor = [0 0 0];
h2.EdgeAlpha = 0;

if iplot == 1
    ax.YLabel.String = '#';
end

if iplot == 4
    ax.XLabel.String = 'Accuracy';
end
ax.FontSize = 14;
ax.XTickLabelRotation = 0;
end