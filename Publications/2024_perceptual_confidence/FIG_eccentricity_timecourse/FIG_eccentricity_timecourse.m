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
for iSubj = 1%:length(sbj_lst)
    
    disp(['Processing subject: ' sbj_lst{iSubj}])
    data_pth                     	= [pth sbj_lst{iSubj} '/raw/'];         % Data path
    
    if isfolder(data_pth)
        cd(data_pth)
        h5_files               	= dir('*.h5');                              % Get all .mat files in directory
        
        % For all files in directory
        for iFile = 3%1:length(h5_files)
            
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
                for iCyc = 2%1:length(cyc.cEnd)
                    
                    % Cycle index
                    cycIdx                  = [];
                    cycIdx                  = d.time >= cyc.cOn(iCyc) & d.time <= cyc.cEnd(iCyc);
                    
                    % Extract nominal stimulus direction
                    tmp_dir_ts              = getTrialData(d.time, cycIdx, idx.RDP_dir);                    % RDP_direction timestamps
                    tmp_dir                 = getTrialData(d.value, cycIdx, idx.RDP_dir);                   % RDP_direction
                    tmp_coh            	= getTrialData(d.value, cycIdx, idx.RDP_coh);                       % RDP coherence
                    tmp_coh_ts        	= getTrialData(d.time, cycIdx, idx.RDP_coh);                        % RDP coherence timestamps
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

                        %%%% TIME WINDOW ANALYSIS %%%
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
            end
        end
    end
end

%%
% figure; cc = 5;
% plot(frme_vec{cc}.resultant_ts, vec_len{cc})
% hold on
% plot(double(tmp_js_ts{cc}-tmp_js_ts{cc}(1))./1e6, tmp_js_ecc{cc})

nSample                     = 99;                                           % Time window size [samples]
roi_ecc = cell2mat(cellfun(@(x) x(end-nSample:end),js_ecc', 'UniformOutput', false));
roi_vec = cell2mat(cellfun(@(x) x(end-nSample:end),vec_len', 'UniformOutput', false));
roi_coh = cell2mat(cellfun(@(x) x(end-nSample:end),act_coh', 'UniformOutput', false));

figure
subplot(2,3,1)
imagesc(roi_ecc)
xlabel('# frames before direction change')
ylabel('# stimulus state')
title('eccentricity')
set(gca, 'fontsize', 14)
colorbar

subplot(2,3,2)
imagesc(roi_vec)
xlabel('# frames before direction change')
ylabel('# stimulus state')
title('vector length')
set(gca, 'fontsize', 14)
colorbar

subplot(2,3,3)
imagesc(roi_coh)
xlabel('# frames before direction change')
ylabel('# stimulus state')
title('actual coherence')
set(gca, 'fontsize', 14)
colorbar

subplot(2,3,4)
imagesc(roi_ecc ./ mean(roi_ecc,2))
xlabel('# frames before direction change')
ylabel('# stimulus state')
set(gca, 'fontsize', 14)
colorbar

subplot(2,3,5)
imagesc(roi_vec ./ mean(roi_vec,2))
xlabel('# frames before direction change')
ylabel('# stimulus state')
set(gca, 'fontsize', 14)
colorbar

subplot(2,3,6)
imagesc(roi_coh ./ mean(roi_coh,2))
xlabel('# frames before direction change')
ylabel('# stimulus state')
set(gca, 'fontsize', 14)
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = extract_resultant_vector(dp,dp_ts,nominal_dir)

good_frme               = 1;
c                       = 0;

for iFrme = 2:length(dp)
    if sum(isnan(dp{iFrme}))
        continue
    end
    
    % Dot position
    c                   = c+1;
    xIdx                = logical(mod([1:size(dp{iFrme},1)],2));
    xpos                = dp{iFrme}(xIdx);
    ypos                = dp{iFrme}(~xIdx);
    xpos_last           = dp{good_frme}(xIdx);                              % Last frame: x-position
    ypos_last         	= dp{good_frme}(~xIdx);                             % Last frame: y-position
    dts(c)              = dp_ts(iFrme);                                     % Dot Timestamp
    
    for iDot = 1:length(ypos)
        vs              = [xpos_last(iDot), ypos_last(iDot)];               % dot position last frame
        ve              = [xpos(iDot), ypos(iDot)];                         % dot position this frame
        dist(iDot)      = pdist([vs;ve],'euclidean');                       % distance/vector length between points
        
        delta               = ve - vs;                                     	% vector that points from x1 to x2
        dot_dir{c}(iDot)	= mod(atan2d(delta(1),delta(2)),360);           % Alternative: rad2deg(cart2pol(delta(2), delta(1)))
        dt{c}(iDot,:)       = delta;
    end
    
    dot_idx{c}              = dist >= median(dist)-.0001 & dist <= median(dist)+.0001; % Exclude new dots appearances
    good_frme               = iFrme;
end

ofs = 1;
clear resultant res_length res_ts rcoh
for iVec = 1:length(dt)
    % Avg frame direction
    mdf                 	= mean(dt{iVec}(dot_idx{iVec},:));          	% Mean x/y of all dots that didn't jump
    resultant(iVec)     	= mod(atan2d(mdf(1),mdf(2)),360);               % Resultant vector [deg]
    res_length(iVec)     	= pdist([[0 0];mdf],'euclidean');               % motion strength - length of vector
    res_ts(iVec)        	= double(dts(iVec)-dts(1))/1e6;
    
    % Check coherence level based on dot direction
    dots                    = dot_dir{iVec}(dot_idx{iVec});
    signal                  = abs(rad2deg(circ_dist(deg2rad(dots),deg2rad(nominal_dir)))) < ofs; % Number of dots moving in nominal direction
    rcoh(iVec)          	= sum(signal)/length(signal);
end

out.resultant_deg           = resultant;
out.resultant_length        = res_length;
out.resultant_ts            = res_ts;
out.actual_coherence        = rcoh;

end