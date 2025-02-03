addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/preprocessing
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/

% cd /Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_social_context/
cd /Users/fschneider/Desktop/randomwalk_pilot/

% Configuration file
cfg_pth = '/Users/fschneider/Documents/GitHub/CPR/Matlab/CFG/felix_context_study.cfg';

% Variables to import
var_import = {
    'INFO_', ...
    'TRIAL_', ...
    'IO_joystickDirection', ...
    'IO_joystickStrength',...
    'IO_fixation_flag',...
    'EYE_x_dva',...
    'EYE_y_dva',...
    '#stimDisplay'};

file_list = {'20250121_Sebastian_CPRsolo_randomwalk.mwk2';
    '20250123_AnnKathrin_CPRsolo_randomwalk.mwk2';
    '20250128_Pinar_CPRsolo_randomwalk.mwk2';
    '20250128_Anahita_CPRsolo_randomwalk.mwk2';
    '20250128_Felix_CPRsolo_randomwalk.mwk2';
    '20250128_Max_CPRsolo_randomwalk.mwk2';};

nSample                     = 29;                                           % Time window size [samples]
nLag                        = 150;                                          % Cross-correlation lag


% Convert each CPR .mwk2 file to .h5
for iFile = 1:length(file_list)
    % Import and write file
    d = CPR_import_mwk2(file_list{iFile}, var_import, true, cfg_pth);
    exp_info                    = strsplit(file_list{iFile},'_');
    
    % Index variables of interest
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
    idx.cum_score            	= d.event == 'INFO_Score';
    idx.trg_score              	= d.event == 'INFO_TargetScore';
    idx.bonus1                  = d.event == 'INFO_bonus_ply1_cents';
    idx.performance             = d.event == 'INFO_performance_percent';
    
    % Extract subject behavior
    out                         = extract_response(d, idx, exp_info, 1);
    out                         = get_outcome_summary(d, idx,1,out);
    out                         = sort_coherence_blocks(out);
    
    % Save summary file
    save(['/Users/fschneider/Desktop/randomwalk_pilot/' exp_info{2} '.mat'], 'out', '-v7.3')
    
    acc(iFile,:)                = out.coh_sorted.acc_mean;
    ecc(iFile,:)                = out.coh_sorted.ecc_mean;
    o_corr(iFile,:)             = out.coh_sorted.o_corr;
    lag(iFile,:)                = out.coh_sorted.lag;
end

%%
snr_lst = round(unique(out.coherence).*100);
col = turbo(size(acc,1))
figure; hold on
for iFile = 1:size(acc,1)
    plot(snr_lst,acc(iFile,:),'LineStyle','-', 'Marker','x', 'LineWidth',2,'Color', col(iFile,:))
%     scatter(snr_lst,acc(iFile,:), 'Marker','x', 'SizeData',30,'MarkerFaceColor', col(iFile,:))
%     f=fit(snr_lst',acc(iFile,:)','poly2');
%     plot(snr_lst',acc(iFile,:)')   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function in_struct = get_outcome_summary(d, idx, player_flag,in_struct)
if player_flag == 1
    score = d.value(idx.cum_score);
    bonus = d.value(idx.bonus1);
    outc = d.value(idx.outcome);
    n_hit = sum(cellfun(@(x) strcmp(x,'hit'), outc));
    n_miss = sum(cellfun(@(x) strcmp(x,'miss'), outc));
elseif player_flag == 2
    score = d.value(idx.cum_score2);
    bonus = d.value(idx.bonus2);
    outc = d.value(idx.outcome2);
    n_hit = sum(cellfun(@(x) strcmp(x,'hit'), outc));
    n_miss = sum(cellfun(@(x) strcmp(x,'miss'), outc));
end

in_struct.score = score{end};
in_struct.bonus = bonus{end};
in_struct.hir = n_hit/(n_hit+n_miss);
end

function out = extract_response(d, idx, exp_info, player_flag)

% Get cycle timestamps
cyc.cOn                     = d.time(idx.cOn);
cyc.cEnd                    = d.time(idx.cEnd);

ccnt                        = 0;
nLag                        = 150;
out                         = add_exp_info(exp_info,player_flag);

for iCyc = 1:length(cyc.cEnd)-1
    % Trial index
    cycIdx                  = [];
    cycIdx                  = d.time >= cyc.cOn(iCyc) & d.time <= cyc.cEnd(iCyc);
    
    % Extract frame times
    frame_ts                    = getTrialData(d.time, cycIdx, idx.frame);                      % Frame timestamps
    
    % Extract coherence blocks
    tmp_coh                     = getTrialData(d.value, cycIdx, idx.RDP_coh);                   % RDP coherence
    tmp_coh_ts                  = getTrialData(d.time, cycIdx, idx.RDP_coh);                    % RDP coherence timestamps
    
    % Extract stimulus direction
    tmp_rdp_dir                 = getTrialData(d.value, cycIdx, idx.RDP_dir);                   % RDP_direction
    tmp_rdp_dir_ts              = getTrialData(d.time, cycIdx, idx.RDP_dir);                    % RDP_direction timestamps
    
    % Joystick data
    if player_flag == 1
        tmp_js_dir_ts        	= getTrialData(d.time, cycIdx, idx.JS_dir);                     % JS player1
        tmp_js_dir            	= getTrialData(d.value, cycIdx, idx.JS_dir);
        tmp_js_ecc           	= getTrialData(d.value, cycIdx, idx.JS_str);
    elseif player_flag == 2
        tmp_js_dir_ts        	= getTrialData(d.time, cycIdx, idx.JS2_dir);                    % JS player2
        tmp_js_dir            	= getTrialData(d.value, cycIdx, idx.JS2_dir);
        tmp_js_ecc           	= getTrialData(d.value, cycIdx, idx.JS2_str);
    end
    
    % Initialise variables
    clear rdp_dir js_dir js_dev js_ecc js_acc js_dff rdp_dff
    rdp_dir                 = nan(1,length(frame_ts));
    js_dir                  = nan(1,length(frame_ts));
    js_ecc                  = nan(1,length(frame_ts));
    
    % Create frame-wise data
    for iFrame = 1:length(frame_ts)
        if sum(tmp_rdp_dir_ts < frame_ts(iFrame)) ~=0
            rdp_dir(iFrame) = tmp_rdp_dir(find(tmp_rdp_dir_ts < frame_ts(iFrame),1,'last'));
            js_dir(iFrame)  = tmp_js_dir(find(tmp_js_dir_ts < frame_ts(iFrame),1,'last'));
            js_ecc(iFrame)  = tmp_js_ecc(find(tmp_js_dir_ts < frame_ts(iFrame),1,'last'));
        end
    end
    
    % Calculate accuracy
    js_dev                  = rad2deg(circ_dist(deg2rad(js_dir),deg2rad(rdp_dir)));         % Get circular distance to RDP direction
    js_acc                	= abs(1 - abs(js_dev / 180));                                   % Calculate accuracy
    
    % Sample-by-sample difference (Derivative)
    js_dff                  = [0 rad2deg(circ_dist(deg2rad(js_dir(1:end-1)),deg2rad(js_dir(2:end))))];
    rdp_dff                 = [0 rad2deg(circ_dist(deg2rad(rdp_dir(1:end-1)),deg2rad(rdp_dir(2:end))))];
    
    ex                    	= isnan(js_dff) | isnan(rdp_dff);
    js_dff(ex)            	= 0;
    rdp_dff(ex)            	= 0;
    
    % Save raw signals
    out.raw.ts{iCyc}     	= frame_ts;
    out.raw.rdp_dir{iCyc} 	= rdp_dir;
    out.raw.js_dir{iCyc}  	= js_dir;
    out.raw.js_dev{iCyc}   	= js_dev;
    out.raw.js_acc{iCyc}   	= js_acc;
    out.raw.js_ecc{iCyc}   	= js_ecc;
    
    % Loop through all coherence blocks in this cycle
    for iCoh = 1:length(tmp_coh_ts)
        
        % Coherence block index
        if iCoh == length(tmp_coh_ts)
            cIdx        	= frame_ts >= tmp_coh_ts(iCoh);
        else
            cIdx         	= frame_ts >= tmp_coh_ts(iCoh) & frame_ts < tmp_coh_ts(iCoh+1);
        end
        
        if sum(cIdx) < 100
            continue
        end
        
        %%% Calculate correlations between stimulus and joystick %%%
        % Circular correlation
        cc                  = circ_corrcc(deg2rad(js_dir(cIdx)),deg2rad(rdp_dir(cIdx)));
        % Cross correlation
        xc                  = xcorr(abs(js_dff(cIdx)),abs(rdp_dff(cIdx)),nLag,'normalized');
        sxc                 = smoothdata(xc,'gaussian',20);                 % Smooth correlation output with gaussian kernel
        max_corr_coef    	= max(sxc(:,nLag:end));                        	% Max cross-correlation coefficient
        peak_sample         = find(sxc(:,nLag:end) == max_corr_coef);     	% Peak position of cross-correlation
        
        % Store output values for each coherence block
        ccnt                = ccnt+1;
        out.block_num(ccnt) = ccnt;
        out.coherence(ccnt) = tmp_coh(iCoh);
        out.acc_raw{ccnt}   = js_acc(cIdx);
        out.ecc_raw{ccnt}   = js_ecc(cIdx);
        out.acc_mean(ccnt)  = nanmean(js_acc(cIdx));
        out.ecc_mean(ccnt)  = nanmean(js_ecc(cIdx));
        out.o_corr(ccnt)   	= cc;
        out.x_corr{ccnt}    = sxc;
        out.lag(ccnt)    	= peak_sample * double(median(diff(frame_ts)))/1e3; % [ms]
        
    end
end
end

function out = sort_coherence_blocks(in)
% Average for coherence level
snr                         = unique(in.coherence);

for iCoh = 1:length(snr)
    cindx                          = in.coherence == snr(iCoh);
    in.coh_sorted.acc_mean(iCoh)   = nanmean(in.acc_mean(cindx));
    in.coh_sorted.ecc_mean(iCoh)   = nanmean(in.ecc_mean(cindx));
    in.coh_sorted.o_corr(iCoh)     = nanmean(in.o_corr(cindx));
    in.coh_sorted.lag(iCoh)        = nanmean(in.lag(cindx));
end

out = in;

end

function out = add_exp_info(exp_info,player_flag)

out.date                    = exp_info{1};
out.dyad                    = exp_info{2};

if player_flag == 1
    out.player_id        	= exp_info{2}(1:3);
    out.setup               = 'psy4';
elseif player_flag == 2
    out.player_id        	= exp_info{2}(4:end);
    out.setup               = 'psy3';
end

out.condition               = exp_info{3};
out.block                   = exp_info{4}(1:end-3);


if strcmp(out.condition, 'CPRcoopration')
    out.condition           = 'CPRcooperation';
end
end