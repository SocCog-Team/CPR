function [out_p1, out_p2] = CPR_extract_response_profile(d,idx,exp_info)

out_p1                  = extract_response(d, idx, exp_info, 1);
out_p2                  = extract_response(d, idx, exp_info, 2);

out_p1                  = get_outcome_summary(d, idx,1,out_p1);
out_p2                  = get_outcome_summary(d, idx,2,out_p2);

out_p1                  = sort_coherence_blocks(out_p1);
out_p2                  = sort_coherence_blocks(out_p2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

cnt                         = 0;
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
   
    % Loop through all coherence blocks in this cycle
    for iCoh = 2:length(tmp_coh_ts)
        
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
        cnt                 = cnt+1;
        out.block_num(cnt) 	= cnt;
        out.coherence(cnt) 	= tmp_coh(iCoh);
        out.acc_raw{cnt}    = js_acc(cIdx);
        out.ecc_raw{cnt}    = js_ecc(cIdx);
        out.acc_mean(cnt)   = mean(js_acc(cIdx));
        out.ecc_mean(cnt)   = mean(js_ecc(cIdx));
        out.o_corr(cnt)   	= cc;
        out.x_corr{cnt}     = sxc;
        out.lag(cnt)    	= peak_sample * double(median(diff(frame_ts)))/1e3; % [ms]
        
        % Save raw signals
        out.raw.ts{cnt}     	= frame_ts;
        out.raw.ts{cnt}     	= frame_ts;
        out.raw.rdp_dir{cnt} 	= rdp_dir;
        out.raw.js_dir{cnt}  	= js_dir;
        out.raw.js_dev{cnt}   	= js_dev;
        out.raw.js_acc{cnt}   	= js_acc;
        out.raw.js_ecc{cnt}   	= js_ecc;
    end
end
end

function out = sort_coherence_blocks(in)
% Average for coherence level
snr                         = unique(in.coherence);
for iCoh = 1:length(snr)
    cindx                          = in.coherence == snr(iCoh);
    in.coh_sorted.acc_mean(iCoh)   = mean(in.acc_mean(cindx));
    in.coh_sorted.ecc_mean(iCoh)   = mean(in.ecc_mean(cindx));
    in.coh_sorted.o_corr(iCoh)     = mean(in.o_corr(cindx));
    in.coh_sorted.lag(iCoh)        = mean(in.lag(cindx));
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