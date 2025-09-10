function [STIM, AGNT] = CPR_create_replay_agnt(solo_summary_file,STIM,cycIdx)
%%% Reproduce human performance as replay for dyadic experiments %%%

% Extract info from upcoming cycle
iCyc = cycIdx + 1;

% Import frame-wise behavior from solo session
load(solo_summary_file)

% Initialise variables
n = 7250; % Length/cycle (1min 7200 samples + 50 buffer for safety)
replay.rdp_dir = nan([n 1]);
replay.rdp_coh = nan([n 1]);
replay.js_dir = nan([n 1]);
replay.js_tlt = nan([n 1]);

% Extract stimulus and response
% RDP direction
replay.rdp_dir(1:length(out.raw.js_dir{iCyc}),1) = out.raw.rdp_dir{iCyc};
replay.rdp_dir(:,1) = replace_nan_zero(replay.rdp_dir);

% RDP coherence
coh_vec = ((iCyc * 6)-5):1:(iCyc * 6);
replay.coherence_lst = out.coherence(coh_vec);

% Joystick direction
replay.js_dir(1:length(out.raw.js_dir{iCyc}),1) = out.raw.js_dir{iCyc};
replay.js_dir(:,1) = replace_nan_zero(replay.js_dir);

%  Joystick tilt
replay.js_tlt(1:length(out.raw.js_dir{iCyc}),1) = out.raw.js_ecc{iCyc};
replay.js_tlt(:,1) = replace_nan_zero(replay.js_tlt);

% Feedback timestamps
replay.feedback_ts = out.raw.trg_smple{iCyc};
replay.feedback_ts(replay.feedback_ts > 7200) = 7150;

% Add 3 missing targets manually here to show 300 targets in total
if iCyc == 2
    replay.feedback_ts = [replay.feedback_ts 1000 2000 5000]; 
end
replay.feedback_ts = sort(replay.feedback_ts);

% tmp_coh = [];
% for i = coh_vec
%     tmp_coh = [tmp_coh; repmat(out.coherence(i),1200,1)];
% end
% 
% replay.rdp_coh(1:length(tmp_coh),iCyc) = tmp_coh;
% replay.rdp_coh(7201:n,:) = repmat(replay.rdp_coh(7200,:),50,1);

% Export to respective format. Overwrite prior stimulus params
STIM.RDP_direction_deg  = replay.rdp_dir;
STIM.RDP_direction_rad  = deg2rad(replay.rdp_dir);
STIM.RDP_coherence      = [replay.coherence_lst 999e3]; % last entry to avoid MWorks list index error
STIM.feedback_ts        = [replay.feedback_ts 999e3]; % last entry to avoid MWorks list index error
AGNT.dir_smooth         = replay.js_dir;
AGNT.str_smooth         = replay.js_tlt;

% close all
% figure
% plot(replay.rdp_dir)
% hold on
% plot(replay.js_dir)
% figure
% plot(replay.js_tlt)

end

function out = replace_nan_zero(in)
v = in;

% Find the first valid entry (non-NaN, non-zero)
first_valid_idx = find(~isnan(v) & v ~= 0, 1, 'first');
first_valid_value = v(first_valid_idx);

% Replace leading NaNs with first valid value
v(1:first_valid_idx-1) = first_valid_value;

% Find the last valid entry (non-NaN, non-zero)
last_valid_idx = find(~isnan(v) & v ~= 0, 1, 'last');
last_valid_value = v(last_valid_idx);

% Replace trailing zeros with last valid value
v(last_valid_idx+1:end) = last_valid_value;
out = v;
end