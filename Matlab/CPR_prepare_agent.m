function AGNT = CPR_prepare_agent(STATE, AGNT)

% Initialise
tmp_rdp_dir             = [];
tmp_rdp_str             = [];

% Create samples-wise data vector
for iState = 1:length(STATE.RDP_direction)
    nSamples            = ceil(STATE.state_duration_ms(iState)/10);
    tmp_rdp_dir         = [tmp_rdp_dir repmat(STATE.RDP_direction_raw(iState),1,nSamples)];
    tmp_rdp_str         = [tmp_rdp_str repmat(STATE.RDP_coherence(iState),1,nSamples)];
end

% Zero-pad, add noise, smooth
pad_vec                 = zeros(1,AGNT.lag);
buffer                  = 49; % to avoid lack of samples during task
AGNT.dir_clean          = [tmp_rdp_dir tmp_rdp_dir(end-buffer:end)];
AGNT.str_clean          = [tmp_rdp_str tmp_rdp_str(end-buffer:end)];
AGNT.dir_nse            = [pad_vec AGNT.dir_clean + normrnd(0,AGNT.dir_sigma, [1,length(AGNT.dir_clean)])];
AGNT.str_nse            = [pad_vec AGNT.str_clean + normrnd(0,AGNT.str_sigma, [1,length(AGNT.str_clean)])];
AGNT.dir_smooth         = smoothdata(AGNT.dir_nse, AGNT.smooth_kernel, AGNT.win);
AGNT.str_smooth         = smoothdata(AGNT.str_nse, AGNT.smooth_kernel, AGNT.win);

% Keep strength between normalised limits
AGNT.str_smooth(AGNT.str_smooth > 1) = 1;
AGNT.str_smooth(AGNT.str_smooth < 0) = 0;

% Crop after zero-padding, return normalised values
% AGNT.dir_clean          = mod(AGNT.dir_clean,360);
% AGNT.str_clean          = mod(AGNT.str_clean,360);
% AGNT.dir_nse            = mod(AGNT.dir_nse(1:length(tmp_rdp_dir)),360);
% AGNT.str_nse            = mod(AGNT.str_nse(1:length(tmp_rdp_str)),360);
% AGNT.dir_smooth         = mod(AGNT.dir_smooth(1:length(tmp_rdp_dir)),360);
% AGNT.str_smooth         = mod(AGNT.str_smooth(1:length(tmp_rdp_str)),360);dots_only

end

% % Log-space between samples - doesn't work as expected
% if strcmp(x,'log')
%     AGNT.dir_log = [];
%     idx = [1 find(diff(AGNT.dir_clean))+1];
%     
%     for i = 2:length(idx)
%         AGNT.dir_log = [AGNT.dir_log exp(linspace(log(AGNT.dir_clean(idx(i-1))),log(AGNT.dir_clean(idx(i))),idx(i)-idx(i-1)))];
%     end
% end
