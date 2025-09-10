function out = CPR_stimulus_control(param)

if nargin < 1
    % Set parameters
    param.trial                   	= 0;
    param.NoStates                  = 10;
    param.NoCoherenceStates         = 5;
    param.state_min_ms              = 1250;
    param.state_max_ms              = 2500;
    param.snr_list                  = [.2 .4 .6 .8];
    param.directionChange_list      = [90,135,180];
    param.target_ITI_ms             = 300;
    param.target_blocked_ms         = 500;
    param.target_prob               = .005;
    param.Fs                        = 100;
end

% Draw trial duration
out.trl_no                          = param.trial;
out.state_duration_ms               = randi([param.state_min_ms param.state_max_ms],1,param.NoStates);

% Randomise sign and magnitude of state-wise direction change
ddir_sign                           = randi([0 1],1,param.NoStates-1);
ddir_magn                           = randi([1 3],1,param.NoStates-1);

% Determine first stimulus state parameter
out.RDP_direction(1)                = randi([0 359]);

if param.trial == 0
    out.RDP_coherence(1)        	= param.snr_list(end);                                                                  % Get last entry of list => highest value
    out.RDP_coherence_cnt(1,:)   	= param.snr_list;                                                                       % Set up counter array
    out.RDP_coherence_cnt(2,:)      = zeros(1, length(param.snr_list));
    out.RDP_coherence_cnt(2,:)      = out.RDP_coherence_cnt(1,:) == out.RDP_coherence(1);
else
    tmp                             = load([param.pth 'RDP_coherence_cnt.txt']);                                            % Import .txt file
    out.RDP_coherence_cnt           = reshape(tmp, [2 length(tmp)/2]);                                                      % Reshape to matrix  
    [val,ind]                       = sort(out.RDP_coherence_cnt(2,:));                                                     % Sort according to count
    coh_pool                       	= ind(val == min(val));                                                                 % Extract indizes with minimum value
    out.RDP_coherence(1)        	= param.snr_list(coh_pool(randi([1 length(coh_pool)])));                                % Draw from this coherence pool
    out.RDP_coherence_cnt(2,:)      = out.RDP_coherence_cnt(2,:) + (out.RDP_coherence_cnt(1,:) == out.RDP_coherence(1));    % Add to coherence count
end

% Calculate state-wise RDP parameters
for iState = 2:param.NoStates
    if ddir_sign(iState-1) == 1
        out.RDP_direction(iState)   = out.RDP_direction(iState-1) + param.directionChange_list(ddir_magn(iState-1));
    else
        out.RDP_direction(iState)   = out.RDP_direction(iState-1) - param.directionChange_list(ddir_magn(iState-1));
    end
    
    if mod(iState-1,param.NoCoherenceStates) == 0
        [val,ind]                  	= sort(out.RDP_coherence_cnt(2,:));
        coh_pool                   	= ind(val == min(val));   
        out.RDP_coherence(iState) 	= param.snr_list(coh_pool(randi([1 length(coh_pool)])));
        out.RDP_coherence_cnt(2,:) 	= out.RDP_coherence_cnt(2,:) + (out.RDP_coherence_cnt(1,:) == out.RDP_coherence(iState));
    else
        out.RDP_coherence(iState)   = out.RDP_coherence(iState-1);
        out.RDP_coherence_cnt(2,:)  = out.RDP_coherence_cnt(2,:) + (out.RDP_coherence_cnt(1,:) == out.RDP_coherence(iState));
    end      
end

% Normalise values to circular space
out.RDP_direction_raw               = out.RDP_direction;
out.RDP_direction                   = mod(out.RDP_direction,360);

% % Draw timesteps of target appearances from trial onset, no target during  
% % first and last X ms of trial
% sample_dur                          = 1000 / param.Fs;
% sample_vec                          = param.target_blocked_ms:sample_dur:sum(out.state_duration_ms)-param.target_blocked_ms; 
% 
% for iSmaple = sample_vec
%     dice                            = rand;
%     if dice < param.target_prob
%         trg(iSmaple)               	= 1;
%     else
%         trg(iSmaple)             	= 0;
%     end
% end
% 
% % Exclude targets that are within ITI
% out.target_ts                       = find(trg);
% excl_Idx                            = [false diff(out.target_ts) < param.target_ITI_ms];
% out.target_ts(excl_Idx)             = [];

end