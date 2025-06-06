function out = CPR_stimulus_control(param)

if nargin < 1
    % Set parameters
    param.trial                   	= 0;
    param.NoStates                  = 50;
    param.NoCoherenceStates         = 10;
    param.state_min_ms              = 1250;
    param.state_max_ms              = 2500;
    param.snr_list                  = [.2 .4 .6 .8];
    param.directionChange_list      = [30 60 90];
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
    out.RDP_coherence(1)        	= param.snr_list(end);
else
    out.RDP_coherence(1)        	= param.snr_list(randi([1 length(param.snr_list)]));
end

% Calculate state-wise RDP parameters
for iState = 2:param.NoStates
    if ddir_sign(iState-1) == 1
        out.RDP_direction(iState)   = out.RDP_direction(iState-1) + param.directionChange_list(ddir_magn(iState-1));
    else
        out.RDP_direction(iState)   = out.RDP_direction(iState-1) - param.directionChange_list(ddir_magn(iState-1));
    end
    
    if mod(iState-1,param.NoCoherenceStates) == 0
        out.RDP_coherence(iState) 	= param.snr_list(randi([1 length(param.snr_list)]));
    else
        out.RDP_coherence(iState)   = out.RDP_coherence(iState-1);
    end      
end

% Normalise values to circular space
out.RDP_direction_raw               = out.RDP_direction;
out.RDP_direction                   = mod(out.RDP_direction,360);

% Draw timesteps of target appearances from trial onset, no target during  
% first and last X ms of trial
sample_dur                          = 1000 / param.Fs;
sample_vec                          = param.target_blocked_ms:sample_dur:sum(out.state_duration_ms)-param.target_blocked_ms; 

for iSmaple = sample_vec
    dice                            = rand;
    if dice < param.target_prob
        trg(iSmaple)               	= 1;
    else
        trg(iSmaple)             	= 0;
    end
end

% Exclude targets that are within ITI
out.target_ts                       = find(trg);
excl_Idx                            = [false diff(out.target_ts) < param.target_ITI_ms];
out.target_ts(excl_Idx)             = [];

end