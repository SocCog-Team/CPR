function out = CPR_create_random_walk(param)

if nargin < 1
    % Set parameters
    param.trial                     = 0;
    param.coh_duration_ms       	= 20000;                % Duration of coherence block
    param.walk_duration_ms      	= 60000;                % Duration of random walk (stimulus cycle)
    param.Fs                        = 1000 / 120;           % Screen sampling rate
    param.snr_list                  = [.2 .4 .6 .8];        % Stimulus coherence
    param.polar_step_size           = 0.1;                  % Step size in polar space
    param.reward_probability        = 0.01;                 % Probability of reward target appearance
    param.feedback_interval_ms      = param.Fs * 30;        % Interval between feedback presentation
end

% Initialize variables
out.trial                           = param.trial;          % Trial number
out.RDP_direction                   = rand * 2 * pi;      	% Random stimulus direction
out.feedback_ts                     = 1;                    % Initialise variable
out.RDP_coherence                   = [];                   % Initialise variable    
cnt                                 = 0;

% Calculate number of samples
nSamples                            = param.walk_duration_ms / param.Fs;                % Number of samples for random walk
nCohSamples                         = param.coh_duration_ms / param.Fs;                 % Number of samples for coherence block
nCohBlocks                          = param.walk_duration_ms/param.coh_duration_ms;     % Number of coheren block in stimulus cycle
out.RDP_coherence_ts                = [0 nCohSamples.*(1:nCohBlocks-1)] .* param.Fs;    % Timstamps of coherence change
out.RDP_coherence                   = param.snr_list(randi([1 length(param.snr_list)],1,3)); % Coherence values

% Avoid repeating coherence values
while sum(diff(out.RDP_coherence)==0) > 0
    out.RDP_coherence            	= param.snr_list(randi([1 length(param.snr_list)],1,3));
end

% Perform random walk in polar space
for iSample = 2:nSamples
    % Update stimulus direction with random step
    out.RDP_direction(iSample)  	= out.RDP_direction(iSample-1) + randn * param.polar_step_size;
    out.RDP_direction(iSample)    	= mod(out.RDP_direction(iSample), 2 * pi); % Wrap around if exceeding 2*pi
      
    % Check for reward target appearance
    if rand < param.reward_probability && (iSample-out.feedback_ts(end)) > param.feedback_interval_ms
        cnt                         = cnt+1;
        out.feedback_ts(cnt)        = iSample; % Frame of feedback presentation
    end
end
end