function out = CPR_create_random_walk_v2(param)

if nargin < 1
    % Set parameters
    param.trial                     = 0;
    param.coh_duration_ms       	= 10000;                % Duration of coherence block 
    param.walk_duration_ms      	= 60000;                % Duration of random walk
    param.Fs                        = 1000 / 120;           % Screen sampling rate [Hz]
    param.snr_list                  = [.7 .75 .8 .85 .9  .95 .99];    % Stimulus coherence
    param.jump_probability          = 0.0025;               % Probability of stimulus direction jump
    param.min_jump_interval_ms      = param.Fs * 100;       % Min interval between direction jumps
    param.feedback_probability    	= 0.0025;               % Probability of reward
    param.min_feedback_interval_ms  = param.Fs * 100;       % Min interval between reward administration
end

% Initialize variables
out.trial                           = param.trial;          % Trial number
out.RDP_direction                   = rand * 2 * pi;      	% Random seed stimulus direction
out.feedback_ts                     = 1;                    % Initialise variable
out.jump_ts                         = 1;                    % Initialise variable
out.RDP_coherence                   = [];                   % Initialise variable
cnt_jump                            = 0;                    % Reset counter
cnt_fb                              = 0;                    % Reset counter
sign_change                         = [-1 1];               % Sign of direction change

% Calculate number of samples
nSamples                            = param.walk_duration_ms / param.Fs;                % Number of samples for random walk
nCohSamples                         = param.coh_duration_ms / param.Fs;                 % Number of samples for coherence block
nCohBlocks                          = param.walk_duration_ms/param.coh_duration_ms;     % Number of coherence block in stimulus cycle
out.RDP_coherence_ts                = [0 nCohSamples.*(1:nCohBlocks-1)] .* param.Fs;    % Timstamps of coherence change [s]
out.RDP_coherence_smple             = [0 nCohSamples.*(1:nCohBlocks-1)] ;               % Samples of coherence change [#]

if param.trial == 0
    out.RDP_coherence              	= fliplr(param.snr_list);   % Ordered coherence values
else
    out.RDP_coherence            	= param.snr_list(randperm(length(param.snr_list))); % Shuffled coherence values
end

% Perform random walk in polar space
for iSample = 2:nSamples %*2 % Generate double the amount of samples to be on the safe side
    
    % Calculate random direction step: uniform distribution * constant (max 1.8deg/sample)
    direction_step                      = (pi/50) * rand * sign_change(randi([1,2])); 
    
    % Update stimulus direction
    out.RDP_direction(iSample)          = out.RDP_direction(iSample-1) + direction_step;
    
    % Random new seed
    if rand < param.jump_probability && (iSample-out.jump_ts(end)) > (param.min_jump_interval_ms/param.Fs)
        cnt_jump                        = cnt_jump+1;
        out.jump_ts(cnt_jump)           = iSample;                          % Frame of feedback presentation
        out.dir_seed(cnt_jump)          = (rand * pi * sign_change(randi([1,2]))) + out.RDP_direction(iSample-1);                    % New random seed
        out.RDP_direction(iSample)  	= out.dir_seed(cnt_jump);
    end
    
    % Random reward
    if rand < param.feedback_probability && (iSample-out.feedback_ts(end)) > (param.min_feedback_interval_ms/param.Fs)
        cnt_fb                      	= cnt_fb+1;
        out.feedback_ts(cnt_fb)     	= iSample;                          % Frame of feedback
    end
end

% Smooth direction values with Gaussian window
out.RDP_direction                       = smoothdata(out.RDP_direction,"gaussian",100);

% Convert to degree
out.RDP_direction                       = mod(out.RDP_direction, 2 * pi); % Wrap around if exceeding 2*pi
out.RDP_direction                       = rad2deg(out.RDP_direction);

end

% figure
% plot(out.RDP_direction, 'k','LineWidth', 1.5) % Plot direction
% scatter(out.feedback_ts, out.RDP_direction(out.feedback_ts), 'ro', 'filled') % Plot reward
% scatter(out.jump_ts, rad2deg(out.dir_seed), 'b^', 'filled') % Plot direction of random seeds
