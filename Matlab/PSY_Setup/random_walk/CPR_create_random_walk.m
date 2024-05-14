function out = CPR_create_random_walk(param)

if nargin < 1
    % Set parameters
    param.trial                     = 0;
    param.coh_duration_ms       	= 20000;                % Duration of coherence block
    param.walk_duration_ms      	= 60000;                % Duration of random walk (stimulus cycle)
    param.Fs                        = 1000 / 120;           % Screen sampling rate
    param.snr_list                  = [.2 .4 .6 .8];        % Stimulus coherence
    param.polar_step_size           = 0.05;                 % Step size in polar space
    param.reward_probability        = 0.0025;                % Probability of reward target appearance
    param.feedback_duration_ms      = param.Fs * 6;         % Duration of target presentation
    param.feedback_interval_ms      = param.Fs * 60;        % Interval between target presentations
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
out.RDP_coherence_ts                = [0 nCohSamples.*(1:nCohBlocks-1)] .* param.Fs;    % Timstamps of coherence change [s]
out.RDP_coherence_smple             = [0 nCohSamples.*(1:nCohBlocks-1)] ;               % Samples of coherence change [#]

if param.trial == 0
    out.RDP_coherence              	= fliplr(param.snr_list);   % Ordered coherence values
else
    out.RDP_coherence            	= param.snr_list(randperm(length(param.snr_list))); % Shuffled coherence values
end

% Perform random walk in polar space
for iSample = 2:nSamples*2 % Generate double the amount of samples to be on the safe side    
    % Reset flag
    feedback_flag                   = false;
    
    % Update stimulus direction with random step from noise distribution: mean=0; std=param.polar_step_size
    out.RDP_direction(iSample)  	= out.RDP_direction(iSample-1) + (randn * param.polar_step_size);
      
    % Check for reward target appearance: Roll dice but only accept after certain time interval
    if rand < param.reward_probability && (iSample-out.feedback_ts(end)) > (param.feedback_interval_ms/param.Fs)
        cnt                         = cnt+1;
        feedback_flag              	= true;
        out.feedback_ts(cnt)        = iSample; % Frame of feedback presentation
           
        % Extract direction values between the end of ITI (last target) and current target
        if cnt == 1
            sample_idx                      = 1:out.feedback_ts(cnt);
        else
            sample_idx                      = (out.feedback_ts(cnt-1)+(param.feedback_duration_ms/param.Fs)):out.feedback_ts(cnt);
        end
        
        % Smooth direction values between two targets: Gaussian window - 40 samples wide
        out.RDP_direction(sample_idx)    	= smoothdata(out.RDP_direction(sample_idx),"gaussian",40);
    end
    
    %%% Keep RDP direction stable after target for target presentation time
    if ~feedback_flag && (iSample-out.feedback_ts(end)) < (param.feedback_duration_ms/param.Fs)
        out.RDP_direction(iSample)      = out.RDP_direction(iSample-1); % Fix stimulus direction for duration of feedback
    end
end

% Convert to degree
out.RDP_direction                   = mod(out.RDP_direction, 2 * pi); % Wrap around if exceeding 2*pi
out.RDP_direction                   = rad2deg(out.RDP_direction);

end

% figure; hold on
% plot(out.RDP_direction)
% scatter(out.feedback_ts, out.RDP_direction(out.feedback_ts), 'rx')
