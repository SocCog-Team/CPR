function out = CPR_create_random_walk_v3(param)

if nargin < 1
    % Set parameters
    param.trial                     = 0;
    param.coh_duration_ms       	= 10000;                                % Duration of coherence block
    param.cycle_duration_ms      	= 60000;                                % Duration of random walk
    param.Fs                        = 1000 / 120;                           % Screen sampling rate [ms]
    param.snr_list                  = [.7 .75 .8 .85 .9  .95 .99];          % Stimulus coherence
    param.feedback_probability    	= 0.0035;                                % Probability of reward
    param.min_feedback_interval_ms  = param.Fs * 100;                       % Min interval between reward administration
end

% Prepare calculations and output
nSamples                            = round(param.cycle_duration_ms / param.Fs);        % Number of samples for random walk (7200 frames, 1min)
nCohSamples                         = param.coh_duration_ms / param.Fs;                 % Number of samples for coherence block
nCohBlocks                          = param.cycle_duration_ms/param.coh_duration_ms;    % Number of coherence block in stimulus cycle
out.RDP_coherence                   = [];                                               % Initialise variable
out.RDP_coherence_ts                = [0 nCohSamples.*(1:nCohBlocks-1)] .* param.Fs;    % Timstamps of coherence change [ms]
out.RDP_coherence_smple             = [0 nCohSamples.*(1:nCohBlocks-1)] ;               % Samples of coherence change [#]
out.feedback_ts                     = 1;                                                % Initialise variable
out.trial                           = param.trial;                                      % Trial number
cnt_fb                              = 0;

%%% This requires SNR inputs to match die duration of the cycle %%%

% Shuffle coherence values (except for first trial)
if param.trial == 0
    out.RDP_coherence              	= fliplr(param.snr_list);                           % Ordered coherence values
else
    tmp = [1 1];
    while sum(diff(tmp) == 0) ~= 0
        tmp                         = param.snr_list(randperm(length(param.snr_list))); % Shuffled coherence values
    end
    out.RDP_coherence               = tmp;
end

% Set parameters for random walk
dt                                  = 1/120;                                % Time step (sec @ 120Hz)
tau_sec                             = 2;                                    % Time constant for the angular velocity decay
sigma_rpm                           = .1;                                   % Standard deviation for the random noise in omega
omega                               = zeros(1, nSamples);                  	% Array to store angular velocities
angle                               = zeros(1, nSamples);                 	% Array to store angles

% Initial seed
omega(1)                            = 0;                % Initial random angular velocity
% omega(1)                            = 2 * pi * (rand - 0.5);                % Initial random angular velocity
angle(1)                            = 2 * pi * rand;                        % Initial random angle between 0 and 2*pi

%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tau:      A larger value keeps the angular velocity persistent longer, 
%           while a smaller value makes it decay faster.
% sigma:    A larger value increases the randomness in the angular velocity.
% dt:       Controls the resolution of the time steps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Random walk in polar space
for i = 2:nSamples%*2 % Double number for safety reasons
    % Update omega [Revolutions per minute] - Gaussian random noise with standard deviation scaled by sqrt(dt)
    omega(i) = omega(i-1) - (omega(i-1) / tau_sec) * dt + sigma_rpm * randn * sqrt(2*dt/tau_sec);
    
    % Update angle(t+dt)
    angle(i) = angle(i-1) + omega(i) * dt;
    
    % Constrain angle within [0, 2*pi]
    angle(i) = mod(angle(i), 2 * pi);
    
    % Random reward
    if rand < param.feedback_probability && (i - out.feedback_ts(end)) > (param.min_feedback_interval_ms/param.Fs)
        cnt_fb                      	= cnt_fb+1;
        out.feedback_ts(cnt_fb)     	= i;                                % Frame of feedback
    end
end

% Convert to degree
out.RDP_direction_rad               	= angle;
out.RDP_direction_deg               	= rad2deg(angle);

% close all
% figure;hold on
% plot(omega(1:7200));
% line([0 7200],[0 0])
% xlabel('time [frames]')
% ylabel('ang. velocity [rpm]')
% title('omega (ang. vel.)')
% 
% figure
% plot(angle(1:7200));
% xlabel('time [frames]')
% ylabel('angle [rad]')
% title('RDP direction')
end

% figure; hold on
% plot(out.RDP_direction_deg, 'k','LineWidth', 1.5) % Plot direction
% scatter(out.feedback_ts, out.RDP_direction_deg(out.feedback_ts), 'ro', 'filled') % Plot reward
