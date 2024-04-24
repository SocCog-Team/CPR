function out = simulate_social_payoff(social_context,random_subj, plot_flag)

% Parameters
num_steps           = 1e4;                  % Number of steps in the random walk
polar_step_size     = 0.1;                  % Step size in polar space
lag                 = 50;                   % Lag in tracking the stimulus direction
reward_probability  = 0.01;                 % Probability of reward target appearance
nscale              = pi/10;                % Noise scaling

% Initialize variables
stimulus_direction  = rand * 2 * pi;        % Initial stimulus direction
subject1_direction  = stimulus_direction;   % Initial subject1 direction
subject2_direction  = stimulus_direction;   % Initial subject2 direction
reward_targets      = [];

% Perform random walk in polar space
for step = 2:num_steps
    % Update stimulus direction with random step
    stimulus_direction(step)            = stimulus_direction(step-1) + randn * polar_step_size;
    stimulus_direction(step)         	= mod(stimulus_direction(step), 2 * pi); % Wrap around if exceeding 2*pi
    
    if random_subj
        subject1_direction(step)        = pi + nscale*randn;
        subject2_direction(step)        = pi + nscale*randn;
    else
        if step > lag
            % Subject direction lags behind stimulus direction
            subject1_direction(step)    = stimulus_direction(step-lag) + (nscale*randn);
            subject2_direction(step)    = stimulus_direction(step-lag/2) + (nscale*randn);
        else
            subject1_direction(step)    = pi + nscale*randn;
            subject2_direction(step)    = pi + nscale*randn;
        end
    end
    
    % Check for reward target appearance
    if rand < reward_probability
        reward_targets  = [reward_targets; step, stimulus_direction(step), subject1_direction(step), subject2_direction(step)];
    end
end

% Initialize variables for rewards
rewards = zeros(size(reward_targets, 1), 3);

% Calculate rewards based on tracking accuracy
for i = 1:size(reward_targets, 1)
    % Calculate tracking errors
    tracking_error_subject1     = circ_dist(reward_targets(i, 2), reward_targets(i, 3));
    tracking_error_subject2     = circ_dist(reward_targets(i, 2), reward_targets(i, 4));
    
    % Calculate accuracies
    accuracy_subject1       = abs(1 - abs(rad2deg(tracking_error_subject1)) / 180);
    accuracy_subject2       = abs(1 - abs(rad2deg(tracking_error_subject2)) / 180);
    
    % Determine rewards based on behavior - dyadic hit
    if accuracy_subject1 > 0.75 && accuracy_subject2 > 0.75   % Dyadic hit
        if strcmp(social_context, 'coop')
            % Cooperative context: average accuracy between players if both hit the target
            reward1         = sum([accuracy_subject1 accuracy_subject2]) / 2;
            reward2         = sum([accuracy_subject1 accuracy_subject2]) / 2;
        elseif strcmp(social_context, 'comp')
            % Competitive context: winner takes all if both hit the target
            if accuracy_subject1 > accuracy_subject2
                reward1     = sum([accuracy_subject1 accuracy_subject2]);
                reward2     = 0;
            else
                reward1     = 0;
                reward2     = sum([accuracy_subject1 accuracy_subject2]);
            end
        end
    else
        reward1             = 0;
        reward2             = 0;
    end
    rewards(i, :)           = [reward_targets(i, 1), reward1, reward2];
end

% Sum of reward for each player
out = sum(rewards(:,2:3)); 

if plot_flag
    % Plot results
    figure;
    hold on;
    plot(stimulus_direction, 'k', 'LineWidth', 2);
    plot(subject1_direction, 'r', 'LineWidth', 1);
    plot(subject2_direction, 'g', 'LineWidth', 1);
    scatter(rewards(:, 1), ones(length(rewards(:, 1)),1), 100, 'k', 'Marker', '*');
    xlabel('Samples');
    ylabel('Direction (radians)');
    title('Stimulus and Subject Tracking with Reward Targets');
    legend('Stimulus Direction', 'Subject1 Direction', 'Subject2 Direction', 'Reward Target');
    set(gca,'fontsize',16)
    hold off;
    
    figure
    pol_diff = circ_dist(stimulus_direction(1:end-1),stimulus_direction(2:end));
    h = histogram(rad2deg(pol_diff));
    h.FaceColor = [.5 .5 .5];
    h.EdgeColor = [.5 .5 .5];
    h.FaceAlpha = 1;
    xlabel('Sample difference [deg]');
    ylabel('#');
    set(gca,'fontsize',16)
    hold off;
end