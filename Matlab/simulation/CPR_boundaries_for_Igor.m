addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/PSY_Setup/random_walk
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/

close all
clear all

nRep = 10000;
for iTrial = 1:nRep
    % Generate CPR stimulus (120Hz, 1 min == 7200 samples)
    STIM                    = CPR_create_random_walk_v3();
    
    % Random joystick response
    js_dir                  = rand(size(STIM.RDP_direction_deg)).*360;
    
    % Computer joystick response (shifted RDP - no noise)
%     AGNT.dir_sigma       	= 0;                                        % Direction sigma [deg]
%     AGNT.str_sigma       	= 0;                                        % Cursor width sigma [%]
%     AGNT.lag              	= 200;                                      % Delay to reference point [samples]
%     AGNT.win              	= 50;                                       % Smoothing window size [samples]
%     AGNT.smooth_kernel     	= 'gaussian';                               % Smoothing kernel
%     AGNT                   	= CPR_create_agent_random_walk(STIM,AGNT);
    
    % Leaky integrator
    stim_leaky_int = leaky_stimulus_integration(STIM);

    % Random response: Calculate circular distance to RDP direction
    js_ang_dev_random      	= rad2deg(circ_dist(deg2rad(js_dir),deg2rad(STIM.RDP_direction_deg)));
    js_acc_random          	= abs(1 - abs(js_ang_dev_random / 180));    % normalised accuracy
    
    % Computer player: Calculate circular distance to RDP direction
    idx                     = 200:7200; % to remove transient from random starting point of computer player
%     js_ang_dev_agnt        	= rad2deg(circ_dist(deg2rad(AGNT.dir_smooth(idx)),deg2rad(STIM.RDP_direction_deg(idx))));
    js_ang_dev_agnt        	= rad2deg(circ_dist(deg2rad(stim_leaky_int(idx)),deg2rad(STIM.RDP_direction_deg(idx))));
    js_acc_agnt            	= abs(1 - abs(js_ang_dev_agnt / 180));
    
    % Average of trial
    js_ang_dev_mean_random(iTrial) 	= mean(abs(js_ang_dev_random));
    js_acc_mean_random(iTrial)     	= mean(js_acc_random);
    js_ang_dev_mean_agnt(iTrial)    = mean(abs(js_ang_dev_agnt));
    js_acc_mean_agnt(iTrial)    	= mean(js_acc_agnt);
end

%% PLOT
nBin = 20;
figure
subplot(2,2,1)
histogram(js_ang_dev_mean_random,nBin)
title('Ang. error: Random')
xlabel('Avg. angular distance [deg]')
ylabel('Trials [#]')

subplot(2,2,3)
histogram(js_ang_dev_mean_agnt,nBin)
title('Ang. error: LeakyInt')
xlabel('Avg. angular distance [deg]')
ylabel('Trials [#]')

subplot(2,2,2)
histogram(js_acc_mean_random,nBin)
title('Accuracy: Random')
xlabel('Avg. accuracy [norm]')
ylabel('Trials [#]')

subplot(2,2,4)
histogram(js_acc_mean_agnt,nBin)
title('Accuracy: LeakyInt')
xlabel('Avg. accuracy [norm]')
ylabel('Trials [#]')

%%
figure
plot(STIM.time_vector,STIM.RDP_direction_deg)
hold on
plot(STIM.time_vector,stim_leaky_int)
xlabel('time [s]')
ylabel('Direction [deg]')
set(gca,'Fontsize',16)
legend('RDP','Leaky integrator')
%%
%%% FUNCTIONS
function theta_filt_deg = leaky_stimulus_integration(in)
% Convert to Cartesian coordinates
x            	= cos(in.RDP_direction_rad);
y             	= sin(in.RDP_direction_rad);

% Parameters for leaky integrator
k             	= 0.2;                  % decay constant
t            	= in.time_vector;       % time vector (must match direction trace)
d             	= 0.5;              	% delay (set as needed)

% Apply your leakyIntegrator to both components
x_filt          = leakyIntegrator(x, k, t, d);
y_filt          = leakyIntegrator(y, k, t, d);

% Convert back to angle (in degrees)
theta_filt_rad  = atan2(y_filt, x_filt);
theta_filt_deg  = mod(rad2deg(theta_filt_rad), 360);

end
