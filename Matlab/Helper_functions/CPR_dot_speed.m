addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab

% fname = '20211119_cla_cpr.mwk2';
% d = MW_readFile(fname, 'include', {'#stimDisplay'}, 'dotPositions');
% load('/Volumes/userinterchange/Felix Schneider/forKate/dots_only.mat')

idx                 = [];
idx.dot             = d.event == 'STIM_RDP_dot_positions';
idx.dot_speed       = d.event == 'STIM_RDP_speed';
idx.dot_lifetime    = d.event == 'STIM_RDP_lifetime';

speed               = unique(cell2mat(d.value(idx.dot_speed)));              	% Predefined dot speed
lifetime            = unique(cell2mat(d.value(idx.dot_lifetime)));            	% Predefined dot lifetime
dp                  = d.value(idx.dot);                                         % Actual dot position
dt                  = d.time(idx.dot);                                          % Actual dot time

% Determine if actual dot speed matches predefined setting
pool                = [];
xIdx                = logical(mod([1:size(dp{1},2)],2));                        % Index for x/y position

for iFrme = 2:1000%size(dp,2)
    clear xpos_last ypos_last xpos ypos dist                                    % Clear temporary variables
    
    xpos_last     	= dp{iFrme-1}(xIdx);                                    	% Last frame: x-position
    ypos_last      	= dp{iFrme-1}(~xIdx);                                      	% Last frame: y-position
    xpos        	= dp{iFrme}(xIdx);                                        	% This frame: x-position
    ypos          	= dp{iFrme}(~xIdx);                                         % This frame: y-position

    for iDot = 1:length(ypos)
        dist(iDot)	= pdist([xpos_last(iDot),ypos_last(iDot);xpos(iDot),ypos(iDot)],'euclidean'); % distance/vector length between points
    end
    
    pool            = [pool, dist];
end

%% VISUALISE

% Plot frame-wise speed of all dots
figure
histogram(pool)
title(['Median: ' num2str(median(pool)) ' dva']);
xlabel('Dot speed [dva]')
ylabel('No. dots [#]')
set(gca,'FontSize',20)

figure
ofs = .0005;
lower_lim = median(pool) - ofs;
upper_lim = median(pool) + ofs;
new = pool(pool > lower_lim & pool < upper_lim); % Exclude dots that reappear elsewhere after lifetime expired
histogram(new)
xlim([lower_lim upper_lim])
xlabel('Dot speed [dva]')
ylabel('No. dots [#]')
set(gca,'FontSize',20)

% Plot timestamp difference
figure
df = diff(dt)/1e3;
histogram(df)
title(['Median: ' num2str(median(df)) ' ms']);
xlabel('Inter-frame interval [ms]')
ylabel('No. frames [#]')
set(gca,'FontSize',20)

figure
df = diff(dt)/1e3;
ofs = 5;
lower_lim = median(df) - ofs;
upper_lim = median(df) + ofs;
dff = df(df > lower_lim & df < upper_lim); % Exclude timestamp differences that are enforced by external task demands (i.e. inter-trial interval)
histogram(dff)
xlim([lower_lim upper_lim])
xlabel('Inter-frame interval [ms]')
ylabel('No. frames [#]')
set(gca,'FontSize',20)