addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab

pth         = '/Users/fschneider/Desktop/';
fname       = '20250616_nil_CPRdyadic_block1_phy4_fxs.h5';
d           = MW_readH5([pth fname]);
%%

var_import = {
    'INFO_', ...
    'TRIAL_', ...
    'IO_joystickDirection', ...
    'IO_joystickStrength',...
    'IO_fixation_flag',...
    'EYE_x_dva',...
    'EYE_y_dva'};

% Initialise variables
idx                         = [];
cyc                         = [];

% Create variable-specific indices
idx.cOn                     = d.event == 'TRIAL_start';
idx.cEnd                    = d.event == 'TRIAL_end';
idx.fixation              	= d.event == 'IO_fixation_flag';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.eye_x_dva              	= d.event == 'EYE_x_dva';
idx.eye_y_dva             	= d.event == 'EYE_y_dva';

% Get trial timestamps
cyc.cOn                     = d.time(idx.cOn);
cyc.cEnd                    = d.time(idx.cEnd);

nSample                     = 10;
tcnt                        = 0;
break_xy                    = [];

% Stimulus cycle (i.e. trial) loop
for iCyc = 1:length(cyc.cEnd)
    clear eye_x eye_y
    % Trial index
    cycIdx                  = [];
    cycIdx                  = d.time >= cyc.cEnd(iCyc)-1e6 & d.time <= cyc.cEnd(iCyc)+2.5e5;
%     cycIdx                  = d.time >= cyc.cEnd(iCyc)-1e6 & d.time <= cyc.cEnd(iCyc)+0.5e5;
    
    outcome                 = getTrialData(d.value, cycIdx, idx.outcome);
    eye_x                   = getTrialData(d.value, cycIdx, idx.eye_x_dva);                    % RDP_direction timestamps
    eye_y                   = getTrialData(d.value, cycIdx, idx.eye_y_dva);                    % RDP_direction timestamps

    if strcmp(outcome{end}, 'FixationBreak')
        tcnt                = tcnt + 1;     
        break_xy(tcnt,:)    = [mean(eye_x(end)) mean(eye_y(end))];
    end
end

xmax    = 27.5;
ymax    = 15.5;

figure
rectangle('Position', [-xmax -ymax xmax*2 ymax*2], 'FaceColor', [.2 .2 .2]);
hold on
xlim([-xmax*1.1 xmax*1.1]);
ylim([-ymax*1.1 ymax*1.1]);

% Create circle with radius 6 centered at [-9 -2]
theta = linspace(0, 2*pi, 100);
x = 6*cos(theta) - 9;
y = 6*sin(theta) - 2;

% Plot the circle
plot(x, y, 'k-', 'LineWidth', 2);
fill(x, y, 'w', 'FaceAlpha', 0.3);        % Plot filled circle with transparency

% Ensure equal aspect ratio to make circle appear round
axis equal

scatter(0,0, 200, 'w+');
scatter(break_xy(:,1), break_xy(:,2), 50, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [1 0 0]);

