plot_RF_mapping_schematic()

function plot_RF_mapping_schematic()
% PLOT_RF_MAPPING_SCHEMATIC  Draw a schematic of the stimulus screen as used
%                            during RF mapping, showing:
%
%   - All 90 RDP stimulus positions (grid: x = -24:3:3 dva, y = -12:3:12 dva)
%   - The 2 dva radius random dot pattern footprint at each position
%   - The fixation cross at [0, 0]
%   - The fixation window circle  (r = 1.5 dva)
%   - The cursor aperture circle  (r = 1.75 dva)
%
% No input required — all parameters are hard-coded to match RF_mapping.


% =========================================================================
% 1.  STIMULUS GRID PARAMETERS  (must match RF_mapping exactly)
% =========================================================================

x_position  = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];   % dva, 1×10
y_position  = [-12  -9  -6  -3   0  3  6  9 12];     % dva, 1×9 (ascending)

rdp_radius  = 2.0;    % dva  — RDP aperture radius
cursor_r    = 1.75;   % dva  — cursor aperture radius
fixwin_r    = 1.50;   % dva  — fixation window radius
fix_pos     = [0, 0]; % dva  — fixation spot location


% =========================================================================
% 2.  FIGURE & AXES  (publication-ready)
% =========================================================================

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 14]);

ax = axes('Parent', fig);
set(ax, ...
    'Color',     [0.97 0.97 0.97], ...
    'Box',       'on', ...
    'TickDir',   'out', ...
    'LineWidth',  0.8, ...
    'FontName',  'Arial', ...
    'FontSize',   11, ...
    'DataAspectRatio', [1 1 1], ...   % equal dva scaling on both axes
    'Layer',     'top');
hold(ax, 'on');


% =========================================================================
% 3.  DRAW RDP CIRCLES AT EVERY GRID POSITION
% =========================================================================
% Each circle represents the 2 dva-radius random dot pattern aperture
% centred on a stimulus grid location.  We draw them first so all other
% overlays sit on top.

angle_vec = linspace(0, 2*pi, 361);   % fine angle vector for smooth circles

rdp_face  = [0.75 0.85 0.95];   % light blue fill
rdp_edge  = [0.45 0.60 0.80];   % medium blue edge

for ix = 1:numel(x_position)
    for iy = 1:numel(y_position)
        cx = x_position(ix);
        cy = y_position(iy);

        % Filled circle for the RDP aperture
        fill(ax, ...
             cx + rdp_radius .* cos(angle_vec), ...
             cy + rdp_radius .* sin(angle_vec), ...
             rdp_face, ...
             'EdgeColor', rdp_edge, ...
             'LineWidth', 0.5, ...
             'FaceAlpha', 0.40);

        % Centre dot to mark the grid position precisely
        plot(ax, cx, cy, '.', ...
             'Color',      rdp_edge, ...
             'MarkerSize', 4);
    end
end


% =========================================================================
% 4.  FIXATION WINDOW  (dashed dark-grey circle, r = 1.5 dva)
% =========================================================================

plot(ax, ...
     fix_pos(1) + fixwin_r .* cos(angle_vec), ...
     fix_pos(2) + fixwin_r .* sin(angle_vec), ...
     '--', ...
     'Color',     [0.25 0.25 0.25], ...
     'LineWidth',  1.4, ...
     'DisplayName', sprintf('Fixation window  (r = %.2g dva)', fixwin_r));


% =========================================================================
% 5.  CURSOR APERTURE  (solid red circle, r = 1.75 dva)
% =========================================================================

plot(ax, ...
     fix_pos(1) + cursor_r .* cos(angle_vec), ...
     fix_pos(2) + cursor_r .* sin(angle_vec), ...
     '-', ...
     'Color',     [0.85 0.10 0.10], ...
     'LineWidth',  1.4, ...
     'DisplayName', sprintf('Cursor aperture  (r = %.2g dva)', cursor_r));


% =========================================================================
% 6.  FIXATION CROSS  (black + marker at [0,0])
% =========================================================================
% Drawn last so it is never occluded by any circle.

cross_arm = 0.8;   % dva — half-length of the cross arms

plot(ax, [fix_pos(1)-cross_arm, fix_pos(1)+cross_arm], ...
         [fix_pos(2),           fix_pos(2)], ...
         'k-', 'LineWidth', 2.0);
plot(ax, [fix_pos(1), fix_pos(1)], ...
         [fix_pos(2)-cross_arm, fix_pos(2)+cross_arm], ...
         'k-', 'LineWidth', 2.0);


% =========================================================================
% 7.  AXIS LIMITS, LABELS & LEGEND
% =========================================================================

margin = rdp_radius + 1;   % dva of white space beyond the outermost circles
xlim(ax, [min(x_position) - margin,  max(x_position) + margin]);
ylim(ax, [min(y_position) - margin,  max(y_position) + margin]);

% Tick at every grid line
ax.XTick = x_position;
ax.YTick = y_position;

xlabel(ax, 'Horizontal position  (dva)', 'FontSize', 13, 'FontName', 'Arial');
ylabel(ax, 'Vertical position  (dva)',   'FontSize', 13, 'FontName', 'Arial');

% Append the RDP aperture manually to the legend (fill patches aren't
% auto-captured by legend when drawn in a loop)
h_rdp = fill(ax, NaN, NaN, rdp_face, ...
             'EdgeColor', rdp_edge, 'LineWidth', 0.5, 'FaceAlpha', 0.40, ...
             'DisplayName', sprintf('RDP aperture  (r = %.g dva)', rdp_radius));

% Draw grid lines at every stimulus position for spatial reference
ax.XGrid      = 'on';
ax.YGrid      = 'on';
ax.GridColor  = [0.8 0.8 0.8];
ax.GridAlpha  = 0.6;
ax.GridLineStyle = ':';

leg = legend(ax, 'Location', 'northeastoutside');
leg.FontSize  = 10;
leg.FontName  = 'Arial';
leg.Box       = 'off';

title(ax, 'RF mapping schematic', ...
      'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'normal');

hold(ax, 'off');
drawnow;

end