% addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
% cd /Users/fschneider/Documents/MWorks/Data
% fname = 'fschneider-PHY_CPR_Training_20210820-20210825-185052.mwk2';
% d = MW_readFile(fname, 'include', {'IO_joystickDirection','IO_joystickStrength','#stimDisplay'}, 'dotPositions');

% addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
% cd /Users/fschneider/Desktop
% fname = 'fschneider-CPR_dyatic_fxs_20210901-20210901-092758.mwk2';
% d = MW_readFile(fname, 'include', {'IO_joystickDirection','IO_joystickStrength','#stimDisplay'}, 'dotPositions');

idx             = [];
idx.dot         = d.event == 'STIM_RDP_dot_positions';
idx.rdp_dir     = d.event == 'STIM_RDP_direction';
idx.rdp_coh     = d.event == 'STIM_RDP_coherence';
idx.js_dirA     = d.event == 'IO_joystickDirection_A';
idx.js_strA     = d.event == 'IO_joystickStrength_A';
idx.js_dirB     = d.event == 'IO_joystickDirection_B';
idx.js_strB     = d.event == 'IO_joystickStrength_B';
idx.tx          = d.event == 'STIM_target_x_position';
idx.ty          = d.event == 'STIM_target_y_position';

dp              = d.value(idx.dot);
dp_ts           = d.time(idx.dot);
js_dirA         = d.value(idx.js_dirA);
js_strA         = d.value(idx.js_strA);
js_dirB         = d.value(idx.js_dirB);
js_strB         = d.value(idx.js_strB);
js_ts           = d.time(idx.js_dirA);
rdp_dir         = d.value(idx.rdp_dir);
rdp_ts          = d.time(idx.rdp_dir);
rdp_coh         = d.value(idx.rdp_coh);
rdp_coh_ts      = d.time(idx.rdp_coh);
trg_x          	= d.value(idx.tx);
trg_y         	= d.value(idx.ty);
trg_ts        	= d.time(idx.tx);

idx.rm_samples          = cell2mat(cellfun(@(x) isa(x,'int64'),js_dirA,'UniformOutput',false));
js_dirA(idx.rm_samples) = [];
js_strA(idx.rm_samples) = [];
js_dirB(idx.rm_samples) = [];
js_strB(idx.rm_samples) = [];
js_ts(idx.rm_samples)   = [];

js_strA(cellfun(@(x) x>1, js_strA)) = {1};
js_strB(cellfun(@(x) x>1, js_strB)) = {1};

f               = figure;
f.Units         = 'centimeters';
f.Position(3:4) = [15 15];
ax              = axes('Units', 'centimeters', 'Position', [1 1 13 13]); hold on
ax.XLim         = [-10 10];
ax.YLim         = [-10 10];
ax.XLabel.String= '[dva]';
ax.YLabel.String= '[dva]';
ax.FontSize     = 20;

th              = 0:pi/50:2*pi;
x_circle        = .5 * cos(th);
y_circle        = .5 * sin(th);
r               = 8.5;
r2              = 9.5;

for iFrame = 1000:4000
    
    tic
    cla
    js_idx    	= js_ts >= dp_ts(iFrame) & js_ts <= dp_ts(iFrame+1);
    coh_idx    	= rdp_coh_ts >= dp_ts(iFrame) & rdp_coh_ts <= dp_ts(iFrame+1);
    rdp_idx    	= rdp_ts >= dp_ts(iFrame) & rdp_ts <= dp_ts(iFrame+1);
    trg_idx    	= trg_ts >= dp_ts(iFrame) & trg_ts <= dp_ts(iFrame+1);

    % Joystick position
    if sum(js_idx) > 0
        smpl  	= find(js_idx,1, 'last');
         ang   	= mod((90-js_dirA{smpl}) + 360,360);
         [x,y] 	= pol2cart(deg2rad(ang),js_strA{smpl}*8);
         wdth   = deg2rad(180 - (180 * js_strA{smpl}));
                  
         ang2   = mod((90-js_dirB{smpl}) + 360,360);
         [x2,y2]= pol2cart(deg2rad(ang2),js_strB{smpl}*8);
         wdth2  = deg2rad(180 - (180 * js_strB{smpl}));
    end
    
    if (180 - (180 * js_strA{smpl})) < 13
        wdth	= deg2rad(13);
    end
    
    if (180 - (180 * js_strB{smpl})) < 13
        wdth2	= deg2rad(13);
    end
    
    th_str      = linspace(wdth/2, -wdth/2, 100) + deg2rad(ang);
    xxx         = r * cos(th_str);
    yyy         = r * sin(th_str);
    th_str2      = linspace(wdth2/2, -wdth2/2, 100) + deg2rad(ang2);
    xxx2         = r2 * cos(th_str2);
    yyy2         = r2 * sin(th_str2);
    
    % Stimulus coherence
    if sum(coh_idx) > 0
        smple_c	= find(coh_idx,1, 'last');
        coh   	= rdp_coh{smple_c};
    end
    
    % Stimulus direction
    if sum(rdp_idx) > 0 
        smple  	= find(rdp_idx,1, 'last');
        ang     = mod((90-rdp_dir{smple}) + 360,360);
        [xx,yy]	= pol2cart(deg2rad(ang),coh*8);
    end
    
    xIdx       	= logical(mod([1:size(dp{iFrame},2)],2));
    xpos        = dp{iFrame}(xIdx);
    ypos        = dp{iFrame}(~xIdx);
    
    p_dots      = plot(xpos,ypos,'.w','MarkerSize', 2); set(gca,'Color',[0 0 0]); hold on;
    p_crcl      = fill(x_circle, y_circle, 'k');
    p_arc       = plot(xxx, yyy, 'LineWidth', 15, 'Color', [1 0 0]);
    p_arc2      = plot(xxx2, yyy2, 'LineWidth', 15, 'Color', [0 0 1]);
    p_rdp      	= quiver(0,0,xx,yy,0, 'LineWidth', 3, 'Color', [.7 .7 .7]);
    p_js      	= quiver(0,0,x,y,0, 'LineWidth', 3, 'Color', [1 0 0]);
    p_js2      	= quiver(0,0,x2,y2,0, 'LineWidth', 3, 'Color', [0 0 1]);
    
%     p_ln        = plot([0 x],[0 y],'k','LineWidth',2);
%     p_js        = plot(x,y,'r.','MarkerSize',40);
%     p_rdp       = plot(xx,yy,'b.','MarkerSize',40);

    % Target appearance
    if sum(trg_idx) > 0
        smple_t	= find(trg_idx,1, 'last');
        tx      = trg_x{smple_t};
        ty      = trg_y{smple_t};
        
        tpos_x 	= mean([r r2]) * cos(cart2pol(tx,ty));
        tpos_y 	= mean([r r2]) * sin(cart2pol(tx,ty));
        p_trg  	= plot(tpos_x,tpos_y,'w.','MarkerSize',150);
    end
        
    t_loop      = toc;
    pause(1/10000)
end
