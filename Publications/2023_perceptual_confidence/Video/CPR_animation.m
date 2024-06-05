load('/Volumes/T7_Shield/CPR_psychophysics/AaA/summary/20230117_aaa_CPRsolo_block2_tbl.mat')
close all

data_idx    = t.cyc_no == 11;
dots        = t.rdp_dot(data_idx);
coh         = t.rdp_coh(data_idx);
rdp_dir     = t.rdp_dir(data_idx);
rdp_coh     = t.rdp_coh(data_idx);
js_dir      = t.js_dir(data_idx);
js_ecc      = t.js_ecc(data_idx);
dt_x        = [];
dt_y        = [];

% Create a VideoWriter object
writerObj                   = VideoWriter('CPR_animation_paper.mp4', 'MPEG-4');
writerObj.Quality        	= 20; % Adjust quality (default is 75)
writerObj.FrameRate         = 100; % Set the frame rate [Hz]
open(writerObj);

% Create a figure
f                           = figure('units','centimeters','position',[0 0 15 15]);
ax                          = gca;

for iState = 1:length(dots)
    for iFrame = 1:length(dots{iState})
         
        % Clear axes
        cla
        
        % Dot index
        xIdx        = 1:2:length(dots{iState}{iFrame});
        yIdx        = 2:2:length(dots{iState}{iFrame});
         
        % Plot dots
        sc          = scatter(dots{iState}{iFrame}(xIdx),dots{iState}{iFrame}(yIdx),'filled','MarkerFaceColor','w');
        sc.SizeData = 10;
        hold on
        
        % Plot nominal stimulus direction - vector scaled by coherence
        [x,y]       = pol2cart(pi/2-deg2rad(rdp_dir(iState)),rdp_coh(iState));
        ln         	= line([0 x*8],[0 y*8],'Color','w','LineWidth',2);
        
        % Joystick direction and cursor width
        js_rad      = pi/2-deg2rad(js_dir{iState}(iFrame)); 
        width       = deg2rad((180 - (180 * js_ecc{iState}(iFrame))));
      
        % Plot patch to indicate angular width of cursor
        th          = linspace(js_rad-(width/2),js_rad+(width/2), 100);
        xx          = 8*cos(th);
        yy          = 8*sin(th);
        pt          = patch([0 xx],[0 yy],'r','EdgeColor','none', 'FaceAlpha',.3);

        % Plot line to indicate joystick eccentricity
        [x,y]       = pol2cart(js_rad,js_ecc{iState}(iFrame));
        ln          = line([0 x*8],[0 y*8],'Color','r','LineWidth',2);
        
        % Calculate and show frame-wise score
        accuracy    = abs(1 - abs(rdp_dir(iState) - js_dir{iState}(iFrame)) / 180);
        score       = accuracy * js_ecc{iState}(iFrame); 
        tx          = text(-7.8,-7.5, ['Score: ' num2str(round(score,2))]);
        tx.Color    = [1 1 1];
        tx          = text(5.5,-7.5, ['Coh: ' num2str(round(rdp_coh(iState),2))]);
        tx.Color    = [1 1 1];
        tx          = text(-7.8,-7, ['Acc: ' num2str(round(accuracy,2))]);
        tx.Color    = [1 1 1];
        tx          = text(-7.8,-6.5, ['Ecc: ' num2str(round(js_ecc{iState}(iFrame),2))]);
        tx.Color    = [1 1 1];
        
        % Adjust axes
        ax          = gca;
        ax.Color    = [0 0 0];
        xlim([-8 8])
        ylim([-8 8])
        ax.XTick = [];
        ax.YTick = [];
        axis square
        
        % Use drawnow to update the figure
        drawnow;
        
        % Write the current frame to the video file
        writeVideo(writerObj, getframe(gcf));
%         pause(1/120)
    end
end


% Close the video file
close(writerObj);
