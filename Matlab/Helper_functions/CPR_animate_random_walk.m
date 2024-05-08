close all

% Create the data
out                         = CPR_create_random_walk();
theta                       = deg2rad(out.RDP_direction);
out.RDP_coherence_smple(1) 	= 1;
out.RDP_coherence_smple    	= [out.RDP_coherence_smple 7200];
r                           = [];
col                         = jet(7200);
lw                          = 1.5;
fs                          = 16;
frameRate                   = 120;

for iCoh = 2:length(out.RDP_coherence_smple)
    indx                    = out.RDP_coherence_smple(iCoh-1):out.RDP_coherence_smple(iCoh);
    r(indx)                 = repmat(out.RDP_coherence(iCoh),1,length(indx));
end

% Create a VideoWriter object
writerObj                   = VideoWriter('polar_animation.mp4', 'MPEG-4');
writerObj.Quality        	= 50; % Adjust quality (default is 75)
writerObj.FrameRate         = frameRate; % Set the frame rate [Hz]
open(writerObj);

% Create a figure
figure;
ax = gca;

% Loop through each frame
for i = 1:7200
    % Clear the previous plot
    cla(ax);
    
    % Plot stimulus history up to the current frame
    polarplot(theta(1:i), r(1:i), 'k'); % Plot the entire curve
    hold on;
    
    % Plot current position
    polarplot(theta(i), r(i), 'ro', 'MarkerSize', 10, 'Linewidth', lw); 
    
    % Plot a line from the origin to the current position
    polarplot([0 theta(i)], [0 r(i)], 'r', 'Linewidth', lw); % Plot the red line
    
    % Customize axes
    ax                      = gca;
    ax.ThetaZeroLocation    = 'top';
    ax.ThetaDir             = 'clockwise';
    ax.RLim                 = [0 1];
    ax.FontSize             = fs;
    tx                      = text(1,1.3,['Frame: ' num2str(i)]);
    
    % Use drawnow to update the figure
    drawnow;
        
    % Write the current frame to the video file
    writeVideo(writerObj, getframe(gcf));
    
    % Pause for a short duration to create the animation effect
    pause(1/frameRate)
end

% Close the video file
close(writerObj);