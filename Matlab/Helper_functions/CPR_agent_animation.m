% ANIMATION

close all
lw                  = 3;
for i = 1:length(AGNT.dir_clean)
    cla
    pl              = polarplot([deg2rad(AGNT.dir_clean(i)) deg2rad(AGNT.dir_clean(i))],[0 AGNT.str_clean(i)]);
    pl.LineWidth    = lw;
    pl.Color        = [0 0 0];
    hold on
    
    pl              = polarplot([deg2rad(AGNT.dir_smooth(i)) deg2rad(AGNT.dir_smooth(i))],[0 AGNT.str_smooth(i)]);
    pl.LineWidth    = lw;
    pl.Color        = [1 0 0 .75];
    
    rlim([0 1])
    pause(0.0005)
end



% negative weights
% monkey needs to fill bar to get reward
