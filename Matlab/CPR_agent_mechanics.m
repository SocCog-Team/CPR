%% Generate toy data

% Time specifications:
Fs              = 1000;                 	% samples per second
dt              = 1/Fs;                 	% seconds per sample
tEnd            = 1;                       	% seconds
t               = (0:dt:tEnd-dt)';          % seconds
frq             = 3;                        % hertz
y               = cos(2*pi*frq*t);          % sine wave:
yy              = sin(2*pi*frq*t);          % cosine wave:

% Add slope to data
slp             = (1:length(y))/500;        % slope
dat_rdp       	= y+slp';                   % RDP data
dat_js         	= yy+slp';                  % JS data
dat_coh         = [ones(1,length(dat_js)/4)*.4,...
    ones(1,length(dat_js)/4)*.8 ...
    ones(1,length(dat_js)/4)*.6 ...
    ones(1,length(dat_js)/4)*.8];           % RDP Coherence

%% Model agent response

% Generate normally distributed noise
r               = random('Normal',0,std(dat_rdp),[1 length(y)]);
r2              = random('Normal',0,std(dat_rdp),[1 length(y)]);
r3              = random('Normal',0,std(dat_rdp),[1 length(y)]);
r4              = random('Normal',0,std(dat_rdp),[1 length(y)])/3;

% Add noise to data ans smooth with moving window
win             = 50;
dat_js_pad      = [zeros(win,1); dat_js(1:end-win)];
agent_rdp      	= smoothdata(dat_rdp+r','gaussian',win);
agent_rdp_js    = smoothdata(mean([dat_rdp,dat_js],2)+r2','gaussian',win);
agent_js        = smoothdata(dat_js_pad+r3','gaussian',win);
agent_str       = smoothdata(dat_coh+r4,'gaussian',win);

%% Plot
close all

lw              = 2;
alp             = .75;
dest_dir        = '/Users/fschneider/Desktop';

f               = figure;
p_rdp           = plot(t,dat_rdp,'LineWidth', lw, 'LineStyle','--', 'Color',[0 0 0]); hold on
p_js            = plot(t,dat_js,'LineWidth', lw, 'LineStyle','--', 'Color',[1 0 0]);
p_agt_rdp       = plot(t,agent_rdp,'LineWidth', lw, 'LineStyle','-', 'Color',[.5 .5 .5 alp]);
p_agt_js        = plot(t,agent_rdp_js,'LineWidth', lw, 'LineStyle','-', 'Color',[1 .5 .5 alp]);

xlabel('Time')
ylabel('Direction signal')
set(gca,'FontSize', 20)
legend('Signal: RDP','Signal: JS', 'Agent: RDP','Agent: RDP+JS','Interpreter','none','Location','southeast','FontSize',12)
print(f, [dest_dir '/Agent'], '-r300', '-dpng');

%% Animate

close all
lw                  = 3;
pl                  = polarplot([dat_rdp(1) dat_rdp(1)],[0 dat_coh(1)]);

for i = 2:length(dat_rdp)
    cla
    pl              = polarplot([dat_rdp(i) dat_rdp(i)],[0 dat_coh(i)]);
    pl.LineWidth    = lw*2;
    pl.Color        = [0 0 0];
    hold on   
    
    pl              = polarplot([dat_js(i) dat_js(i)],[0 agent_str(i)]);
    pl.LineWidth    = lw;
    pl.Color        = [1 0 0];
    
%     pl              = polarplot([agent_rdp_js(i) agent_rdp_js(i)],[0 agent_str(i)]);
%     pl.LineWidth    = lw;
%     pl.Color        = [0 1 1];
    
    pl              = polarplot([agent_js(i) agent_js(i)],[0 agent_str(i)]);
    pl.LineWidth    = lw;
    pl.Color        = [0 1 1];
    
%     pl             = polarplot([agent_rdp(i) agent_rdp(i)],[0 agent_str(i)]);
%     pl.LineWidth   = lw;
%     pl.Color       = [0 1 1];
        
    rlim([0 1])
    pause(0.0005)

end