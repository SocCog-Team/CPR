%% Sample difference between expected and recorded state duration and 

expected = round(t.ss_dur/10);
recorded = cellfun(@length,t.js_dir);

idx = abs(expected-recorded) > 1;
tt = t(idx,:);

histogram([round(t.ss_dur/10)-cellfun(@length,t.js_dir)],50)
% histogram([round(tt.ss_dur/10)-cellfun(@length,tt.js_dir)],50)
xlabel('Sample Difference [#]')
set(gca,'FontSize',16)

%% Joystick sample difference

histogram(diff(d.time(d.event=='IO_joystickDirection')))
xlim([9900 10100])
xlabel('Inter-sample interval [us]')
set(gca,'FontSize',16)