tbl1 = load('/Users/fschneider/ownCloud/CPR_data/TestFiles/20210624_nes_cpr_tbl.mat');
tbl2 = load('/Users/fschneider/ownCloud/CPR_data/TestFiles/20210624_sut_cpr_tbl.mat');
clear t 

tbl1 = tbl1.t;
tbl2 = tbl2.t;
trl = 1;

% Check number of samples
tmp = cellfun(@size, tbl2.trl_js_ts, 'UniformOutput', false);
idx = cell2mat(cellfun(@(x) x(2)>1, tmp, 'UniformOutput', false));

for iTrl = find(idx)
    for iSmple = 1:size(tbl2.trl_js_ts{iTrl},2)
    end
end

figure
plot(tbl1.trl_js_ts{trl},tbl1.trl_js_dir{trl})
hold on
plot(tbl2.trl_js_ts{trl},tbl2.trl_js_dir{trl})

% Determine similarity between joystick signals:
% Accuracy
% Displacement
% Circular correlation
% Assessment of psychometric function

% Determine if Leader-Follower present:
% Crosscorrelation -> Lag distribution
% Granger causality
% Shifted, Kalman-filtered signal to predict timecourse of other joystick

