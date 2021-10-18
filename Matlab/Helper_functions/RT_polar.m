function RT_polar(fname, pth, ID)

addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
cd(pth)

%% Import file

var_import = {
    'ML_', ...
    'CTRL_', ...
    'RDP_', ...
    'INFO_', ...
    'TRIAL_'};

d = MW_readFile(fname, 'include',var_import);

%% Extract data & calculate RT

% Get variable index
idx.tOn                     = d.event == 'TRIAL_start';
idx.tEnd                    = d.event == 'TRIAL_end';
idx.trg                     = d.event == 'TRIAL_reactionTrigger';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.ttype                   = d.event == 'TRIAL_type';
idx.res                     = d.event == 'CTRL_start_area_flag';

% Extract timestamps
trl.tOn                     = d.time(idx.tOn);                               	% Trial onset
trl.tEnd                    = d.time(idx.tEnd);                              	% Trial end
c                           = 0;

for iTrl = 1:length(trl.tEnd)                                                   % For all trials...
    
    clear trlIdx
    trlIdx                  = d.time >= trl.tOn(iTrl) & ...
        d.time <= trl.tEnd(iTrl);                         % Build trial index
    outcome                 = getTrialData(d.value, trlIdx, idx.outcome);
    ttype                   = getTrialData(d.value, trlIdx, idx.ttype);
    
    if strcmp(outcome, 'hit')
        c = c+1;
        trg_dir(c)          = str2num(ttype{:}(4:end));
        trg_value{c}      	= getTrialData(d.time, trlIdx, idx.trg);          	% Reaction trigger [stimulus] values
        res_value{c}       	= getTrialData(d.time, trlIdx, idx.res);          	% Reaction trigger [stimulus] values
        rt(c)             	= (res_value{c}(end) - trg_value{c}(1))/1e3;
    end
end

%% Bin RT based on target direction

k                   = cellfun(@size, res_value, 'UniformOutput', false);
ex                  = cellfun(@(x) x(2)>1,k);
rt(ex)              = [];
trg_dir(ex)         = [];
n_step              = 30;
vec                 = 0:n_step:360;

clear dat stddat tdir
for iDir = 1:size(vec,2)-1
    didx            = trg_dir > vec(iDir) & trg_dir <= vec(iDir+1);
    dat(iDir)       = median(rt(didx));
    stddat(iDir)    = std(rt(didx));
    tdir(iDir)      = mean([vec(iDir) vec(iDir+1)]);
end

%% PLOT

f                   = figure;
pp                  = polarscatter(deg2rad(trg_dir), rt);
pp.Marker           = 'o';
pp.SizeData         = 15;
pp.MarkerFaceColor  = [0 0 0];
pp.MarkerEdgeColor  = 'none';
pp.MarkerFaceAlpha  = .5;
ax                  = gca;
ax.RLim             = [0 1000];
ax.FontSize         = 16;
hold on

pl                  = polarplot(deg2rad([tdir tdir(1)]), [dat dat(1)]);
pl.LineWidth        = 3;
pl.Color            = [1 0 0 .75];

pv                  = polarplot(deg2rad([tdir tdir(1)]), [stddat stddat(1)]);
pv.LineWidth        = 3;
pv.Color            = [.5 .5 1 .75];

lg                  = legend('Trial RT','Binned Median','Binned Std');
lg.Position         = [0.03    0.07    0.2473    0.0940];
lg.Box              = 'off';

dest_dir = '/Users/fschneider/Documents/MWorks/Plots';
print(f, [dest_dir '/RT_' ID], '-r300', '-dpng');

end