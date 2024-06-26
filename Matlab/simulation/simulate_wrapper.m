addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/simulation/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab/
addpath('/Users/fschneider/Documents/MATLAB/CircStat2012a/')

nRep                    = 10000;
random_subj             = false;
plot_flag               = false;

for iRep = 1:nRep
    reward_coop(iRep,:) = simulate_social_payoff('coop',random_subj, plot_flag);
    reward_comp(iRep,:) = simulate_social_payoff('comp',random_subj, plot_flag);
end

figure
avg_dff = median(mean(reward_coop,2))-median(mean(reward_comp,2));
[p,h] = signrank(mean(reward_coop,2),mean(reward_comp,2));

% Display average reward across players
violinplot([mean(reward_coop,2) mean(reward_comp,2)])

ax                          = gca;
ax.FontSize               	= 16;
ax.YLabel.String           	= 'Total reward in simulated block';
ax.XLabel.String           	= 'Social context';
ax.Title.String             = ['Diff: ' num2str(avg_dff) '; p: ' num2str(p)];
ax.XTick                 	= 1:2;
ax.XTickLabel            	= {'Cooperative','Competitive'};
shg