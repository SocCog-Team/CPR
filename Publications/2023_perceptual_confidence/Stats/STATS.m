%%% STATS %%%

% Import data
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/dyad_human_human_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/dyad_human_human_performance.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/dyad_pairwise_correlation.mat')
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/dyad_pairwise_performance.mat')

%% Scores increase with motion coherence
tbl_scr             = generate_table(solo_perf,'trg_score','perf');
[lmeStat,comp]    	= LME_coherence_model(tbl_scr);

%% Hit rate increases with motion coherence
tbl_hir             = generate_table(solo_perf,'hir','hir');
stats_hir           = LME_coherence_model(tbl_hir);

%% Average response lag
for iSubj = 1:length(solo_cr)
    mlag(iSubj)     = median(solo_cr{iSubj}.lag); % coherence pooled
end

avg_lag_population = round(mean(mlag));
std_lag_population = round(std(mlag));

%% Faster responses for high motion coherence
tbl_lag             = generate_table(solo_cr,'lag','cr');
stats_lag           = LME_coherence_model(tbl_lag);

%% Accuracy increases with motion coherence
tbl_acc             = generate_table(solo_perf,'trg_score','perf');
stats_acc           = LME_coherence_model(tbl_acc);

%% Eccentricity increases with motion coherence
tbl_ecc             = generate_table(solo_perf,'ecc','perf');
stats_ecc           = LME_coherence_model(tbl_ecc);

% Lag solo vs dyadic different

% Avg lag difference

% Response accuracy solo vs dyadic

% Response eccentricity solo vs dyadic

% Correlation: score difference - AUROC difference
%%% see FIG2_v3 line 485
%%% [r,p] = corrcoef(abs(dff),abs(effect_dff));

% N subjects with significantly different eccentricity
%%% see FIG2_v3 line 498

% Coherence wise differences

% Hit rate agent vs solo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = generate_table(in, str, flag)
% Intialise
dat             = [];
id              = [];
coh             = [];

for iSubj = 1:length(in)
    if strcmp(flag, 'cr')
        tmp_dat     = in{iSubj}.(str);
        tmp_coh     = in{iSubj}.coh;
        tmp_dat     = check_dimensions(tmp_dat);
        tmp_coh     = check_dimensions(tmp_coh);
        dat         = [dat; tmp_dat];
        id          = [id; repmat(iSubj,length(tmp_dat),1)];
        coh         = [coh; tmp_coh];
    elseif strcmp(flag, 'perf')
        for iCoh = 1:length(in{iSubj}.carr)
            tmp_dat = in{iSubj}.(str){iCoh};
            tmp_dat = check_dimensions(tmp_dat);
            dat     = [dat; tmp_dat];
            id      = [id; repmat(iSubj,length(tmp_dat),1)];
            coh     = [coh; repmat(in{iSubj}.carr(iCoh),length(tmp_dat),1)];
        end
    elseif strcmp(flag, 'hir')
        dat         = [dat; in{iSubj}.hir'];
        coh         = [coh; in{iSubj}.carr'];
        id          = [id; repmat(iSubj,length(in{iSubj}.hir),1)];
    end
end

% Shuffle coherence labels -> control model
coh_shuffled    = coh(randperm(length(coh)));

% Build table
out             = table(id,coh,coh_shuffled,dat, 'VariableNames',{'ID','Coh','Coh_shuffled','Dat'});
end

function [lmeStat,comp] = LME_coherence_model(in)

% Random Intercept and Random Slope for One Predictor
lme1            = fitlme(in, 'Dat ~ Coh_shuffled + ( Coh_shuffled | ID)');
lme2            = fitlme(in, 'Dat ~ Coh + ( Coh | ID)'); % Coherence as predictor
[~,~,lmeStat]   = fixedEffects(lme2);

% Compare models
comp          	= compare(lme1,lme2);
end

function out = check_dimensions(in)
if size(in,1) == 1
    out = in';
else
    out = in;
end
end