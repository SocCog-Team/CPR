addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Pilot
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab

% Import 4AFC data
var_import = {
    'RDP_', ...
    'TRIAL_'};

fname_4afc = {
    '4afc_fxs_20210317.mwk2';
    '4afc_sut_20210317.mwk2';
    '4afc_kan_20210319.mwk2';
    '4afc_mac_20210319.mwk2';
    '4afc_nes_20210319.mwk2';
    '4afc_pas_20210322.mwk2';
    '4afc_stm_20210323.mwk2';
    '4afc_sem_20210324.mwk2';
    '4afc_lac_20210324.mwk2';
    };

fname_cpr = {
    'fxs_20210204_tbl.mat';
    'sut_20210318_tbl.mat';
    'kan_20210204_tbl.mat';
    'mac_20210204_tbl.mat';
    'nes_20210210_tbl.mat';
    'pas_20210205_tbl.mat';
    'stm_20210205_tbl.mat';
    'sem_20210205_tbl.mat';
    'lac_20210205_tbl.mat';
    };

for iSubj = 1:size(fname_cpr,1)  
    % Import 4AFC data file
    cd /Users/fschneider/ownCloud/CPR_data/4AFC
    d = MW_readFile(fname_4afc{iSubj}, 'include', var_import);              % Import .mwk2 sesion file
    
    % Load CPR data table
    cd /Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/Pilot_free_viewing/
    load(fname_cpr{iSubj})
    
    % Fit psychometric function
    ax = subplot(3,3,iSubj);
    hold on
    CPR_fit_function(t,d)
    ax.Title.String = fname_4afc{iSubj}(6:8);
end

lg = legend('CPR bin','', '4AFC','','CPR acc','', 'Location', 'southeast');
lg.FontSize = 12;
