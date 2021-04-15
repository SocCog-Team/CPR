addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/ownCloud/Documents/Matlab/CPR

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
    'fxs_cpr_20210204_fxs.mwk2_tbl.mat';
    'fxs_cpr_20210318_sut.mwk2_tbl.mat';
    'fxs_cpr_20210204_kan.mwk2_tbl.mat';
    'fxs_cpr_20210204_mac.mwk2_tbl.mat';
    'fxs_cpr_20210210_nes.mwk2_tbl.mat';
    'fxs_cpr_20210205_pas.mwk2_tbl.mat';
    'fxs_cpr_20210205_stm.mwk2_tbl.mat';
    'fxs_cpr_20210205_sem.mwk2_tbl.mat';
    'fxs_cpr_20210205_lac.mwk2_tbl.mat';
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

lg = legend('CPR','', '4AFC','', 'Location', 'southeast');



