addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab

%% Import files

fname_cpr = {   
    '20210204_bew_CPR_fxs.mwk2';
    '20210204_fxs_CPR_fxs.mwk2';
    '20210204_kan_CPR_fxs.mwk2';
    '20210204_mac_CPR_fxs.mwk2';
    '20210205_lac_CPR_fxs.mwk2';
    '20210205_pas_CPR_fxs.mwk2';
    '20210205_sem_CPR_fxs.mwk2';
    '20210205_stm_CPR_fxs.mwk2';
    '20210210_ilv_CPR_fxs.mwk2';
    '20210210_nes_CPR_fxs.mwk2';
    '20210224_piy_CPR_fxs.mwk2';
    '20210318_sut_CPR_fxs.mwk2';
    
    '20210407_piy_CPR_eye_sut.mwk2';
    '20210407_stm_CPR_eye_sut.mwk2';
    '20210408_mac_CPR_eye_sut.mwk2';
    '20210408_nes_CPR_eye_sut.mwk2';
    '20210409_kan_CPR_eye_sut.mwk2';
    '20210413_pas_CPR_eye_sut.mwk2';
    '20210415_lac_CPR_eye_sut.mwk2';    
    };

for iSubj = 1:size(fname_cpr,1)
    if iSubj < 13
        pth     = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/Pilot_free_viewing/';
    else
        pth     = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/Pilot_eye_fixation/';
    end
    
    id                  = fname_cpr{iSubj}(10:12);
    [t,d,out{iSubj}]    = CPR_psych_import(id, pth, fname_cpr{iSubj});
end

%% Group to population table

tbl_free            = [];
tbl_eye             = [];

for iSubj = 1:size(fname_cpr,1)
    if iSubj < 13
        pth         = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/Pilot_free_viewing/';
    else
        pth         = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/Pilot_eye_fixation/';
    end
    
    id              = fname_cpr{iSubj}(10:12);
    dte             = fname_cpr{iSubj}(1:8); 
    
    load([pth id '_' dte '_tbl.mat'])
    
    if iSubj < 13
        tbl_free    = [tbl_free; t];
    else
        tbl_eye     = [tbl_eye; t];
    end
end

%% Save data

% Free-viewing
t           = [];
t           = tbl_free;
pth         = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/Pilot_free_viewing/';
save([pth 'CPR_pop_tbl.mat'], 't', '-v7.3') 
tmp         = out;
out         = [];
out         = tmp(1:12);
save('pop_summary.mat', 'out', '-v7.3')

% Gaze-controlled
t           = [];
t           = tbl_eye;
pth      	= '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/Pilot_eye_fixation/';
save([pth 'CPR_pop_tbl.mat'], 't', '-v7.3')  % Save as .mat file
out         = [];
out         = tmp(13:end);
save('pop_summary.mat', 'out', '-v7.3')


