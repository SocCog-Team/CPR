% Connect to PSY4, then change into data folder
cd /Volumes/cnl/Documents/MWorks/Data

% Import file
d = CPR_import_mwk2('20240605_hekane_CPRneutral_block5.mwk2', {'INFO_'}, false);

% P1 bonus
tmp = d.value(d.event == 'INFO_bonus_ply1_cents');
tmp{end}

% P2 bonus
tmp = d.value(d.event == 'INFO_bonus_ply2_cents');
tmp{end}