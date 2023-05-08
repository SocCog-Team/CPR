spth = '/Volumes/T7_Shield/CPR_psychophysics/';
% spth = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/';
% spth = '/Users/fschneider/Documents/CPR_psychophysics/';
dpth = '/Users/fschneider/ownCloud/CPR_files/tbl/';

cd(spth);
folders = dir;

for i = 1:length(folders)
    disp(['Copy: ' folders(i).name])
    if length(folders(i).name) == 3 || length(folders(i).name) == 6
        new_dir = [dpth folders(i).name  '/summary/'];
        mkdir(new_dir)
        copyfile([spth folders(i).name '/summary/'], new_dir)
    end
    disp('Done!')
end