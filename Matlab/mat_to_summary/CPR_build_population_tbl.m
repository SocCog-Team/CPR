dpth                        = '/Users/fschneider/Documents/CPR_psychophysics/'; cd(dpth)
folder                      = dir;
cpr                         = [];

for iFolder = 1:size(folder,1)
    if numel(folder(iFolder).name) == 3 || (contains(folder(iFolder).name, 'Dyad') && numel(folder(iFolder).name) == 6)
        
        if strcmp(folder(iFolder).name, 'RoH')
            continue
        end
        
        cd([dpth folder(iFolder).name '/summary/'])
        fname               = dir('*.mat');
        
        for iFile = 1:size(fname,1)
            f_split         = strsplit(fname(iFile).name,'_');
            
            if length(f_split{1}) == 8 && (strcmp(f_split{2},lower(folder(iFolder).name)) || contains(folder(iFolder).name, 'Dyad'))
                in          = load(fname(iFile).name);
                cpr         = [cpr; in.t];
            end
        end
    end
end

save('/Users/fschneider/Desktop/cpr_data.mat', 'cpr', '-v7.3')