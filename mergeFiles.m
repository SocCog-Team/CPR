function d = mergeFiles(fname, var_import)

% If existing file...
if isfile([fname{1}(1:12) '_merged.mat'])
    tmp         = load([fname{1}(1:34) '_merged.mat']);
    d           = tmp.d;
else
    % Otherwise, import & merge...
    d.time      = [];
    d.event     = [];
    d.value     = [];
    
    for i = 1:size(fname,2)
        tmp    	= MW_readFile(fname{i}, 'include', var_import);
        
        d.time  = [d.time tmp.time];
        d.event = [d.event tmp.event];
        d.value = [d.value tmp.value];
    end
    
    disp('Save struct...')
    save([fname{1}(1:12) '_CPR_merged.mat'], 'd', '-v7.3')
    disp('Done!')
end
end
