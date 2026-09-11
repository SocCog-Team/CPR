function fn_addpath_plexon()
% FN_ADDPATH_PLEXON  Put the Plexon Offline SDK on the MATLAB path if needed.
%
% Silent if the SDK is already reachable; otherwise adds the first known
% location that exists (lab machine or local laptop). Guarded so it never
% raises the "nonexistent directory" warning from addpath.
%
% Add your machine's SDK folder to the candidate list if it differs.
%
% Felix Schneider, CNL  (drafted with Claude)

if exist('PL2GetFileIndex', 'file') == 2
    return   % already on the path
end

candidates = { ...
    '/Users/cnl/Desktop/CPR/PlexonMatlabOfflineFilesSDK', ...
    '/Users/fschneider/Documents/MATLAB/PlexonMatlabOfflineFilesSDK'};

for i = 1:numel(candidates)
    if isfolder(candidates{i})
        addpath(candidates{i});
        return
    end
end

warning('fn_addpath_plexon:notFound', ...
    'Plexon Offline SDK not found in known locations — add yours to fn_addpath_plexon.');

end % fn_addpath_plexon
