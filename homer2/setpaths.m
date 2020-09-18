function setpaths(add, fluence_simulate)

%
% USAGE: 
%
%   paths_str = setpaths(add, fluence_simulate)
%
% DESCRIPTION:
%
%   Sets all the paths needed by the homer2 and atlasViewer source code. 
%
% INPUTS:
%
%   add: if false or ommitted then setpaths adds all the homer2 paths,
%           specified in getpaths. If true, adds all the homer2 paths 
%           specified in getpaths. 
%
%   fluence_simulate: If true, then setpaths searches for precalculated fluence
%                       profile files used by AtlasViewer. When found it tries
%                       to add a second fluence profile to the fluence file
%                       to simulate a second wavelength if an actual
%                       fluence profile at a second wavelength doesn't
%                       already exist.
%
% OUTPUTS:
%
%   path homer2  
%    
% EXAMPLES:
%
%   setpaths;
%   paths = setpaths(0);
%   paths = setpaths(1);
%   paths = setpaths(0,0);
%   paths = setpaths(0,1);
%   paths = setpaths([],1);
%



if ~exist('add','var') | isempty(add)
    add = true;
end
if ~exist('fluence_simulate','var') | isempty(fluence_simulate)
    fluence_simulate = false;
end

paths = getpaths();

rootpath = pwd;
k = find(rootpath=='\');
rootpath(k)='/';

paths_str = '';
err = false;
for ii=1:length(paths)
    paths{ii} = [rootpath, paths{ii}];
    if ~exist(paths{ii}, 'dir')
        err = true;
        continue;
    end
    if isempty(paths_str)
        paths_str = paths{ii};
    else
        paths_str = [paths_str, ';', paths{ii}];
    end
end

if err
    menu('WARNING: The current folder does NOT look like a homer2 root folder. Please change current folder to the root homer2 folder and rerun setpaths.', 'OK');
    paths = {};
    return;
end

if add
    fprintf('ADDED homer2 paths to matlab search paths:\n');
    addpath(paths_str, '-end')
    
    fprintf('\n');

    if isunix()
        idx = findExePaths(paths);
        for ii=1:length(idx)
            fprintf(sprintf('chmod 755 %s/*\n', paths{idx(ii)}));
            files = dir([paths{idx(ii)}, '/*']);
            if ~isempty(files)
                system(sprintf('chmod 755 %s/*', paths{idx(ii)}));
            end
        end
    end
    if fluence_simulate
        genMultWlFluenceFiles_CurrWorkspace;
    end
    
    [r(1), toolboxes1] = checkToolboxes_Homer2();    
    [r(2), toolboxes2] = checkToolboxes_AtlasViewer();
    
    fprintf('\n');
    if all(r==1)
        fprintf('All required toolboxes are installed.\n');
    elseif ismember(3, r) 
        fprintf('Unable to verify if all required toolboxes are installed ...\n');
    elseif ismember(4, r) 
        fprintf('Unable to verify if all required toolboxes are installed with older Matlab release...\n');
    else
        fprintf('Some required toolboxes are missing...\n');
    end
    
    % Check if wavelet data db2.mat is available in toolbox.
    % If no then create it from known data
    fullpathhomer2 = fileparts(which('Homer2_UI.m'));
    if fullpathhomer2(end)~='/' & fullpathhomer2(end)~='\'
        fullpathhomer2(end+1)='/';
    end
    findWaveletDb2([fullpathhomer2, 'UTILITIES/Wavelet/']);
    
    pause(2);
    % open([fullpathhomer2, 'PACKAGES/Test/Testing_procedure.pdf']);
    msg{1} = sprintf('For instructions to perform basic tests of Homer2_UI and AtlasViewerGUI, open the PDF file %s', ...
                      [fullpathhomer2, 'PACKAGES/Test/Testing_procedure.pdf']);
    fprintf('\n\n*** %s ***\n\n', [msg{:}]);
    
else
    fprintf('REMOVED homer2 paths from matlab search paths:\n');
    rmpath(paths_str);
end



% ---------------------------------------------------
function idx = findExePaths(paths)

idx = [];
kk = 1;
for ii=1:length(paths)
    if ~isempty(findstr(paths{ii}, '/bin'))
        idx(kk) = ii;
        kk=kk+1;
    end
end

