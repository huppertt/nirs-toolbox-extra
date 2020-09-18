function setup()

global h
global nSteps
global iStep


h = waitbar(0,'Installation Progress ...');
nSteps = 100;
iStep = 1;

if ismac()
    dirnameSrc = '~/Downloads/homer2_install/';
else
	dirnameSrc = [pwd, '/'];
end
dirnameDst = getAppDir();

% Uninstall
try
    if exist(dirnameDst,'dir')
        rmdir(dirnameDst, 's');
    end
catch ME
    msg{1} = sprintf('Error: Could not remove old installation folder. It might be in use by other applications.\n');
    msg{2} = sprintf('Try closing and reopening file browsers or any other applications that might be using the\n');
    msg{3} = sprintf('installation folder and then retry installation.');
    menu([msg{:}], 'OK');
    close(h);
    rethrow(ME)
end

platform = setplatformparams(dirnameSrc);

fprintf('Platform params:\n');
fprintf('  arch: %s\n', platform.arch);
fprintf('  mc_exe: %s\n', platform.mc_exe);
fprintf('  homer2_exe: %s\n', platform.homer2_exe{1});
fprintf('  atlasviewer_exe: %s\n', platform.atlasviewer_exe{1});
fprintf('  setup_exe: %s\n', platform.setup_exe{1});
fprintf('  setup_script: %s\n', platform.setup_script);
fprintf('  dirnameApp: %s\n', platform.dirnameApp);
fprintf('  mcrpath: %s\n', platform.mcrpath);
fprintf('  iso2meshmex: %s\n', platform.iso2meshmex{1});
fprintf('  iso2meshbin: %s\n\n', platform.iso2meshbin);

try
    if ispc()
        cmd = sprintf('IF EXIST %%userprofile%%\\desktop\\%s.lnk (del /Q /F %%userprofile%%\\desktop\\%s.lnk)', ...
            platform.atlasviewer_exe{1}, platform.atlasviewer_exe{1});
        system(cmd);
        
        cmd = sprintf('IF EXIST %%userprofile%%\\desktop\\%s.lnk (del /Q /F %%userprofile%%\\desktop\\%s.lnk)', ...
            platform.homer2_exe{1}, platform.homer2_exe{1});
        system(cmd);
        
        cmd = sprintf('IF EXIST %%userprofile%%\\desktop\\Test.lnk (del /Q /F %%userprofile%%\\desktop\\Test.lnk)');
        system(cmd);
    elseif islinux()
        if exist('~/Desktop/Homer2_UI.sh','file')
            delete('~/Desktop/Homer2_UI.sh');
        end
        if exist('~/Desktop/AtlasViewerGUI.sh','file')
            delete('~/Desktop/AtlasViewerGUI.sh');
        end
        % For symbolic links exist doesn't work if the file/folder that is
        % pointed to does not exist. Delete symbolic link unconditionally. 
        % If the link itself isn't there then you get only a warning from
        % matlab when you try to delete it.
        delete('~/Desktop/Test');
        if ~exist(platform.mcrpath,'dir') | ~exist([platform.mcrpath, '/mcr'],'dir') | ~exist([platform.mcrpath, '/runtime'],'dir')
            menu('Error: Invalid MCR path under ~/libs/mcr. Terminating installation...\n','OK');
        end
    elseif ismac()
        if exist('~/Desktop/Homer2_UI.command','file')
            delete('~/Desktop/Homer2_UI.command');
        end
        if exist('~/Desktop/AtlasViewerGUI.command','file')
            delete('~/Desktop/AtlasViewerGUI.command');
        end
        % For symbolic links exist doesn't work if the file/folder that is
        % pointed to does not exist. Delete symbolic link unconditionally. 
        % If the link itself isn't there then you get only a warning from
        % matlab when you try to delete it.
        delete('~/Desktop/Test');

        if ~exist(platform.mcrpath,'dir') | ~exist([platform.mcrpath, '/mcr'],'dir') | ~exist([platform.mcrpath, '/runtime'],'dir')
            menu('Error: Invalid MCR path under ~/libs/mcr. Terminating installation...\n','OK');
        end
    end
catch
    menu('Warning: Could not delete Desktop icons Homer2_UI and AtlasViewerGUI. They might be in use by other applications.', 'OK');
end

pause(2);

% Create destination folders
try 
    mkdir(dirnameDst);
catch ME
    msg{1} = sprintf('Error: Could not create installation folder. It might be in use by other applications.\n');
    msg{2} = sprintf('Try closing and reopening file browsers or any other applications that might be using the\n');
    msg{3} = sprintf('installation folder and then retry installation.');
    menu([msg{:}], 'OK');
    close(h);
    rethrow(ME)
end

try 
    mkdir([dirnameDst, 'Colin']);
    mkdir([dirnameDst, 'Colin/anatomical']);
    mkdir([dirnameDst, 'Colin/fw']);
    mkdir([dirnameDst, 'tMCimg']);
    mkdir([dirnameDst, 'tMCimg/bin']);
    mkdir([dirnameDst, 'tMCimg/bin/Win']);
    mkdir([dirnameDst, 'tMCimg/bin/Linux']);
    mkdir([dirnameDst, 'tMCimg/bin/Darwin']);
    mkdir([dirnameDst, 'Test']);
catch ME
    msg{1} = sprintf('Error: Could not create installtion subfolder. Installtion folder might be in use by other applications.\n');
    msg{2} = sprintf('Try closing and reopening file browsers or any other applications that might be using the\n');
    msg{3} = sprintf('installation folder and then retry installation.');
    menu([msg{:}], 'OK');
    close(h);
    rethrow(ME)
end

% Copy all the AtlasViewerGUI app folder files

% Important: For next 2 copyfile calls make sure to keep the destination
% the way it is, with the destination file name specified. This is important
% for mac installation because the executable is actually a directory.
% Copyfile only copies the contents of a folder so to copy the whole thing
% you need to specify the root foder same as the source.
for ii=1:length(platform.atlasviewer_exe)
    copyFileToInstallation([dirnameSrc, platform.atlasviewer_exe{ii}],  [dirnameDst, platform.atlasviewer_exe{ii}]);
end
for ii=1:length(platform.homer2_exe)
    copyFileToInstallation([dirnameSrc, platform.homer2_exe{ii}], [dirnameDst, platform.homer2_exe{ii}]);
    file = [dirnameSrc, platform.atlasviewer_exe{ii}];
end

% Copy all the Colin atlas folder files
copyFileToInstallation([dirnameSrc, 'headsurf.mesh'],         [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, 'headsurf2vol.txt'],      [dirnameDst, 'Colin/anatomical']);
if (copyFileToInstallation([dirnameSrc, 'headvol.vox.gz'],    [dirnameDst, 'Colin/anatomical']) < 0)
    copyFileToInstallation([dirnameSrc, 'headvol.vox'],       [dirnameDst, 'Colin/anatomical']);
end
copyFileToInstallation([dirnameSrc, 'headvol2ras.txt'],       [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, 'headvol_dims.txt'],      [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, 'headvol_tiss_type.txt'], [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, 'labelssurf.mat'],        [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, 'labelssurf2vol.txt'],    [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, 'pialsurf.mesh'],         [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, 'pialsurf2vol.txt'],      [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, 'refpts.txt'],            [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, 'refpts2vol.txt'],        [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, 'refpts_labels.txt'],     [dirnameDst, 'Colin/anatomical']);
copyFileToInstallation([dirnameSrc, platform.mc_exe],         [dirnameDst, 'tMCimg/bin/', platform.arch]);
copyFileToInstallation([dirnameSrc, 'db2.mat'],               dirnameDst);

% Check if there a fluence profile to load in this particular search path
fluenceProfFnames = dir([dirnameSrc, 'fluenceProf*.mat']);
for ii=1:length(fluenceProfFnames)
    copyFileToInstallation([dirnameSrc, fluenceProfFnames(ii).name],  [dirnameDst, 'Colin/fw']);
    genMultWavelengthSimInFluenceFiles([dirnameSrc, fluenceProfFnames(ii).name], 2, [dirnameDst, 'Colin/fw']);
end

copyFileToInstallation([dirnameSrc, 'projVoltoMesh_brain.mat'], [dirnameDst, 'Colin/fw']);
copyFileToInstallation([dirnameSrc, 'projVoltoMesh_scalp.mat'], [dirnameDst, 'Colin/fw']);

for ii=1:length(platform.iso2meshmex)
    % Use dir instead of exist for mex files because of an annoying matlab bug, where a
    % non existent file will be reported as exisiting as a mex file (exist() will return 3)
    % because there are other files with the same name and a .mex extention that do exist.
    % dir doesn't have this problem.
    if ~isempty(dir([dirnameSrc, platform.iso2meshmex{ii}]))
        fprintf('Copying %s to %s\n', [dirnameSrc, platform.iso2meshmex{ii}], dirnameDst);
        copyFileToInstallation([dirnameSrc, platform.iso2meshmex{ii}], dirnameDst);
        if isunix()
            system(sprintf('chmod 755 %s', [dirnameDst, '', platform.iso2meshmex{ii}]'));
        end
    else
        fprintf('ERROR: %s does NOT exist...\n', [dirnameSrc, platform.iso2meshmex{ii}]);
    end
end

copyFileToInstallation([dirnameSrc, 'Test'], [dirnameDst, 'Test'], 'dir');

% Create desktop shortcuts to Homer2_UI and AtlasViewerGUI
try
    if ispc()
        k = dirnameDst=='/';
        dirnameDst(k)='\';
        
        cmd = sprintf('call %s\\createShortcut.bat %s AtlasViewerGUI.exe', dirnameSrc(1:end-1), dirnameDst);
        system(cmd);
        
        cmd = sprintf('call %s\\createShortcut.bat %s Homer2_UI.exe', dirnameSrc(1:end-1), dirnameDst);
        system(cmd);
        
        cmd = sprintf('call %s\\createShortcut.bat %s Test', dirnameSrc(1:end-1), dirnameDst(1:end-1));
        system(cmd);
        
    elseif islinux()
        cmd = sprintf('sh %s/createShortcut.sh sh', dirnameSrc(1:end-1));
        system(cmd);
    elseif ismac()
        cmd = sprintf('sh %s/createShortcut.sh command', dirnameSrc(1:end-1));
        system(cmd);
    end
catch
    msg{1} = sprintf('Error: Could not create Homer2 shortcuts on Desktop. Exiting installtion.');
    menu([msg{:}], 'OK');
    return;    
end

waitbar(iStep/nSteps, h); iStep = iStep+1;
pause(2);


% Check that everything was installed properly
r = finishInstallGUI();

waitbar(iStep/nSteps, h); iStep = iStep+1;
pause(2);

if r==0
    try
        open([dirnameDst, 'Test/Testing_procedure.pdf']);
    catch ME
        msg{1} = sprintf('Warning at line 225 in setup.m: %s', ME.message);
        menu([msg{:}], 'OK');
        close(h);
        fprintf('Error at line 225 in setup.m: %s\n', ME.message); 
        rethrow(ME);
    end
end

waitbar(nSteps/nSteps, h);
close(h);

% cleanup();


% -----------------------------------------------------------------
function cleanup()

% Cleanup
if exist('~/Desktop/homer2_install/','dir')
    rmdir('~/Desktop/homer2_install/', 's');
end
if exist('~/Desktop/homer2_install.zip','file')
    delete('~/Desktop/homer2_install.zip');
end
if exist('~/Downloads/homer2_install/','dir')
    rmdir('~/Downloads/homer2_install/', 's');
end
if exist('~/Downloads/homer2_install.zip','file')
    delete('~/Downloads/homer2_install.zip');
end



% -------------------------------------------------------------------
function r = copyFileToInstallation(src, dst, type, errtype)

global h
global nSteps
global iStep

if ~exist('type', 'var')
    type = 'file';
end
if ~exist('errtype', 'var')
    errtype = 'Error';
end

r = 0;
try
    assert(logical(exist(src, type)));
    copyfile(src, dst);
    waitbar(iStep/nSteps, h); iStep = iStep+1;
    pause(1);
catch ME
    msg{1} = sprintf('Error at line 274 in setup.m: Could not copy %s ', src);
    msg{2} = sprintf('to installation folder. %s', ME.message);
    menu([msg{:}], 'OK');
    if strcmpi(errtype, 'Warning')
        r = -1;
    else
        close(h);
        fprintf('Error at line 274 in setup.m: %s\n', ME.message);
        rethrow(ME);
    end
end




