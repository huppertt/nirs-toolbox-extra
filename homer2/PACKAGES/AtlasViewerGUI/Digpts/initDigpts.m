function digpts = initDigpts(handles)

digpts(1) = struct( ...
   'name', 'digpts', ...
   'pathname', '', ...
   'handles',struct( ...
        'hSrcpos',[], ...
        'hDetpos',[], ...
        'hOptodes',[], ...
        'hPcpos',[], ...
        'hRefpts',[], ...
        'radiobuttonShowDigpts', [], ...
        'menuItemRegisterAtlasToDigpts', [], ...
        'axes',[] ...
   ), ...
   'refpts', initRefpts(), ...
   'srcpos',[], ...
   'srcmap',[], ...
   'detpos',[], ...
   'detmap',[], ...
   'pcpos',[], ...
   'T_2mc',eye(4), ...
   'center',[], ...    
   'orientation', '', ...
   'checkCompatability',@checkDigptsCompatability, ...
   'isempty',@isempty_loc, ...
   'prepObjForSave',[], ...
   'digpts',[] ...
);

if exist('handles','var')
    digpts.handles.radiobuttonShowDigpts = handles.radiobuttonShowDigpts;
    digpts.handles.menuItemRegisterAtlasToDigpts = handles.menuItemRegisterAtlasToDigpts;
    set(digpts.handles.radiobuttonShowDigpts,'enable','off');
    set(digpts.handles.radiobuttonShowDigpts,'value',0);
    digpts.handles.axes = handles.axesSurfDisplay;    
end




% ----------------------------------------------------------------
function digpts = checkDigptsCompatability(digpts)

if ~isstruct(digpts.refpts) & isfield(digpts, 'labels')
    refpts = digpts.refpts;
    labels = digpts.labels;
    
    digpts.refpts = initRefpts();
    digpts = rmfield(digpts, 'labels');
    
    digpts.refpts.pos = refpts;
    digpts.refpts.labels = labels;
end



% --------------------------------------------------------------
function b = isempty_loc(digpts)

b = true;
if ~isempty(digpts.refpts.pos) & ~isempty(digpts.refpts.labels)
    b = false;
end 
if ~isempty(digpts.srcpos)
    b = false;
end 
if ~isempty(digpts.detpos)
    b = false;
end
if ~isempty(digpts.pcpos)
    b = false;
end


