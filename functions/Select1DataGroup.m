% Jason Otterstrom, MATLAB 2013b
%
% Simple GUI that will ask the user to select files pertaining to two
% "types" that can be specified by the user as string inputs. The filenames
% + paths of these data-groups are returned by the program
%   
% Call function like this:
%   OutputFileLists = Select1DataGroup(type1str,filetypes,PathName)
%
% Inputs (all are optional)
%   type: string of characters displayed above the  listbox 
%       this string is associated with the output ".data" field
%   filetypes: a string specifying the file extension allowed during file
%       selection. If no input is specified, then the default value is
%       '*.*' to permit any file type to be chosen.
%       Example values could be 'bin' or '*.bin'
%   PathName: string specifying the location where the user would like to
%       start looking for files
%
% Output
%   OutputFileLists: structure with fields 'data' associated with the 
%       related field 'dataname' that includes the optionally
%       user-specified file group name.
%       The data field is a cell having the format:   
%           {filename, file_location}


function OutputFileLists = Select1DataGroup(type1str,filetypes,PathName)
% function Compare2ClusterResults(type1str,type2str)

%  Create and then hide the GUI as it is being constructed.
guisize = [360,500,400,350];
gui = figure('Visible','off','Position',guisize);
set(gui,'Name','Select data files in separate folders',...
    'Units','normalized',...
    'NumberTitle','off',... 
    'MenuBar','none',... 
    'Resize','off')

% initialize & check input
if nargin >= 2 && ~isempty(type1str) 
    type1str = [type1str ' Datasets'];
else
    type1str = 'Datasets';
end
type1field = type1str;
type1field( type1field==' ' ) = '_';
if ~exist('filetypes','var') || isempty(filetypes) 
    filetypes = '*.*';
elseif iscell(filetypes)
    % could be a column vector of only file extensions  (no change needed)
    % or could be a row vector of only file extensions
    if isrow(filetypes)
        % need to flip to be a column vector
        filetypes = transpose(filetypes);
    % or a cell array with file extensions in the first column and
     % a corresponding description in the second
     % more than two columns is not permitted
    elseif ismatrix(filetypes)
        if size(filetypes,2) > 2
            warning('Filetype specifier may only be a [m 2] sized cell array, removing columns > 2')
            filetypes = filetypes(:,1:2);
        end
        allFile = false;
        for rw = 1:size(filetypes,1)
            if strcmp(filetypes{rw,1},'*.*')
                allFile = true;
            end
        end
        if ~allFile
            filetypes{rw+1,1} = '*.*';
            filetypes{rw+1,2} = 'All Files';
        end
    end
    % ensure that all of the file extensions in the first column are
    % properly formatted
    cl = 1;
    for rw = 1:size(filetypes,1)
        filetypes{rw,cl} = checkFileTypeString( filetypes{rw,cl} );
    end
elseif ischar(filetypes) %~strcmp(filetypes(1:2),'*.')
    filetypes = checkFileTypeString( filetypes );
end
if ~exist('PathName','var') || isempty(PathName) || ~isdir(PathName)
    disp('Starting path set to current working directory')
    PathName = pwd;
end

% initialize output
FileNamesDefault = {};
OutputFileLists = struct('data',[],'dataname',type1field);
% setappdata(gui,'OutputFileLists',OutputFileLists);

% default sizes of things
Listboxsize = [300 200];
menuFontSize = 10;
titleFontSize = 13;
Textsep = 5;
ListboxTitlesize = [230,20];
SaveLoadButtonSize = [100 40];
closeButtonSize = [250 50];
UDsep = [10 80 120];
UDbuttonsize = [50,25]; 
grey = get(gui,'Color');

% define internal GUI variables
setappdata(gui,'Datasets',OutputFileLists);
setappdata(gui,'PathName',PathName)
setappdata(gui,'filetype',filetypes)
setappdata(gui,'val',[])


% make the 'type 1' file listbox
listposn_T1 = [25 80+SaveLoadButtonSize(2) Listboxsize];
h.Listbox_T1 = uicontrol(gui,'Style','listbox',...
    'String',FileNamesDefault,...OutputFileLists.data1,...
    'HorizontalAlignment','left',...
    'Position',listposn_T1,...    'Value',[],... no initial selection        'Max',T1max,'Min',0,...);% Max-Min>1 permits multi-selection; use value to identify handle
    'Callback',{@FileListbox_Callback},...
    'FontSize',menuFontSize);

h.text_T1 = uicontrol(gui,'Style','text',...
    'String',type1str,...
    'Position',[listposn_T1(1),listposn_T1(2)+listposn_T1(4)+Textsep,ListboxTitlesize],...
    'BackgroundColor',grey,...
    'HorizontalAlignment','left',...
    'FontSize',titleFontSize,...
    'FontWeight','bold');
% construct buttons to add/remove type 1 listbox items
h.AddFile1 = uicontrol(gui,'Style','pushbutton',...
    'String','Add',...
    'HorizontalAlignment','Center',...
    'Position',[listposn_T1(1)+Listboxsize(1)+UDsep(1),listposn_T1(2)+UDsep(3),UDbuttonsize],...
    'Callback',{@AddFile_Callback,'Listbox_T1','data'});%,...    'Max',T1max);

h.RemFile1 = uicontrol(gui,'Style','pushbutton',...
    'String','Remove',...
    'HorizontalAlignment','Center',...
    'Position',[listposn_T1(1)+Listboxsize(1)+UDsep(1),listposn_T1(2)+UDsep(2),UDbuttonsize],...
    'Callback',{@RemFile_Callback,'Listbox_T1','data'});%,...    'Max',T1max);

% % make the 'type 2' file listbox
% listposn_T2 = [listposn_T1([1 2])+[400 0] Listboxsize];
% h.Listbox_T2 = uicontrol(gui,'Style','listbox',...
%     'String',FileNamesDefault,...OutputFileLists.data2,...
%     'HorizontalAlignment','left',...
%     'Position',listposn_T2,...    'Value',[],... no initial selection        'Max',T2max,'Min',0,...);% Max-Min>1 permits multi-selection; use value to identify handle
%     'Callback',{@FileListbox_Callback},...
%     'FontSize',menuFontSize);
% h.text_T2 = uicontrol(gui,'Style','text',...
%     'String',type2str,...
%     'Position',[listposn_T2(1),listposn_T2(2)+listposn_T2(4)+Textsep,ListboxTitlesize],...
%     'BackgroundColor',grey,...
%     'HorizontalAlignment','left',...
%     'FontSize',titleFontSize,...
%     'FontWeight','bold');
% % construct buttons to add/remove type 2 listbox items
% h.AddFile2 = uicontrol(gui,'Style','pushbutton',...
%     'String','Add',...
%     'HorizontalAlignment','Center',...
%     'Position',[listposn_T2(1)+Listboxsize(1)+UDsep(1),listposn_T2(2)+UDsep(3),UDbuttonsize],...
%     'Callback',{@AddFile_Callback,'Listbox_T2','data2'});%,...    'Max',T2max);
% h.RemFile2 = uicontrol(gui,'Style','pushbutton',...
%     'String','Remove',...
%     'HorizontalAlignment','Center',...
%     'Position',[listposn_T2(1)+Listboxsize(1)+UDsep(1),listposn_T2(2)+UDsep(2),UDbuttonsize],...
%     'Callback',{@RemFile_Callback,'Listbox_T2','data2'});%,...    'Max',T2max);

% make buttons to allow users to save currently selected files or load a
% previously saved file list
SLbuttonY = 65;
h.SaveButton = uicontrol(gui,'Style','pushbutton',...
    'Position',[0.5*(guisize(3)-closeButtonSize(1)-SaveLoadButtonSize(1)), SLbuttonY,...
        SaveLoadButtonSize],...
    'HorizontalAlignment','Center',...
    'String','Save Lists',...
    'Callback',{@SaveButton_Callback},...
    'FontSize',titleFontSize,...
    'FontWeight','bold');
h.LoadButton = uicontrol(gui,'Style','pushbutton',...
    'Position',[0.5*(guisize(3)+closeButtonSize(1)-SaveLoadButtonSize(1)), SLbuttonY,...
        SaveLoadButtonSize],...
    'HorizontalAlignment','Center',...
    'String','Load Lists',...
    'Callback',{@LoadButton_Callback},...
    'FontSize',titleFontSize,...
    'FontWeight','bold');

% make a button to close the GUI once files are selected
h.closeButton = uicontrol(gui,'Style','pushbutton',...
    'String','Done Selecting Files',...
    'HorizontalAlignment','Center',...
    'Position',[ 0.5*(guisize(3)-closeButtonSize(1)),5, closeButtonSize],...
    'Callback',{@closeButton_Callback},...
    'FontSize',titleFontSize,...
    'FontWeight','bold');

% save the handles in guidata
guidata(gui,h);

% Place the GUI on the screen.
movegui(gui,'center')%'north')%
% Make the GUI visible.
set(gui,'Visible','on');

% output the files structure: 
uiwait(gui)
if nargout > 0 && ishandle(gui)
    OutputFileLists = getappdata(gui,'Datasets');
    close(gui)
else
    OutputFileLists = [];
end
if ishandle(gui)
    close(gui)
end




% callback functions

    % get the value of the item selected in the list box
    function FileListbox_Callback(source,~) 
        val = get(source,'Value');
        setappdata(gui,'val',val);
    end

    % add files to the listbox via "Browse" functionality
    function AddFile_Callback(~,~,ListboxName,dataname)
        handles = guidata(gcbf);
        filetype = getappdata(gui,'filetype');
        currPath = getappdata(gui,'PathName');
        if ischar(filetype)
            asktext = ['Select desired ' filetype(2:end) ' files'];
            [Files, PathNameIN] = uigetfile(fullfile(currPath,filetype),asktext,...
                'MultiSelect','on');
        elseif iscell(filetype)
            item = str2double(ListboxName(end));
            asktext = ['Select desired ' filetype{item}(2:end) ' files'];
            disp('Starting path definition not available when specifying multiple file extensions')
            [Files, PathNameIN] = uigetfile(filetype,asktext,...
                'MultiSelect','on');
        else
            error('Unrecognized format for filetype input')
        end
            
        if length(Files) > 1
                       
            if iscell(Files)
                newFiles = [Files', repmat({PathNameIN},size(Files,2),1)];
            else
                newFiles = {Files, PathNameIN};
            end
            % update the stored datasets
            setappdata(gui,'PathName',PathNameIN)
            currFiles = getappdata(gui,'Datasets');
            currFiles.(dataname) = [currFiles.(dataname); newFiles];
            % check this list for duplicates
            currDispFiles = get(handles.(ListboxName),'String');
            [newFilesPaths,displayFiles] = CheckRepeats(currFiles.(dataname),currDispFiles);
            % update stored datasets again along with displayed files
            currFiles.(dataname) = newFilesPaths;
            setappdata(gui,'Datasets',currFiles)
            set(handles.(ListboxName),'String',displayFiles);
        end
    end

    % Remove the selected item from the listbox
    function RemFile_Callback(~,~,ListboxName,dataname)
        handles = guidata(gui);
        currDispFiles = get(handles.(ListboxName),'String');
        if ~isempty(currDispFiles)
            val = getappdata(gui,'val');
            if isempty(val) || val==0
                val = 1;
            end
            idx = (1:size(currDispFiles,1))';
            idxkeep = setxor(idx,val);
            if val > length(idxkeep)
                val = length(idxkeep);
                setappdata(gui,'val',length(idxkeep))
                set(handles.(ListboxName),'Value',val)
            end
            % update displayed & stored datasets
            set(handles.(ListboxName),'String',currDispFiles(idxkeep))
            currFiles = getappdata(gui,'Datasets');
            currFiles.(dataname) = currFiles.(dataname)(idxkeep,:);
            setappdata(gui,'Datasets',currFiles)
            % check for removal of duplicates
            currDispFiles = get(handles.(ListboxName),'String');
            [newFilesPaths,displayFiles] = CheckRepeats(currFiles.(dataname),currDispFiles);
            % update stored datasets again along with displayed files
            currFiles.(dataname) = newFilesPaths;
            setappdata(gui,'Datasets',currFiles)
            set(handles.(ListboxName),'String',displayFiles);
        else
            defaultFile = getappdata(gui,'FileNamesDefault');
            set(handles.(ListboxName),'String',defaultFile)
        end
        currval = get(handles.(ListboxName),'Value');
        if currval == 0
            set(handles.(ListboxName),'Value',1)
        end
    end

    % function to save the currently selected file list
    function SaveButton_Callback(~,~)
        currPathName = getappdata(gui,'PathName');
        if ~strcmp(currPathName(end),'\'), 
            currPathName = [currPathName '\'];
        end
        Currently_Selected_Files = getappdata(gui,'Datasets');
        Save_Date = datestr(now,'mmmm dd, yyyy HH:MM:SS');
        uisave({'Currently_Selected_Files','Save_Date'},...
            [currPathName 'Selected Files.mat']);
    end

    % function to select a .mat file containing a previously saved file list
    function LoadButton_Callback(~,~)
        handles = guidata(gcbf);
        currPathName = getappdata(gui,'PathName');
        [fl,pth] = uigetfile([currPathName,'\*.mat'],'Choose MAT file with file-list',...
            'MultiSelect','off');
        prevFiles = load( fullfile(pth,fl) );
        if ~isfield(prevFiles,'Currently_Selected_Files')
            warndlg({'The MATLAB file selected does not have the appropriate data.';...
                'Please try again.'})
        else
            currFiles = prevFiles.Currently_Selected_Files;
            ListboxName = {'Listbox_T1'};
            dataname = {'data'};
            for lb = 1
                currDispFiles = get(handles.(ListboxName{lb}),'String');
                [newFilesPaths,displayFiles] = CheckRepeats(currFiles.(dataname{lb}),currDispFiles);
                % update stored datasets again along with displayed files
                currFiles.(dataname{lb}) = newFilesPaths;
                setappdata(gui,'Datasets',currFiles)
                set(handles.(ListboxName{lb}),'String',displayFiles);
            end
        end
    end

    % Button to close the GUI and return the values of the selected files
    function closeButton_Callback(~,~)
        uiresume(gui)
    end

% subfunctions
    % check that the file types input variable is properly formatted
    function outputstr = checkFileTypeString( inputstr )
        if ~strcmp(inputstr(1:2),'*.')
            idxdot = strfind(inputstr,'.');
            if isempty(idxdot)
                outputstr = ['*.' inputstr];
            else
                outputstr = ['*' inputstr(idxdot:end)];
            end
        else
            outputstr = inputstr;
        end
    end

    % function to see if there are repeats in the files selected
    function [newFilesPaths,displayFiles] = CheckRepeats(currFilesPaths,currDispFiles)
        % compare filename & path, if two names same, add part of path to 
        % displayed value; if name & path equal, remove one
%         currFilesPaths = getappdata(gui,'Datasets');
        % 1) find replicated file names
        repfiles = nan(size(currFilesPaths,1)^2,2);
        m = 0;
        for f1 = 1:size(currFilesPaths,1)
            for f2 = f1+1:size(currFilesPaths,1)
                if strcmp(currFilesPaths(f1,1),currFilesPaths(f2,1))
                    m = m+1;
                    repfiles(m,:) = [f1, f2];
                end
            end
        end
        repfiles = repfiles(~isnan(repfiles(:,1)),:);
        % 2) if there are replicated files, check to see if the paths are different
        if isempty(repfiles)
            newFilesPaths = currFilesPaths;
            if isempty(currDispFiles)
                displayFiles = currFilesPaths(:,1);
            else
                displayFiles = currDispFiles;
            end
            % check to make sure all current FilesPaths are displayed
            for fp = 1:size(currFilesPaths,1)
                if isempty( cell2mat( strfind( displayFiles,currFilesPaths{fp,1} ) ) )
                    displayFiles = [displayFiles; currFilesPaths(fp,1)];
                end
            end
        else
            % initialize guess that they are all from different paths
            displayFiles = cell(size(currFilesPaths,1),1);
            newFilesPaths = cell(size(currFilesPaths));
            % find which files do not have repeats & make no changes
            nonrep = setxor((1:size(currFilesPaths,1))',unique(repfiles(:)));
            for f = 1:size(currFilesPaths,1)
                if ismember(f,nonrep)
                    newFilesPaths(f,:) = currFilesPaths(f,:);
                    if isempty(currDispFiles) || f>size(currDispFiles,1)
                        displayFiles(f,1) = currFilesPaths(f,1);
                    else
                        displayFiles(f,1) = currDispFiles(f);
                    end
                end
            end
            
            elim = nan(size(repfiles,1),1);
            m = 0;
            for r = 1:size(repfiles,1)
                if strcmp(currFilesPaths(repfiles(r,1),2),currFilesPaths(repfiles(r,2),2))
                    % then paths are identical, keep the first instance only.
                    m = m+1;
                    elim(m) = repfiles(r,2);
                end
            end
            elim = elim(~isnan(elim));
            for r = 1:size(repfiles,1)
                for rr = 1:2
                    idx = repfiles(r,rr);
                    if ~ismember(idx,elim)
                        % the paths are different, keep both file+path adapting
                        % displayed version
                        newFilesPaths(idx,:) = currFilesPaths(idx,:);
                        if rr==1 && ismember(repfiles(r,2),elim)
                            % no need to change displayed filename
                            displayFiles{idx,1} = currDispFiles{idx};
%                             displayFiles{idx,1} = currFilesPaths{idx,1};
                        else
                            % keep the last 5 characters of path, adding "..." before
                            displayFiles{idx,1} = ['...' currFilesPaths{idx,2}(end-5:end) ...
                                currFilesPaths{idx,1}];
                        end
                    end
                end
            end
            % remove empty cells
            displayFiles = displayFiles(~cellfun('isempty',displayFiles),:);
            newFilesPaths = newFilesPaths(~cellfun('isempty',newFilesPaths(:,1)),:);
        end
    end
    

end



