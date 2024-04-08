% Reads ND2 images.
%
% This class requires a C++ compiler to be installed, a 64-bit 
% version of Matlab and a Windows OS in order to be able to use 
% the ND2SDK library. 
%
% If either of these three requirements fail then the slower 
% bioformats library is used to read nd2 files.
%
% Read example
% nd2 = ND2('sample.nd2')
% [pixels, t, x, y, z] = nd2.read(1) % read the acquisition information for the first frame
% nd2.close() % must close the connection to the ND2 SDK library when you are done with this file
%
classdef ND2 < handle    
    properties
        numFrames; % int, number of frames
        width; % int, frame width
        height; % int, frame height
        filename; % string, the path of the nd2 file
        textInfo = struct(); % struct of the text information, LIMTEXTINFO
        metadata_Desc = struct(); % struct of the text information, METADATA_DESC        
        attributes = struct(); % struct of file attributes, LIMATTRIBUTES
        picture = struct(); % struct of the picture information, LIMPICTURE
        localMetadata = struct(); % struct of the local metadata, LIMLOCALMETADATA
        handle = 0; % the dll handle
    end
    properties(Hidden)
        dll = 'v6_w32_nd2ReadSDK'; % nd2sdk library name
        usingSDK; % boolean, indicates if using the nd2sdk dll
        verbose; % boolean, indicates whether to display message to the command window
        BFmeta = struct(); % bioformats metadata struct
        libname = '';
        hfile = '';
    end    
    methods        
        function self = ND2(filename, verbose)
        % function self = ND2(filename, verbose)
        % Class constructor
        %
        % Inputs
        % ------
        % filename : string
        %   then path of the nd2 file
        %
        % verbose (optional argument) : boolean
        %   display messages to the screen, default is true
            
            % default values
            self.numFrames = 0;
            self.width = 0;
            self.height = 0;
            self.filename = '';
            
            if nargin < 2
                self.verbose = true;
            else
                self.verbose = verbose;
            end
        
            % make sure that the file exists
            if exist(filename, 'file') ~= 2
                error('ERROR ND2() :: ND2 file does not exist\n\t%s\n',filename);
            end
            
            % load the library if it was not already loaded
            self.usingSDK = true;
            if ~libisloaded(self.dll)
                dir = [strrep(mfilename('fullpath'), 'ND2', ''), 'nd2lib\'];
                warning off MATLAB:loadlibrary:TypeNotFoundForStructure          
                try
                    self.libname = fullfile(dir, self.dll);
                    self.hfile = fullfile(dir, 'nd2ReadSDK');
                    loadlibrary(self.libname, self.hfile);
                catch
                    self.usingSDK = false;
                    msg = ['ND2 WARNING! Using the slower bioformats library to read nd2 files.', char(10), ...
                          'You need a 64-bit version of Matlab, a Windows OS and a C++ compiler', char(10), ...
                          'installed in order to be able to use the fast ND2 reader.', char(10), ...
                          'If your version of Matlab and OS meet these requirements then you need', char(10), ...
                          'to run "mex -setup" to tell Matlab which C++ compiler to use.'];
                    if self.verbose
                        fprintf('%s\n\n',msg);
                        str = input('Do you want to run the mex command [y]/n? ', 's');
                        if isempty(str)
                            str = 'y';
                        end
                    end
                end
                warning on MATLAB:loadlibrary:TypeNotFoundForStructure
            end

            self.filename = filename;

            % initialize this nd2 file
            if self.usingSDK
                self.initializeSDK();
            else
                if ~isempty(strfind(str, 'n')) || ~isempty(strfind(str, 'N'))
                    self.initializeBF();
                else
                    mex -setup
                end
            end
            
        end       
        
        function [pixels, timestamp, xPos, yPos, zPos] = read(self, frameNumber)
        % function [pixels, timestamp, xPos, yPos, zPos] = read(frameNumber)
        %
        % Returns the acquisition data for this frame. The timestamp, xPos, yPos
        % and zPos values only make sense if you are using the nd2sdk
        % library to get the frame data. Otherwise these values are always 
        % equal to 0 if using the bioformats library to read the frame.
        %
        % Inputs
        % ------
        % frame : integer
        %   the first frame is 1
        %  
        % Returns
        % -------
        % pixels : 2D array
        %   the image data
        %
        % timestamp : double
        %   the relative time, in seconds, from the first frame
        %
        % xPos : double
        %   the X position of the stage for the current frame
        %
        % yPos : double
        %   the Y position of the stage for the current frame
        %
        % zPos : double
        %   the Z position of the stage for the current frame
        %
            
            timestamp = 0.0;
            xPos = 0.0; 
            yPos = 0.0;
            zPos = 0.0;
            
            % check that the frame number makes sense
            if (frameNumber > 0 && frameNumber <= self.numFrames)
                
                if self.usingSDK
                    
                    if self.handle == 0
                        error('ERROR ND2.read() :: You have already closed the ND2 file. Cannot read frame %d\n', frameNumber);
                    end
                    
                    % get the frame information
                    uiSeqIndex = floor((frameNumber-1) / self.picture.uiComponents);
                    [err, self.picture, self.localMetadata] = calllib(self.dll, 'Lim_FileGetImageData', self.handle, uiSeqIndex, self.picture, self.localMetadata);
                    if err < 0
                        error('ERROR ND2.read() :: Not able to read frame %d\n', frameNumber);
                    end
                    
                    % get the timing and microscope-stage position
                    timestamp = self.localMetadata.dTimeMSec * 1e-3;
                    xPos = self.localMetadata.dXPos;
                    yPos = self.localMetadata.dYPos;
                    zPos = self.localMetadata.dZPos;
                    
                    % get the pointer to the pixel data
                    n = self.attributes.uiWidth * self.attributes.uiHeight * self.picture.uiComponents;
                    if self.picture.uiBitsPerComp == 8
                        setdatatype(self.picture.pImageData,'uint8Ptr',n,1)
                    elseif self.picture.uiBitsPerComp == 32
                        setdatatype(self.picture.pImageData,'singlePtr',n,1)
                    else
                        setdatatype(self.picture.pImageData,'uint16Ptr',n,1)
                    end    
                    data = get(self.picture.pImageData);
                    
                    offset = 1 + mod(frameNumber-1, self.picture.uiComponents);
                    pixels = transpose(reshape(data.Value(offset:self.picture.uiComponents:n), self.width, self.height));

                else
                    if self.verbose
                        fprintf('Reading ND2 frame %d... ', frameNumber);
                    end
                    channel = 1 + mod(frameNumber-1, self.BFmeta.channels);
                    zplane = 1 + mod(floor((frameNumber-1)/self.BFmeta.channels), self.BFmeta.zsize);
                    tframe = 1 + floor((frameNumber-1)/(self.BFmeta.channels * self.BFmeta.zsize));
                    pixels = imreadBF(self.filename, zplane, tframe, channel);
                    if self.verbose
                        fprintf('DONE\n');
                    end
                end
            else
                self.close();
                error('ERROR ND2.read() :: Invalid frame number %d, value can be from 1 to %d\n',frameNumber,self.numFrames)
            end
            
        end
        
        function showMetadata(self)
        % function showMetadata()
        %
        % Prints some of the metadata to the Command Window
        %   
            if ~self.usingSDK                
                fprintf('X: %d\n', self.BFmeta.width);
                fprintf('Y: %d\n', self.BFmeta.height);
                fprintf('CHANNEL: %d\n', self.BFmeta.channels);
                fprintf('Z: %d\n', self.BFmeta.zsize);
                fprintf('TIME: %d\n', self.BFmeta.nframes);
                n = length(self.BFmeta.parameterNames);
                for i=1:n
                    fprintf('%s: %s\n', self.BFmeta.parameterNames(i), self.BFmeta.parameterValues(i));
                end
            else
                [year,month,day,hour,minu,sec,dayweek] = julian2greg(self.metadata_Desc.dTimeStart);            
                fprintf('Acquisition time: %s, %s\n', dayweek, datestr([year,month,day,hour,minu,sec]));
                fprintf('dAngle: %f\n', self.metadata_Desc.dAngle)
                fprintf('dCalibration: %f\n', self.metadata_Desc.dCalibration)
                fprintf('dAspect: %f\n', self.metadata_Desc.dAspect)
                disp(['wszObjectiveName: ' self.metadata_Desc.wszObjectiveName])
                fprintf('dObjectiveMag: %f\n', self.metadata_Desc.dObjectiveMag)
                fprintf('dObjectiveNA: %f\n', self.metadata_Desc.dObjectiveNA)
                fprintf('dRefractIndex1: %f\n', self.metadata_Desc.dRefractIndex1)
                fprintf('dRefractIndex2: %f\n', self.metadata_Desc.dRefractIndex2)
                fprintf('dPinholeRadius: %f\n', self.metadata_Desc.dPinholeRadius)
                fprintf('dZoom: %f\n', self.metadata_Desc.dZoom)
                fprintf('dProjectiveMag: %f\n', self.metadata_Desc.dProjectiveMag)
                fprintf('uiImageType: %d\n', self.metadata_Desc.uiImageType)
                fprintf('uiPlaneCount: %d\n', self.metadata_Desc.uiPlaneCount)
                fprintf('uiComponentCount: %d\n', self.metadata_Desc.uiComponentCount)

                disp(['wszImageID: ' self.textInfo.wszImageID])
                disp(['wszType: ' self.textInfo.wszType])
                disp(['wszGroup: ' self.textInfo.wszGroup])
                disp(['wszSampleID: ' self.textInfo.wszSampleID])
                disp(['wszAuthor: ' self.textInfo.wszAuthor])
                disp('wszDescription:')
                disp(self.textInfo.wszDescription)
                disp(['wszCapturing: ' self.textInfo.wszCapturing])
                disp(['wszSampling: ' self.textInfo.wszSampling])
                disp(['wszLocation: ' self.textInfo.wszLocation])
                disp(['wszDate: ' self.textInfo.wszDate])
                disp(['wszConclusion: ' self.textInfo.wszConclusion])
                disp(['wszInfo1: ' self.textInfo.wszInfo1])
                disp(['wszInfo2: ' self.textInfo.wszInfo2])
                disp(['wszOptics: ' self.textInfo.wszOptics])
            end
        end
        
        function close(self)
        % function close()
        %
        % Close the connection to the ND2SDK library
            if self.usingSDK && self.handle > 0
                %calllib(self.dll, 'Lim_DestroyPicture', self.picture);
                calllib(self.dll, 'Lim_FileClose', self.handle);
                self.handle = 0;
            end
        end
        
        function filename = getFilename(self)
        % function filename = getFilename()
        %
        % Returns the name of the nd2 file
            filename = self.filename;
        end
        
        function width = getWidth(self)
        % function width = getWidth()
        %
        % Returns the frame width
            width = self.width;
        end
        
        function height = getHeight(self)
        % function height = getHeight()
        %
        % Returns the frame height
            height = self.height;
        end
        
        function nframes = getNumberOfFrames(self)
        % function nframes = getNumberOfFrames()
        %
        % Returns the number of frames
            nframes = self.numFrames;
        end 
        
        function setVerbose(self, verbose)
        % function setVerbose(verbose)
        %
        % Set whether to display messages to the command window
        %
        % Inputs
        % ------
        % verbose : boolean
        %  if 'true' then display warning/error messages to the screen
        %
            self.verbose = verbose;
        end
        
        function show(self)
        % function show()
        %
        % Display the movie in a custom-written Matlab canvas
            canvas(self)
        end
        
    end
    
    methods(Hidden)
        
        function initializeBF(self)
        % Load the nd2 file using bioformat's loci_tools.jar file
            if self.verbose
                fprintf('Reading the ND2 metadata... ')
            end
            self.BFmeta = imreadBFmeta(self.filename);
            self.numFrames = self.BFmeta.nframes * self.BFmeta.channels * self.BFmeta.zsize;
            self.width = self.BFmeta.width;
            self.height = self.BFmeta.height;
            fprintf('DONE\n')
        end
        
        function initializeSDK(self)
        % Load the nd2 file using the nd2sdk
            
            % initialize the structs
            self.attributes.uiWidth = uint32(0);             % Width of images
            self.attributes.uiWidthBytes = uint32(0);        % Line length 4-byte aligned
            self.attributes.uiHeight = uint32(0);            % Height of images
            self.attributes.uiComp = uint32(0);              % Number of components
            self.attributes.uiBpcInMemory = uint32(0);       % Bits per component 8, 16 or 32 (for float image)
            self.attributes.uiBpcSignificant = uint32(0);    % Bits per component used 8 .. 16 or 32 (for float image)
            self.attributes.uiSequenceCount = uint32(0);     % Number of images in the sequence
            self.attributes.uiTileWidth = uint32(0);         % If an image is tiled size of the tile/strip 
            self.attributes.uiTileHeight = uint32(0);        % otherwise both zero 
            self.attributes.uiCompression = uint32(0);       % 0 (lossless), 1 (lossy), 2 (None)
            self.attributes.uiQuality = uint32(0);           % 0 (worst) - 100 (best)            

            self.picture.uiWidth = uint32(0);                % Width of images
            self.picture.uiHeight = uint32(0);               % Height of images
            self.picture.uiBitsPerComp = uint32(0);          % BPC 8, 16 or 32
            self.picture.uiComponents = uint32(0);           % Number of components
            self.picture.uiWidthBytes = uint32(0);           % aligned to 4-byte
            self.picture.uiSize = uint32(0);
            self.picture.pImageData = libpointer;            % void* 
            
            self.localMetadata.dTimeMSec = double(0);             % Relative time msec from the first frame
            self.localMetadata.dXPos = double(0);                 % Stage XPos
            self.localMetadata.dYPos = double(0);                 % Stage YPos
            self.localMetadata.dZPos = double(0);                 % Stage ZPos
            
            % wchar_t[256] in C becomes uint16(zeros(1, 256)) in Matlab
            self.textInfo.wszImageID = uint16(zeros(1, 256));
            self.textInfo.wszType = uint16(zeros(1, 256));
            self.textInfo.wszGroup = uint16(zeros(1, 256));
            self.textInfo.wszSampleID = uint16(zeros(1, 256));
            self.textInfo.wszAuthor = uint16(zeros(1, 256));
            self.textInfo.wszDescription = uint16(zeros(1, 4096));
            self.textInfo.wszCapturing = uint16(zeros(1, 4096));
            self.textInfo.wszSampling = uint16(zeros(1, 256));
            self.textInfo.wszLocation = uint16(zeros(1, 256));
            self.textInfo.wszDate = uint16(zeros(1, 256));
            self.textInfo.wszConclusion = uint16(zeros(1, 256));
            self.textInfo.wszInfo1 = uint16(zeros(1, 256));
            self.textInfo.wszInfo2 = uint16(zeros(1, 256));
            self.textInfo.wszOptics = uint16(zeros(1, 256));
            
            self.metadata_Desc.dTimeStart = double(0); % Absolute Time in JDN
            self.metadata_Desc.dAngle = double(0); % Camera Angle
            self.metadata_Desc.dCalibration = double(0); % um/px (0.0 = uncalibrated)
            self.metadata_Desc.dAspect = double(0); % pixel aspect (always 1.0)
            self.metadata_Desc.wszObjectiveName = uint16(zeros(1, 256)); % The name of the objective
            self.metadata_Desc.dObjectiveMag = double(0); % Optional additional information
            self.metadata_Desc.dObjectiveNA = double(0); % dCalibration takes into accont all these
            self.metadata_Desc.dRefractIndex1 = double(0);
            self.metadata_Desc.dRefractIndex2 = double(0);
            self.metadata_Desc.dPinholeRadius = double(0);
            self.metadata_Desc.dZoom = double(0);
            self.metadata_Desc.dProjectiveMag = double(0);
            self.metadata_Desc.uiImageType = int32(0); % 0 (normal), 1 (spectral)
            self.metadata_Desc.uiPlaneCount = int32(0); % Number of logical planes (uiPlaneCount <= uiComponentCount)
            self.metadata_Desc.uiComponentCount = int32(0); % Number of physical components (same as uiComp in LIMFILEATTRIBUTES)
            self.metadata_Desc.pPlanes = libpointer; 
            
            % sometimes matlab wouldn't open the file... strange
            % opening it in a loop, seems to help a bit but does not solve
            % the problem
            cnt = 0;
            while cnt < 5000
                % open this file, must convert the filename to uint16
                wszFileName = uint16(self.filename);
                self.handle = calllib(self.dll, 'Lim_FileOpenForRead', wszFileName);
                if self.handle == 0
                    cnt = cnt + 1;
                    try
                        unloadlibrary self.dll
                        loadlibrary(self.libname, self.hfile);
                    catch
                    end
                    clear wszFileName
                else
                    break;
                end
            end 
            if self.handle == 0
                error('ERROR ND2.initializeSDK() :: Failed to open file');
            end
            
            % read the file attributes
            [err, self.attributes] = calllib(self.dll, 'Lim_FileGetAttributes', self.handle, self.attributes);
            if err ~= 0
                calllib(self.dll, 'Lim_FileClose', self.handle);
                error('ERROR ND2.initializeSDK() :: Failed to get file attributes');
            end
            self.width = self.attributes.uiWidth;
            self.height = self.attributes.uiHeight;
            
            [err, self.textInfo] = calllib(self.dll, 'Lim_FileGetTextinfo', self.handle, self.textInfo);
            if err ~= 0
                calllib(self.dll, 'Lim_FileClose', self.handle);
                error('ERROR ND2.initializeSDK() :: Not able to read Lim_FileGetTextinfo');
            end
            self.textInfo.wszImageID = strrep(char(self.textInfo.wszImageID), char(10), '');
            self.textInfo.wszType = strrep(char(self.textInfo.wszType), char(10), '');
            self.textInfo.wszGroup = strrep(char(self.textInfo.wszGroup), char(10), '');
            self.textInfo.wszSampleID = strrep(char(self.textInfo.wszSampleID), char(10), '');
            self.textInfo.wszAuthor = strrep(char(self.textInfo.wszAuthor), char(10), '');
            self.textInfo.wszDescription = strrep(char(self.textInfo.wszDescription), char(10), '');
            self.textInfo.wszCapturing = strrep(char(self.textInfo.wszCapturing), char(10), '');
            self.textInfo.wszSampling = strrep(char(self.textInfo.wszSampling), char(10), '');
            self.textInfo.wszLocation = strrep(char(self.textInfo.wszLocation), char(10), '');
            self.textInfo.wszDate = strrep(char(self.textInfo.wszDate), char(10), '');
            self.textInfo.wszConclusion = strrep(char(self.textInfo.wszConclusion), char(10), '');
            self.textInfo.wszInfo1 = strrep(char(self.textInfo.wszInfo1), char(10), '');
            self.textInfo.wszInfo2 = strrep(char(self.textInfo.wszInfo2), char(10), '');
            self.textInfo.wszOptics = strrep(char(self.textInfo.wszOptics), char(10), '');

            [err, self.metadata_Desc] = calllib(self.dll, 'Lim_FileGetMetadata', self.handle, self.metadata_Desc);
            if err ~= 0
                calllib(self.dll, 'Lim_FileClose', self.handle);
                error('ERROR ND2.initializeSDK() :: Not able to read Lim_FileGetMetadata');
            end
            self.metadata_Desc.wszObjectiveName = strrep(char(self.metadata_Desc.wszObjectiveName), char(10), '');
            
            % initialize the picture
            [~, self.picture] = calllib(self.dll, 'Lim_InitPicture', self.picture, self.attributes.uiWidth, self.attributes.uiHeight, self.attributes.uiBpcSignificant, self.attributes.uiComp);
            self.numFrames = self.attributes.uiSequenceCount * self.picture.uiComponents;
            
        end
        
    end
    
end