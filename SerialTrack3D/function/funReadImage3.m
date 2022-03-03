function [file_name,Img] = funReadImage3(varargin)
% [file_name,f,g,winsize,winstepsize,gridRange] = ReadImage3(filename) is
% the main function that loads 3D volumetric images
%
% INPUTS
% -------------------------------------------------------------------------
%   filename: string for the filename prefix for the volumetric images in 
%             the current directory.  
%             Input options:
%             --- If image is not within a cell) ---
%             1) 'filename*.mat' or 'filename*' 
% 
%             --- If image is within a cell that contains multichannels ---
%             2) filename{1} = 'filename*.mat' or 'filename*' and
%                filename{2} = channel number containing images you want to
%                              run IDVC on.
%                (if the channel is not provided, i.e. length(filename) = 1
%                , then channel = 1
%                 
% OUTPUTS
% -------------------------------------------------------------------------
%   file_name: Store 3D volumetric images file names
%   Img: reference and deformed image;
%       Img{1} = f: reference image
%       Img{2} = g: deformed image
%    
%
% NOTES
% -------------------------------------------------------------------------
% none
% 
% For more information please see 
%

% fprintf('Choose method to load images:  \n')
% fprintf('     0: Select images folder;  \n')
% fprintf('     1: Use prefix of image names;  \n')
% fprintf('     2: Manually select images.  \n')
% prompt = 'Input here: ';
% LoadImgMethod = input(prompt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---- Opening & Reading the First Image into CPU Memory ----
[fileInfo,ImgSeqNum] = parseInputs(varargin{:});
file_name = fileInfo.filename;
Img{1} = loadFile(fileInfo,1);

%% ---- Opening and Reading Subsequent Images ---
if ImgSeqNum>1
% numImages = length(fileInfo.filename);
for i = ImgSeqNum % 2:numImages % Reads Volumes Starting on the Second Volumes
    Img{2} = loadFile(fileInfo,i);
end
end

end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supp functions
function I = loadFile(fileInfo,idx)
I = load(fileInfo.filename{idx});
fieldName = fieldnames(I);
I = getfield(I,fieldName{1});
if iscell(I)
    if numel(I), I = I{1};
    else
        I = I{fileInfo.dataChannel};
    end
end
end

function varargout = parseInputs(varargin)
%  = parseInputs(filename, sSize, sSizeMin, runMode)
 
ImgSeqNum = varargin{2};

% Parse filenames
filename = varargin{1};
if iscell(filename)
    if length(filename) == 1, fileInfo.datachannel = 1;
    else fileInfo.datachannel = filename{2};
    end
    filename = filename{1};
end


[~,filename,~] = fileparts(filename);
filename = dir([filename,'.mat']);
fileInfo.filename = {filename.name};

if isempty(fileInfo), error('File name doesn''t exist'); end

% % Ensure dimensionality of the subset size
% sSize = varargin{2};
% if numel(sSize) == 1
%     sSize = sSize*[1 1 1];
% elseif numel(sSize) ~=3
%     error('Subset size must be a scalar or a three column array');
% end
% 
% %get minium subset size
% sSizeMin = varargin{3};
%                   
% % Ensure range of subset size
% if min(sSize) < 16 || max(sSize > 128)
%    error('Subset size must be within 16 and 128 pixels');
% end
% 
% % Ensure even subset size
% % if sum(mod(sSize,4)) > 0
% %     error('Subset size must be even');
% % end
% 
% if sum(mod(sSize,32)) ~= 0
%     error('Subset size must be 16, 32, 64, 96, or 128 voxels in each dimension');
% end

% Check run method input
% runMode = 'c';
% runMode  = varargin{4};

% if ~(strcmpi(runMode(1),'c') || strcmpi(runMode(1),'i') || strcmpi(runMode(1),'h'))
%     error('Run method must be incremental or cumulative or hybrid');
% end
% 
% % Initial guess of displacement field = [0 0 0];
% u0 = num2cell(zeros(1,3));

% Outputs
varargout{      1} = fileInfo;
% varargout{end + 1} = sSize;
% varargout{end + 1} = sSizeMin;
varargout{end + 1} = ImgSeqNum;
% varargout{end + 1} = u0;

end


