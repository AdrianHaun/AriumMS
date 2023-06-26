function [PeaksMS1,timesMS1,varargout] = readmzXML(DataPath,options)
%% reads mzXML files and outputs mz, Intensity and time data
% .mzXML files can be 32bit and 64bit encoded can use zlib compression
% Input: Datapath: datapath to .mzML file
% optional inputs:
% StartTime: excludes scans with lower scan time, default is 0
% EndTime: excludes scans with higher scan time, default is Inf
% MSLevel: specifies the maximum MS level to be included and enables additional
% outputs, default is 1
% Outputs: 
% Peaks: cell array each containing a two column matrix, with mz values and corresponding intensity
% times: vector with scan times corresponding to each cell of Peaks.
% optional outputs: two column matrix, column1:MSn Precursor mass; column2: corresponding level
% Peak Data for MS2 then time Data for MS2, repeatable until all MSn levels
% specifided by the MSLevel option
arguments
    DataPath (1,1) string
    options.StartTime (1,1) {mustBeNumeric} = 0
    options.EndTime (1,1) {mustBeNumeric} = Inf
    options.MSLevel (1,1) {mustBeInteger} = 1
end
S = readstruct(DataPath,'FileType','xml');
% get file precision, file compression & PC bitorder
precision = S.msRun.scan(1).peaks.precisionAttribute;
switch precision
    case 32
        precision = 'single';
    case 64
        precision = 'double';
end
compression = S.msRun.scan(1).peaks.compressionTypeAttribute;
switch compression
    case "none"
        compression = false;
    case "zlib"
        compression = true;
end
[~,~,endian] = computer;
%% data extraction
%remove empty scans
Scans=[S.msRun.scan.totIonCurrentAttribute]';
idx=Scans <= 0;
S.msRun.scan(idx)=[];
%get time and transform to numeric
times = [S.msRun.scan.retentionTimeAttribute]';
times = cellfun(@(x) sscanf(x,'PT %f'), times);
%remove scans outside time range
idx=times < options.StartTime | times > options.EndTime;
times(idx)=[];
S.msRun.scan(idx)=[];
%Get Peak Data and msLevel
PeaksString = [S.msRun.scan.peaks];
PeaksString = [PeaksString.Text];
msLevels = [S.msRun.scan.msLevelAttribute]';
%Get PrecursorMass
MS2scanStruct = [S.msRun.scan(msLevels>1)];
%preallocation
Peaks=cell(size(PeaksString,2),1);
PrecursorMass = zeros(size(MS2scanStruct))';
MSL=options.MSLevel;
for nScan = 1:size(PrecursorMass,1)
    PrecursorMass(nScan) = getfield(MS2scanStruct(nScan),'precursorMz','Text');
end
parfor s=1:size(PeaksString,2)
    if msLevels(s) > MSL
        Peaks{s}=0;
    else
    if compression == true
        intlist = decodeCompressed(PeaksString{s},precision,endian);
        Peaks{s} = reshape(intlist,2,[])';
    else
        intlist = decodeUncompressed(PeaksString{s},precision,endian);
        Peaks{s} = reshape(intlist,2,[])';
    end
    end
end
%split MsLevels
timesMS1 = times(msLevels == 1);
PeaksMS1 = Peaks(msLevels == 1);
if options.MSLevel > 1 && options.MSLevel <= max(msLevels)
    msLevels(msLevels==1)=[];
    varargout{1}=[PrecursorMass,msLevels];
    out=cell(options.MSLevel-1,2);
    for n = 2:options.MSLevel
        out{n-1,1} = Peaks(msLevels == n);
        out{n-1,2} = times(msLevels == n);
    end
    out=reshape(out,1,[]);
    varargout = [varargout,out];
elseif options.MSLevel > 1 && max(msLevels) < options.MSLevel
    out=cell(options.MSLevel-1,2);
    varargout{1}=0;
    for n = 2:options.MSLevel
        out{n-1,1} = cell(1,1);
        out{n-1,2} = 0;
    end
    out=reshape(out,1,[]);
    varargout = [varargout,out];
end
end
function out = decodeUncompressed(DataString,precision,endian)
out = typecast(matlab.net.base64decode(DataString),'int8');
out = typecast(out,precision);
if strcmp(endian,'L')
    out = swapbytes(out);
end
end
function out = decodeCompressed(DataString,precision,endian)
out = typecast(matlab.net.base64decode(DataString),'int8');
out = zmat(out,0,'zlib');
out = typecast(out,precision);
if strcmp(endian,'L')
    out = swapbytes(out);
end
end