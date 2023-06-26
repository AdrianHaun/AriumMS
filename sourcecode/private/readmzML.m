function [PeaksMS1,timesMS1,varargout] = readmzML(DataPath,options)
%% reads mzML files and outputs mz, Intensity and time data
% .mzML files can be 32bit and 64bit encoded can use zlib compression, numpress compression is not yet
% supported
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
    DataPath            (1,1) string
    options.StartTime   (1,1) {mustBeNumeric} = 0
    options.EndTime     (1,1) {mustBeNumeric} = Inf
    options.MSLevel     (1,1) {mustBeInteger} = 1
end

S = readstruct(DataPath,"FileType","xml");
S = [S.mzML.run.spectrumList.spectrum];
%get encoding precision and compression
mzPrecision = S(1).binaryDataArrayList.binaryDataArray(1).cvParam(1).nameAttribute;
switch mzPrecision
    case "32-bit float"
        mzPrecision = 'single';
    case "32-bit integer"
        mzPrecision = 'single';
    case "64-bit float"
        mzPrecision = 'double';
    case "64-bit integer"
        mzPrecision = 'double';
end
IntPrecision = S(1).binaryDataArrayList.binaryDataArray(2).cvParam(1).nameAttribute;
switch IntPrecision
    case "32-bit float"
        IntPrecision = 'single';
    case "32-bit integer"
        IntPrecision = 'single';
    case "64-bit float"
        IntPrecision = 'double';
    case "64-bit integer"
        IntPrecision = 'double';
end
mzDecode=S(1).binaryDataArrayList.binaryDataArray(1).cvParam(2).nameAttribute;
switch mzDecode
    case "no compression"
        mzzlib = false;
    case "zlib compression"
        mzzlib = true;
    case "MS-Numpress linear prediction compression"
        error("Files that use numpress compression are not supported at this time.")

    case "MS-Numpress linear prediction compression followed by zlib compression"
        error("Files that use numpress compression are not supported at this time.")
end
IntDecode=S(1).binaryDataArrayList.binaryDataArray(2).cvParam(2).nameAttribute;
switch IntDecode
    case "no compression"
        Intzlib = false;
    case "zlib compression"
        Intzlib = true;
    case "MS-Numpress positive integer compression"
        error("Files that use numpress compression are not supported at this time.")

    case "MS-Numpress positive integer compression followed by zlib compression"
        error("Files that use numpress compression are not supported at this time.")
end


nScan = size(S,2);
%get times and cut to size
timeVec = zeros(nScan,1);
parfor n=1:nScan
    timeVec(n) = S(n).scanList.scan.cvParam.valueAttribute;
end
% convert to seconds
TimeScale=S(1).scanList.scan.cvParam(1).unitNameAttribute;
switch TimeScale
    case "minute"
        timeVec = timeVec*60;
end
idx=timeVec < options.StartTime | timeVec > options.EndTime;
timeVec(idx)=[];
S(idx)=[];
%gather scan data
nScan = size(S,2);
PrecursorMass = zeros(nScan,1);
msLevelVec = zeros(nScan,1);
Peaks = cell(nScan,1);
MSL=options.MSLevel;
parfor n=1:nScan
    findMSLevel=deal([S(n).cvParam.nameAttribute]);
    findMSLevel=strcmp(findMSLevel,'ms level');
    msLevelVec(n) = S(n).cvParam(findMSLevel).valueAttribute;
    if msLevelVec(n) > 1
        PrecursorMass(n) = S(n).precursorList.precursor.selectedIonList.selectedIon.cvParam(1).valueAttribute;
    else
        PrecursorMass(n) = 0;
    end
    if msLevelVec(n) > MSL
        mzbinary=0;
        Intbinary=0;
    else
        [mzbinary,Intbinary] = S(n).binaryDataArrayList.binaryDataArray.binary;
        %base64decode
        mzbinary = typecast(matlab.net.base64decode(mzbinary),'uint8');
        Intbinary = typecast(matlab.net.base64decode(Intbinary),'uint8');
        % zlib decompression
        if mzzlib == true
            mzbinary = zmat(mzbinary,0,'zlib');
        end
        if Intzlib == true
            Intbinary = zmat(Intbinary,0,'zlib');
        end
        % conversion to single or double
        mzbinary = typecast(mzbinary,mzPrecision)';
        Intbinary = typecast(Intbinary,IntPrecision)';
    end
    %store converted data
    Peaks{n} = [mzbinary Intbinary];
end
%split MsLevels
timesMS1 = timeVec(msLevelVec == 1);
PeaksMS1 = Peaks(msLevelVec == 1);
if options.MSLevel > 1 && options.MSLevel <= max(msLevelVec)
    varargout{1}=[PrecursorMass,msLevelVec];
    out=cell(options.MSLevel-1,2);
    for n = 2:options.MSLevel
        out{n-1,1} = Peaks(msLevelVec == n);
        out{n-1,2} = timeVec(msLevelVec == n);
    end
    out=reshape(out,1,[]);
    varargout = [varargout,out];
end
end
