function FileInfo = mzMLinfo(FileName)
%Reads basic file information from mzML files, only mzML files without
%numPress compression are supported.
%Output Fileinfo is a struct with fields: FileName, FileModDate, FileSize,
%NumberOf Scans, StartTime, EndTime. Start and End times are in seconds

arguments
    FileName(1,:) {mustBeFile}
end

FileInfo = struct('Filename',[],...
    'FileModDate',[],...
    'FileSize',[],...
    'NumberOfScans','N/A',...
    'StartTime','N/A',...
    'EndTime','N/A');

% Get filename, filesize and date
fInfo = dir(FileName);
FileInfo.Filename = fInfo.name;
FileInfo.FileModDate = fInfo.date;
FileInfo.FileSize = fInfo.bytes;

%Parse mzXML file
import matlab.io.xml.dom.*
d = parseFile(Parser, FileName);
%get chromatogram Nodes
chromatogramList = d.getElementsByTagName('chromatogramList').item(0);
cvParams = chromatogramList.getElementsByTagName('cvParam');
%get encoding precision
cvParam = cvParams.item(1);
precision = cvParam.getAttribute('name');
switch precision
    case "32-bit float"
        precision = 'single';
    case "32-bit integer"
        precision = 'single';
    case "64-bit float"
        precision = 'double';
    case "64-bit integer"
        precision = 'double';
end
%get compression type
cvParam = cvParams.item(2);
compression = cvParam.getAttribute('name');
switch compression
    case "no compression"
        compression = false;
    case "zlib compression"
        compression = true;
    case "MS-Numpress positive integer compression"
        error("Files that use numpress compression are not supported at this time.")
    case "MS-Numpress positive integer compression followed by zlib compression"
        error("Files that use numpress compression are not supported at this time.")
end
%get binary data
cvParam = cvParams.item(3);
unit = cvParam.getAttribute('unitName');
TimeArrayNode = cvParam.getParentNode();
binaryNode = TimeArrayNode.getElementsByTagName('binary').item(0);
TimeBinary =binaryNode.TextContent;
%decode binary
TimeBinary = typecast(matlab.net.base64decode(TimeBinary),'uint8');
if compression == true
    TimeBinary = zmat(TimeBinary,0,'zlib');
end
TimeBinary = typecast(TimeBinary,precision)';
%conversion to seconds
if strcmp(unit,'minute')
    TimeBinary = TimeBinary.*60;
end
FileInfo.StartTime = min(TimeBinary);
FileInfo.EndTime = max(TimeBinary);
FileInfo.NumberOfScans = length(TimeBinary);
end