function FileInfo = mzMLinfo(FileName)
%Reads basic file information from mzML files, only mzML files without
%numPress compression are supported.
%Output Fileinfo is a struct with fields: FileName, FileModDate, FileSize,
%NumberOf Scans, StartTime, EndTime. Start and End times are in seconds

% BSD 3-Clause License
% 
% Copyright (c) 2022, Adrian Haun
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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