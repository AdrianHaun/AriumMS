function DecodedMS2Data = DecodeStrings(EncodedStrings)
DecodedMS2Data=cell(size(EncodedStrings,1),1);
parfor n=1:size(EncodedStrings,1)
    Decoded = matlab.net.base64decode(EncodedStrings(n));
    Decoded = typecast(Decoded,'double');
    DecodedMS2Data{n} = reshape(Decoded,[],2);
end