function database = generateAMDISDatabase(AMDISfile)

arguments
    AMDISfile (1,1) string {mustBeFile}
end


% Create an empty cell array to store your data
database = struct('ACCESSION', {}, 'AUTHORS', {}, 'NAME', {}, 'FORMULA', {}, ...
    'CAS_NUMBER', {}, 'EXACT_MASS', {}, 'INSTRUMENT_TYPE', {},'FRAGMENTATION_TYPE', {}, ...
    'PRECURSOR_TYPE', {}, 'PRECURSOR_MZ', {}, 'SPECTRUM', {});

%read file
entryCells = readLibraryAMDIS(AMDISfile);

% Loop through the fileList and read the text files
parfor iEntry = 1:height(entryCells)
    currentEntry = entryCells{iEntry,1};
    % Add the extracted data to the database structure
    entry = struct('ACCESSION', "", 'AUTHORS', "", 'NAME', currentEntry(1,2), ...
        'FORMULA', currentEntry(2,2), 'CAS_NUMBER', currentEntry(3,2), 'EXACT_MASS', "", ...
        'INSTRUMENT_TYPE', "GC-MS", 'FRAGMENTATION_TYPE', "EI", 'PRECURSOR_TYPE', "MaxIntensityFragment", ...
        'PRECURSOR_MZ', currentEntry(5,2), 'SPECTRUM', currentEntry(6,2));

    database = [database, entry];
end
database=struct2table(database);
conn = sqlite('AMDIS.db',"create");
sqlwrite(conn,"SpectralDataBase",database);
close(conn)
end
