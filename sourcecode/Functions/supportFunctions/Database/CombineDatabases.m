[file,path] = uigetfile('*.db',MultiSelect="on");
databases = fullfile(path,file);
FinalDB = "D:\MSData\MS_Databases\HMDB_full.db";

FinalCon = sqlite(FinalDB,"connect");
for n = size(databases,2)
    CurrentDBCon = sqlite(databases{1,2},"readonly");
    %Set query to execute on the database
    query = ['SELECT * ' ...
    'FROM SpectralData'];

    %% Execute query and fetch results
    data = fetch(CurrentDBCon,query);
    close(CurrentDBCon)
    sqlwrite(FinalCon,"SpectralData",data)

end
close(FinalCon)