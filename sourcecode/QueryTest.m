tolerance = 0.1;

out = cell(length(Masses),1);

for n=1:length(Masses)

min = Masses(n)-tolerance;
max = Masses(n)+tolerance;


query = ['SELECT DISTINCT NAME, ' ...
    '	FORMULA, ' ...
    '	CAS_NUMBER, ' ...
    '	EXACT_MASS, ' ...
    '	PRECURSOR_TYPE ' ...
    'FROM MassBankData ' ...
    'WHERE EXACT_MASS <= '];


query = append(query,convertStringsToChars(max + "	AND EXACT_MASS >= " + min));

result = fetch(app.DataBase.DataBaseConnection, query);

%calculate difference in ppm
diff = ((abs(Masses(n)-result.EXACT_MASS))./Masses(n))*10^6;

result.DELTA = diff;
result = sortrows(result,"DELTA","ascend");

out{n}=result;
end

outPPM = cell(length(Masses),1);
tolerance = 100;
for n=1:length(Masses)

min = Masses(n)-(Masses(n)*tolerance*10^-6);
max = Masses(n)+(Masses(n)*tolerance*10^-6);


query = ['SELECT DISTINCT NAME, ' ...
    '	FORMULA, ' ...
    '	CAS_NUMBER, ' ...
    '	EXACT_MASS, ' ...
    '	PRECURSOR_TYPE ' ...
    'FROM MassBankData ' ...
    'WHERE EXACT_MASS <= '];


query = append(query,convertStringsToChars(max + "	AND EXACT_MASS >= " + min));

result = fetch(app.DataBase.DataBaseConnection, query);

%calculate difference in ppm
diff = ((abs(Masses(n)-result.EXACT_MASS))./Masses(n))*10^6;

result.DELTA = diff;
result = sortrows(result,"DELTA","ascend");

outPPM{n}=result;
end