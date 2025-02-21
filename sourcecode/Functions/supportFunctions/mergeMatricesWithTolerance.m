function out = mergeMatricesWithTolerance(Spectra,mzTol,mzTolUnit)
for n=2:size(Spectra,1)    
    mergedMatrix = [];
    matrix1 = Spectra{1};
    matrix2 = Spectra{n};
    % Iterate through each row of the first matrix
    for i = 1:size(matrix1, 1)
        mergedRow = [];
        merged = false;

        % Iterate through each row of the second matrix
        for j = 1:size(matrix2, 1)
            % Check if values in the first column are within tolerance
            switch mzTolUnit
                case "Da"
                    Check = abs(matrix1(i, 1) - matrix2(j, 1)) <= mzTol;
                case "ppm"
                    Check = abs((matrix1(i, 1)-matrix2(j, 1))/matrix2(j, 1)*10^6)<=mzTol;
            end
            if Check == true 
                % Average the values in the first column
                averageValue = (matrix1(i, 1) + matrix2(j, 1)) / 2;

                % Add the averaged value to the merged row
                mergedRow = [averageValue, matrix1(i, 2:end), matrix2(j, 2:end)];
                
                merged = true;
                break; % Exit inner loop once a match is found
            end
        end

        % If no match is found, add the row from the first matrix
        if ~merged
            mergedRow = [matrix1(i, :), zeros(1, size(matrix2, 2)-1)];
        end

        % Add the merged row to the result matrix
        mergedMatrix = [mergedMatrix; mergedRow];
    end

    % Add remaining rows from the second matrix
    for j = 1:size(matrix2, 1)
        switch mzTolUnit
                case "Da"
                    Check = abs(mergedMatrix(:, 1) - matrix2(j, 1)) <= mzTol;
                case "ppm"
                    Check = abs((mergedMatrix(:, 1)-matrix2(j, 1))/matrix2(j, 1)*10^6)<=mzTol;
        end
        if ~any(Check)
            % Add the row from the second matrix
            mergedMatrix = [mergedMatrix; [matrix2(j, 1), zeros(1, size(matrix1, 2)-1), matrix2(j, 2:end)]];
        end
    end
    Spectra{1} = mergedMatrix;
end
out = Spectra{1};
end
