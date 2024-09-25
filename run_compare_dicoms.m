clear % save data before starting

% Initial Setup
directory = ''; % Directory of DICOM files (end slash is important)
res = 82; % Image resolution (in micrometers)
calibrate_slope = 0.000357;
calibrate_int = -0.0012625;

excelFileName = 'output.xlsx'; % Output file name
allData = {}; 
columns = ["Anterior", "Posterior", "Medial", "Lateral"]; 
rows = [""; "Scan 1"; "Scan 2"; "Difference"];

% Filenames
LCV_name = 'lcv'; mask1_name = 'scan1'; mask2_name = 'scan2';

% Labels for output
labels = cellstr(["File name"; "Total Volume [cm^3]"; "Bone Volume [cm^3]"; "Bone Mineral Content [g]"; "Bone Mineral Density [g/cm^3]"]);

subs = ["15"]; % List of subjects
medial_left = 0;

% Loop through subjects
for s = 1:length(subs)
    subjects = [strcat(subs(s), " 4 "), strcat(subs(s), " 30 ")];
    angle = 0;
    
    % Prepare file names
    for i = 1:length(subjects)
        output{s,i} = strcat(subjects(i), LCV_name);
        mask1s{s,i} = strcat(subjects(i), mask1_name);
        mask2s{s,i} = strcat(subjects(i), mask2_name);
    end
    
    % Compare DICOMs and store results
    for i = 1:length(output)
        disp(strcat("Running ", subjects(i)))
        [output{s,2,i}, output{s,3,i}, output{s,4,i}, output{s,5,i}, medial_left, angle] = compare_dicoms(directory, ...
            res, output{s,1,i}, mask1s{i}, mask2s{i}, calibrate_slope, calibrate_int, medial_left, angle);
        
        % Add headers to output data
        for j = 2:5
            output{s,j,i} = cat(1, columns, output{s,j,i});
            output{s,j,i} = cat(2, rows, output{s,j,i});
            output{s,j,i}(1,1) = labels(j);
        end
    end
    
    % Append subject's data into output
    output{s} = [labels, output];
end

% Save the results to the Excel file
for s = 1:length(subs)
    for j = 2:length(subjects) + 1
        baseRow = 36 * (s - 1) + (j - 2) * 18 + 1; 
        allData{baseRow, 1} = output{1, j}; 
        
        for i = 2:5
            startRow = baseRow + (i - 2) * 4; 
            dataMatrix = output{i, j}; 
            [numDataRows, numDataCols] = size(dataMatrix); 
            
            % Insert data into allData
            for r = 1:numDataRows
                for c = 1:numDataCols
                    allData{startRow + r - 1, c + 1} = dataMatrix(r, c); 
                end
            end
        end
    end
end

writecell(allData, excelFileName, 'Sheet', 1);