clear;  % Clear workspace

% Initial Setup
subs = ["14","80","81","78","28"];  % List of subjects
directory = '';  % Directory of DICOM files (end slash is important)
res = 82;  % Image resolution (in micrometers)
calibrate_slope = 0.000379;
calibrate_int = -0.00603;
excelFileName = 'output.xlsx';  % Output file name

columns = ["Anterior", "Posterior", "Medial", "Lateral"]; 
rows = [""; "Scan 1"; "Scan 2"; "Difference"];
labels = cellstr(["File name"; "Total Volume [cm^3]"; "Bone Volume [cm^3]"; ...
    "Bone Mineral Content [g]"; "Bone Mineral Density [g/cm^3]"]);

% Initialize
output = cell(length(subs), 2, 5);  % Preallocate for output
mask1s = cell(length(subs), 2);  % Preallocate for mask1s
mask2s = cell(length(subs), 2);  % Preallocate for mask2s

% Loop over subjects
for s = 1:length(subs)
    subjects = strcat(subs(s), [" 4 ", " 30 "]);
    angle_rot = 0;
    medial_left = 0;
    
    % Prepare file names and compare DICOMs
    for i = 1:length(subjects)
        output{s, i, 1} = strcat(subjects(i), 'lcv');  % Scan name
        mask1s{s, i} = strcat(subjects(i), 'scan1');   % Mask 1 name
        mask2s{s, i} = strcat(subjects(i), 'scan2');   % Mask 2 name
         
        disp("Running " + subjects(i));  % Log current subject
        
        % Compare DICOMs and store results
        try
            % Each feature is now a 3x4 table
            [tv, bv, bmc, bmd, medial_left, angle_rot] = ...
                compare_dicoms(directory, res, output{s, i, 1}, mask1s{s, i}, mask2s{s, i}, ...
                calibrate_slope, calibrate_int, medial_left, angle_rot);
            
            % Store 3x4 tables for each feature
            output{s, i, 2} = tv;   % 3x4 table for tv
            output{s, i, 3} = bv;   % 3x4 table for bv
            output{s, i, 4} = bmc;  % 3x4 table for bmc
            output{s, i, 5} = bmd;  % 3x4 table for bmd
            
        catch
            disp(['No data for ' + subjects(i)]);
        end
    end
end

% Initialize data array for output to Excel
% Calculate total number of rows and columns
numSubjects = length(subs);
numScans = size(output, 2);
totalRows = numSubjects * numScans * 36; 
totalCols = 5; % Maximum number of columns used

% Preallocate 'data' cell array
data = cell(totalRows, totalCols);

% Save results to Excel
for s = 1:length(subs)  % Loop over subjects
    for i = 1:size(output, 2) 
        baseRow = 36 * (s - 1) + 18 * (i - 1) + 1; 
        
        % Store subject info (filename of the scan)
        if ~isempty(output{s, i, 1})
            data{baseRow + 1, 1} = output{s, i, 1};  
        else
            data{baseRow + 1, 1} = 'No Data';
        end
        
        % Create labels for each row
        data{baseRow + 1, 1} = 'Scan 1';
        data{baseRow + 2, 1} = 'Scan 2';
        data{baseRow + 3, 1} = 'Difference';
        
        % Labels for the columns (Anterior, Posterior, Medial, Lateral)
        data(baseRow, 2:5) = {'Anterior', 'Posterior', 'Medial', 'Lateral'};
        
        % Loop through each feature (tv, bv, bmc, bmd)
        for j = 2:5
            startRow = baseRow + 4 * (j - 2);
            
            featureLabel = {'Total Volume [cm^3]', 'Bone Volume [cm^3]', ... 
                'Bone Mineral Content [g]', 'Bone Mineral Density [g/cm^3]'};  
            
            % Label for the feature
            data{startRow, 1} = strcat(subs{s},'_',featureLabel{j - 1});  
            data{startRow + 1, 1} = 'Scan 1';
            data{startRow + 2, 1} = 'Scan 2';
            data{startRow + 3, 1} = 'Difference';
            
            % Check if output exists for this subject, scan, and feature
            if ~isempty(output{s, i, j}) && size(output{s, i, j}, 1) >= 3 && size(output{s, i, j}, 2) >= 4
                dataMatrix = output{s, i, j}; 
                for r = 1:3  % Row 1 = Scan1, Row 2 = Scan2, Row 3 = Difference
                    data(startRow + r, 2:5) = num2cell(dataMatrix(r, :));
                end
            else
                % If no data is available, add placeholders
                for r = 1:3
                    data(startRow + r, 2:5) = {'No Data', 'No Data', 'No Data', 'No Data'};
                end
            end
        end
    end
end


% Determine a new file name with incremented version number
fileNumber = 1;
[folder, name, ext] = fileparts(excelFileName);
newFileName = fullfile(folder, [name, '_v', num2str(fileNumber), ext]);

while isfile(newFileName)
    fileNumber = fileNumber + 1;
    newFileName = fullfile(folder, [name, '_v', num2str(fileNumber), ext]);
end

% Write the data to Excel as a new file
writecell(data, newFileName, 'Sheet', 1);
