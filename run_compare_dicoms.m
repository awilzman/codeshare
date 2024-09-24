clear %save data before starting!!!
directory = ''; %end slash is important
%resolution in micrometers
res = 82;
%Ensure slope and intercept are correct!
calibrate_slope = 0.000357;
calibrate_int = -0.0012625; 
excelFileName = 'output.xlsx';
allData = {};
columns = ["Anterior","Posterior","Medial","Lateral"];
rows = ["";"Scan 1";"Scan 2";"Difference"];

LCV_name = 'lcv';
mask1_name = 'scan1';
mask2_name = 'scan2';
labels = ["File name";
    "Total Volume [cm^3]";
    "Bone Volume [cm^3]";
    "Bone Mineral Content [g]";
    "Bone Mineral Density [g/cm^3]"];
labels = cellstr(labels);
    
%list all subjects as comma separated list
subs = ["15"];
medial_left = 0;
for s=1:length(subs)
    subjects = [strcat(subs(s)," 4 "),strcat(subs(s)," 30 ")];
    angle = 0;
    for i=1:length(subjects)
        output{i} = strcat(subjects(i),LCV_name);
        mask1s{i} = strcat(subjects(i),mask1_name);
        mask2s{i} = strcat(subjects(i),mask2_name);
    end
    
    for i=1:length(output)
        disp(strcat("Running ",subjects(i)))
        [output{2,i},output{3,i},output{4,i},output{5,i},medial_left,angle] = compare_dicoms(directory, ...
            res,output{1,i} ,mask1s{i},mask2s{i},calibrate_slope,calibrate_int,medial_left,angle);
        for j=2:5
            output{j,i} = cat(1,columns,output{j,i});
            output{j,i} = cat(2,rows,output{j,i});
            output{j,i}(1,1) = labels(j);
        end
    end
    output = [labels,output];
    
    for j = 2:length(subjects) + 1
        baseRow = 36 * (s - 1) + (j - 2) * 18 + 1; % Adjust baseRow to ensure unique placement
        allData{baseRow, 1} = output{1, j}; % Set the filename
    
        for i = 2:5
            startRow = baseRow + (i - 2) * 4; % Starting row for current data
            dataMatrix = output{i, j}; % Current data matrix
            [numDataRows, numDataCols] = size(dataMatrix); % Get size of the data matrix
    
            % Fill allData with data from the current output
            for r = 1:numDataRows
                for c = 1:numDataCols
                    allData{startRow + r - 1, c + 1} = dataMatrix(r, c); % Fill in data
                end
            end
        end
    end
end

writecell(allData, excelFileName, 'Sheet', 1);