clear %save data before starting!!!
directory = ''; %end slash is important
%resolution in micrometers
res = 82;
%Ensure slope and intercept are correct!
calibrate_slope = 0.000357;
calibrate_int = -0.0012625; 

%list all subjects as comma separated list
subjects = ["15 4 "]; %eg ["1","2","3"]
LCV_name = 'lcv';
mask1_name = 'scan1';
mask2_name = 'scan2';
labels = ["File name";"Total Volume [cm^3]";"Bone Volume [cm^3]";"Bone Mineral Content [g]";"Bone Mineral Density [g/cm^3]"];
labels = cellstr(labels);

% this code requires each lcv, mask1, and mask2 to be named the same with 
% an incrementing suffix. this could be changed to feed from a list of file
% names
for i=1:length(subjects)
    output{i} = strcat(subjects(i),LCV_name);
    mask1s{i} = strcat(subjects(i),mask1_name);
    mask2s{i} = strcat(subjects(i),mask2_name);
end
columns = ["Anterior","Posterior","Medial","Lateral"];
rows = ["";"Scan 1";"Scan 2";"Difference"];
for i=1:length(output)
    [output{2,i},output{3,i},output{4,i},output{5,i}] = compare_dicoms(directory, ...
        res,output{1,i} ,mask1s{i},mask2s{i},calibrate_slope,calibrate_int);
    for j=2:5
        output{j,i} = cat(1,columns,output{j,i});
        output{j,i} = cat(2,rows,output{j,i});
        output{j,i}(1,1) = labels(j);
    end
end
output = [labels,output];
excelFileName = 'output.xlsx';
allData = {};

for j = 2:length(subjects) + 1
    baseRow = 18 * (j - 2) + 3;
    allData{baseRow, 1} = output{1, 2};
    for i = 2:5
        startRow = baseRow + (i - 1) * 4;
        label = output{i, 1};
        dataMatrix = output{i, 2};
        numDataRows = size(dataMatrix, 1);
        numDataCols = size(dataMatrix, 2);
        allData{startRow, 1} = label;
        for r = 1:numDataRows
            for c = 1:numDataCols
                allData{startRow + r - 1, 1 + c} = dataMatrix(r, c);
            end
        end
    end
end
writecell(allData, excelFileName, 'Sheet', 1);