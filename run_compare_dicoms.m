clear %save data before starting!!!
directory = ''; %end slash is important
%resolution in micrometers
res = 82;
%Ensure slope and intercept are correct!
calibrate_slope = 0.000357;
calibrate_int = -0.0012625; 

%list all subjects as comma separated list
subjects = ["15 4 ","28 4 "]; %eg ["1","2","3"]
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
ap_x = [1,0];
ap_y = [0,1];
for i=1:length(output)
    [output{2,i},output{3,i},output{4,i},output{5,i},ap_x,ap_y] = compare_dicoms(directory, ...
        res,output{1,i} ,mask1s{i},mask2s{i},calibrate_slope,calibrate_int,subjects(i),ap_x,ap_y);
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
    baseRow = 18 * (j - 2) + 1;
    allData{baseRow, 1} = output{1, j};
    for i = 2:5
        startRow = baseRow + (i - 2) * 4;
        label = output{1, j};
        dataMatrix = output{i, j};
        numDataRows = size(dataMatrix, 1);
        numDataCols = size(dataMatrix, 2);
        
        for r = 1:numDataRows
            for c = 1:numDataCols
                allData{startRow + r - 1, 1 + c} = dataMatrix(r, c);
            end
        end
    end
end
writecell(allData, excelFileName, 'Sheet', 1);