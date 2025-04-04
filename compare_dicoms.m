% Read DICOM stacks of LCV, Time1, and Time2, and raw image of time1
% raw image time1 is reference for registering time2 earlier
% Rotate scans on user-defined AP axis and
% calculate intensity [unit] differences between t1 and t2 of
% bone volume (bv) [cm^3], bone mineral content (bmc) [g], and 
% bone mineral density [g/cm^3] (bmd) 
%               % Anterior          % Posterior         % Medial           % Lateral
% Scan 1        % (1,1)             % (1,2)             % (1,3)            % (1,4)
% Scan 2        % (2,1)             % (2,2)             % (2,3)            % (2,4)
% Difference    % (3,1)             % (3,2)             % (3,3)            % (3,4)
%
% Written by Andrew Wilzman and Karen Troy 06/2023
% Updated 12/13/2024
% Example run:
% calibrate_slope = 0.00035619;
% calibrate_int = -0.00365584; 
%res in um, voxel edge length
% [bv,bmc,bmd,medial_left,angle_rot] = compare_dicoms(default_directory,res,LCV_name,
% mask1_name,mask2_name,calibrate_slope,calibrate_int,
% medial_left,angle_rot,first_full_slice)

% T. Hildebrand, A. Laib, R. Müller, J. Dequecker, P. Rüegsegger. 
% Direct 3-D morphometric analysis of human cancellous bone: 
% microstructural data from spine, femur, iliac crest and calcaneus. 
% J Bone Miner Res 1999;14(7):1167-74.
function [tv, bv, bmc, bmd, medial_left, angle_rot] = compare_dicoms(default_directory,res, ...
    mask1_name,mask2_name,calibrate_slope,calibrate_int, ...
    medial_left,angle_rot)

    if nargin < 9
        angle_rot = 0;
    end
    if nargin < 8
        medial_left = 0;
    end

    %mask_LCV = get_mask(LCV_name,default_directory);
    mask_1 = get_mask(mask1_name,default_directory);
    mask_2 = get_mask(mask2_name,default_directory);
    
    % un-Pad matrices
    mask_1 = pad_3dmat(mask_1);
    mask_2 = pad_3dmat(mask_2);

    tangle = pi()/4;

    bone_threshold = 0.5/calibrate_slope;
    mask_1(mask_1 < bone_threshold) = 0;
    mask_2(mask_2 < bone_threshold) = 0;

    if angle_rot==0
        %threshold image
        %classify tibia and fibula by size (larger = tibia)
        %define centers of each as above
        %calculate the angle of the line with respect to the X axis CCW+
        dicom_files = dir(fullfile(strcat(chdir, '\', mask1_name, '_raw\'), '*.dcm'));
        dir_path = strcat(chdir,'\',mask1_name,'_raw\');
        % Check if any DICOM files exist and select the first one
        if ~isempty(dicom_files)
            d_name = dicom_files(1).name;  % First DICOM file
            full_dicom_path = fullfile(dir_path, d_name);  % Full path to DICOM file
        else
            error('No DICOM files found in the directory.');
        end
        
        raw_image = dicomread(full_dicom_path);
        info = dicominfo(full_dicom_path);
        raw_image = raw_image * info.RescaleSlope + info.RescaleIntercept;

        raw_image_b = raw_image > bone_threshold;
        raw_image_b = logical(raw_image_b);
        stats = regionprops(raw_image_b, 'Area', 'Centroid');
        min_area_threshold = 5000;
        stats = stats([stats.Area] > min_area_threshold);
        areas = [stats.Area];
        [~, sorted_idx] = sort(areas, 'descend');
        idx_tibia = sorted_idx(1);
        idx_fibula = sorted_idx(2);

        tibia_centroid = stats(idx_tibia).Centroid; 
        fibula_centroid = stats(idx_fibula).Centroid;

        delta_y = fibula_centroid(2) - tibia_centroid(2);
        delta_x = fibula_centroid(1) - tibia_centroid(1);
        % Angle calculation
        
        % Medial side calculation
        if tibia_centroid(1) < fibula_centroid(1)
            medial_left = 1;
            angle_rot = atan2(delta_y, delta_x)-tangle;
        else
            medial_left = 0;
            angle_rot = atan2(delta_y, delta_x)-pi()-tangle;
        end
    end

    z_center = round(size(mask_1, 3) / 2); % Middle slice index
    center_slice = mask_1(:, :, z_center); % 2D slice along Z axis
    non_zero_voxels_2D = center_slice > 0;
    filled_slice = imfill(non_zero_voxels_2D, 'holes');
    [rows, cols] = find(filled_slice); 
    
    x = mean(cols);
    y = mean(rows);
    ccenter = [x, y];
    % Rotate mask to set anterior in quadrant 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                    #
    %                    #
    %                    #
    % Medial/Lateral     #  Anterior
    % ########################################
    %         Posterior  #  Medial/Lateral
    %                    #
    %                    #
    %                    #
    %                    #
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mask_1 = rotate_mask(mask_1,angle_rot,ccenter);
    mask_2 = rotate_mask(mask_2,angle_rot,ccenter);
    
    mask_1 = pad_3dmat(mask_1);
    mask_2 = pad_3dmat(mask_2);

    % Calculate metrics
    tv = zeros(3,4);
    bv = zeros(3,4);
    bmc = zeros(3,4);
    bmd = zeros(3,4);
    
    % Initialize masks
    mask_1_med = zeros(size(mask_1));
    mask_1_ant = zeros(size(mask_1));
    mask_1_post = zeros(size(mask_1));
    mask_1_lat = zeros(size(mask_1));
    
    mask_2_med = zeros(size(mask_2));
    mask_2_ant = zeros(size(mask_2));
    mask_2_post = zeros(size(mask_2));
    mask_2_lat = zeros(size(mask_2));
    
    % Iterate through slices

    med_s = 4;
    lat_s = 1;

    if medial_left
        med_s = 1;
        lat_s = 4;
    end

    for z = 1:size(mask_1, 3)
        % Generate binary wedge masks for the current slice
        wedge_masks_1 = generateWedgeMasks(mask_1(:,:,z));
        wedge_masks_2 = generateWedgeMasks(mask_2(:,:,z));

        mask_1_slice = double(mask_1(:,:,z));
        mask_2_slice = double(mask_2(:,:,z));
    
        % Apply the binary wedge masks to extract the original values
        mask_1_med(:,:,z) = max(mask_1_slice .* double(wedge_masks_1(:,:,med_s)), double(wedge_masks_1(:,:,med_s)));
        mask_1_ant(:,:,z) = max(mask_1_slice .* double(wedge_masks_1(:,:,2)), double(wedge_masks_1(:,:,2)));
        mask_1_post(:,:,z) = max(mask_1_slice .* double(wedge_masks_1(:,:,3)), double(wedge_masks_1(:,:,3)));
        mask_1_lat(:,:,z) = max(mask_1_slice .* double(wedge_masks_1(:,:,lat_s)), double(wedge_masks_1(:,:,lat_s)));
    
        mask_2_med(:,:,z) = max(mask_2_slice .* double(wedge_masks_2(:,:,med_s)), double(wedge_masks_2(:,:,med_s)));
        mask_2_ant(:,:,z) = max(mask_2_slice .* double(wedge_masks_2(:,:,2)), double(wedge_masks_2(:,:,2)));
        mask_2_post(:,:,z) = max(mask_2_slice .* double(wedge_masks_2(:,:,3)), double(wedge_masks_2(:,:,3)));
        mask_2_lat(:,:,z) = max(mask_2_slice .* double(wedge_masks_2(:,:,lat_s)), double(wedge_masks_2(:,:,lat_s)));
    end

    [tv(1,1), bv(1,1), bmc(1,1), bmd(1,1)] = bv_bmc(mask_1_ant,res,calibrate_slope,calibrate_int);
    [tv(1,2), bv(1,2), bmc(1,2), bmd(1,2)] = bv_bmc(mask_1_post,res,calibrate_slope,calibrate_int);
    [tv(1,3), bv(1,3), bmc(1,3), bmd(1,3)] = bv_bmc(mask_1_med,res,calibrate_slope,calibrate_int);
    [tv(1,4), bv(1,4), bmc(1,4), bmd(1,4)] = bv_bmc(mask_1_lat,res,calibrate_slope,calibrate_int);
    [tv(2,1), bv(2,1), bmc(2,1), bmd(2,1)] = bv_bmc(mask_2_ant,res,calibrate_slope,calibrate_int);
    [tv(2,2), bv(2,2), bmc(2,2), bmd(2,2)] = bv_bmc(mask_2_post,res,calibrate_slope,calibrate_int);
    [tv(2,3), bv(2,3), bmc(2,3), bmd(2,3)] = bv_bmc(mask_2_med,res,calibrate_slope,calibrate_int);
    [tv(2,4), bv(2,4), bmc(2,4), bmd(2,4)] = bv_bmc(mask_2_lat,res,calibrate_slope,calibrate_int);
    tv(3,1) = tv(2,1)-tv(1,1);
    tv(3,2) = tv(2,2)-tv(1,2);
    tv(3,3) = tv(2,3)-tv(1,3);
    tv(3,4) = tv(2,4)-tv(1,4);
    bv(3,1) = bv(2,1)-bv(1,1);
    bv(3,2) = bv(2,2)-bv(1,2);
    bv(3,3) = bv(2,3)-bv(1,3);
    bv(3,4) = bv(2,4)-bv(1,4);
    bmc(3,1) = bmc(2,1)-bmc(1,1);
    bmc(3,2) = bmc(2,2)-bmc(1,2);
    bmc(3,3) = bmc(2,3)-bmc(1,3);
    bmc(3,4) = bmc(2,4)-bmc(1,4);
    bmd(3,1) = bmd(2,1)-bmd(1,1);
    bmd(3,2) = bmd(2,2)-bmd(1,2);
    bmd(3,3) = bmd(2,3)-bmd(1,3);
    bmd(3,4) = bmd(2,4)-bmd(1,4);
    
    z_center = round(size(mask_1, 3) / 2); % Middle slice index
    center_slice = mask_1(:, :, z_center); % 2D slice along Z axis
    non_zero_voxels_2D = center_slice > 0;
    filled_slice = imfill(non_zero_voxels_2D, 'holes');
    [rows, cols] = find(filled_slice); 

    x = mean(cols);
    y = mean(rows);

    slice_image = mask_1(:,:,round(size(mask_1,3)/2));
    [rows, cols] = size(slice_image);
    [xGrid, yGrid] = meshgrid(1:cols, 1:rows);
    ap_slope = tan(tangle);

    line_mask = abs(yGrid - (ap_slope * (xGrid - x) + y)) < 1;
    slice_image(line_mask) = max(slice_image(:));

    line_mask = abs(yGrid - ((1/(-ap_slope + 1e-6)) * (xGrid - x) + y)) < 1;
    slice_image(line_mask) = max(slice_image(:));

    slice_image = uint16(slice_image);
    imwrite(slice_image,strcat(default_directory,mask1_name,".png"))

    color_ant = [255, 0, 0];   % Red for anterior
    color_post = [0, 255, 0];  % Green for posterior
    color_med = [0, 0, 255];   % Blue for medial
    color_lat = [255, 255, 0]; % Yellow for lateral
    slice_image_rgb = zeros([size(slice_image), 3], 'uint8');
    
    for z = 1:size(mask_1, 3)
        % Anterior section
        slice_image_rgb(:,:,1) = slice_image_rgb(:,:,1) + uint8(mask_1_ant(:,:,z)) * color_ant(1);
        slice_image_rgb(:,:,2) = slice_image_rgb(:,:,2) + uint8(mask_1_ant(:,:,z)) * color_ant(2);
        slice_image_rgb(:,:,3) = slice_image_rgb(:,:,3) + uint8(mask_1_ant(:,:,z)) * color_ant(3);
    
        % Posterior section
        slice_image_rgb(:,:,1) = slice_image_rgb(:,:,1) + uint8(mask_1_post(:,:,z)) * color_post(1);
        slice_image_rgb(:,:,2) = slice_image_rgb(:,:,2) + uint8(mask_1_post(:,:,z)) * color_post(2);
        slice_image_rgb(:,:,3) = slice_image_rgb(:,:,3) + uint8(mask_1_post(:,:,z)) * color_post(3);
    
        % Medial section
        slice_image_rgb(:,:,1) = slice_image_rgb(:,:,1) + uint8(mask_1_med(:,:,z)) * color_med(1);
        slice_image_rgb(:,:,2) = slice_image_rgb(:,:,2) + uint8(mask_1_med(:,:,z)) * color_med(2);
        slice_image_rgb(:,:,3) = slice_image_rgb(:,:,3) + uint8(mask_1_med(:,:,z)) * color_med(3);
    
        % Lateral section
        slice_image_rgb(:,:,1) = slice_image_rgb(:,:,1) + uint8(mask_1_lat(:,:,z)) * color_lat(1);
        slice_image_rgb(:,:,2) = slice_image_rgb(:,:,2) + uint8(mask_1_lat(:,:,z)) * color_lat(2);
        slice_image_rgb(:,:,3) = slice_image_rgb(:,:,3) + uint8(mask_1_lat(:,:,z)) * color_lat(3);
    end
    
    imwrite(slice_image_rgb, strcat(default_directory, mask1_name, '_segmented.png'));
    % Red anterior
    % Green Posterior
    % Blue Medial
    % Yellow Lateral
end

% Function definitions
function totalArea = calculateFilledArea(array)
    % Convert array to binary mask
    binary_image = array > 0;

    % Perform morphological closing to fill small gaps
    se = strel('disk', 5);  % Adjust structuring element as needed
    closed_image = imclose(binary_image, se);

    % Calculate total filled area directly
    totalArea = sum(closed_image(:)); % Counts nonzero pixels (alternative to polyarea)
end


function [tv, bv, bmc, bmd] = bv_bmc(mask, res, slope, int)
    tv = 0;    
    bv = 0;
    bmc = 0;

    vox_ed = res / 10000.0; % Convert um to cm

    for z = 1:size(mask, 3)
        slice = double(mask(:, :, z));
        area = calculateFilledArea(slice);

        % Skip empty slices
        if area == 0
            continue;
        end

        vol = area * vox_ed^3;
        bv_vol = nnz(slice > 1) * vox_ed^3;

        if bv_vol == 0
            continue;
        end

        slice_density = mean(slice(slice > 1), 'all') * slope + int;
        slice_content = slice_density * bv_vol;

        tv = tv + vol; % Total Volume
        bv = bv + bv_vol; % Bone Volume
        bmc = bmc + slice_content;
    end

    if tv > 0
        bmd = bmc / tv;
    else
        bmd = 0;
    end
end

function [new_mask] = rotate_mask(mask, rotationAngle, centroid, interpolationMethod)
    % Input validation
    assert(isnumeric(mask) && ndims(mask) >= 2, 'Input mask must be a numeric 2D or 3D array.');
    assert(isscalar(rotationAngle), 'Rotation angle must be a scalar.');
    
    % Set default interpolation method
    if nargin < 4
        interpolationMethod = 'linear';
    end
    % Precompute values
    max_length = round(1.2 * sqrt(size(mask, 1)^2 + size(mask, 2)^2));
    pad_LR = round((max_length - size(mask, 2)) / 2);
    pad_TB = round((max_length - size(mask, 1)) / 2);
    % Pad the mask
    mask = padarray(mask, [pad_TB, pad_LR], 0, 'both');
    % Initialize the result
    new_mask = zeros(size(mask));    
    centroid(1) = centroid(1) + pad_TB;
    centroid(2) = centroid(2) + pad_LR;
    for channelIndex = 1:size(mask, 3)
        matrix = double(mask(:, :, channelIndex));
        [rows, cols] = size(matrix);        
        [x, y] = meshgrid(1:cols, 1:rows);
        x = x - centroid(1);
        y = y - centroid(2);        
        x_rot = x * cos(rotationAngle) - y * sin(rotationAngle);
        y_rot = x * sin(rotationAngle) + y * cos(rotationAngle);        
        % Perform interpolation with specified method
        new_mask(:, :, channelIndex) = interp2(x, y, matrix, x_rot, y_rot, interpolationMethod);
    end
end

function red_mask = pad_3dmat(mask)
    % Find non-zero voxel indices along each dimension
    mask(isnan(mask)) = 0;
    [x, y, z] = ind2sub(size(mask), find(mask));
    
    % Compute bounds
    x1 = min(x); x2 = max(x);
    y1 = min(y); y2 = max(y);
    z1 = min(z); z2 = max(z);
    
    % Crop to bounding box
    red_mask = mask(x1:x2, y1:y2, z1:z2);
end

function [mask] = get_mask(name, default_directory)
    chdir = strcat(default_directory, name);
    files = dir(fullfile(chdir, '*.DCM*'));
    slices = length(files);

    % Read first slice to get dimensions
    first_slice = dicomread(fullfile(chdir, files(1).name));
    info = dicominfo(fullfile(chdir, files(1).name));

    % Preallocate 3D matrix
    mask = zeros(size(first_slice, 1), size(first_slice, 2), slices, 'double');

    % Read DICOM slices and convert to Hounsfield Units (HU)
    for i = 1:slices
        img = double(dicomread(fullfile(chdir, files(i).name)));
        % Convert to Hounsfield Units
        mask(:, :, i) = img * info.RescaleSlope + info.RescaleIntercept;
    end
end

function [ind] = pad_dim(dimension,mask,direction)
    if direction > 0
        if dimension == 3
            for i = 1:size(mask,dimension)
                if sum(mask(:,:,i),'all') == 0
                    continue
                end
                ind = i;
                break
            end
        elseif dimension == 2
            for i = 1:size(mask,dimension)
                if sum(mask(:,i,:),'all') == 0
                    continue
                end
                ind = i;
                break
            end 
        elseif dimension == 1
            for i = 1:size(mask,dimension)
                if sum(mask(i,:,:),'all') == 0
                    continue
                end
                ind = i;
                break
            end
        end
    else 
        if dimension == 3
            for i = 1:size(mask,dimension)
                if sum(mask(:,:,end-(i-1)),'all') == 0
                    continue
                end
                ind = i;
                break
            end
        elseif dimension == 2
            for i = 1:size(mask,dimension)
                if sum(mask(:,end-(i-1),:),'all') == 0
                    continue
                end
                ind = i;
                break
            end 
        elseif dimension == 1
            for i = 1:size(mask,dimension)
                if sum(mask(end-(i-1),:,:),'all') == 0
                    continue
                end
                ind = i;
                break
            end
        end
    end
end

function wedge_masks = generateWedgeMasks(mask)
    % Ensure binary mask
    mask = logical(mask);
    
    [rows, cols] = size(mask); % Get image dimensions
    cx = round(cols / 2); % X center
    cy = round(rows / 2); % Y center

    % Define angular segments (4 quadrants)
    num_wedges = 4;
    wedge_masks = false(rows, cols, num_wedges);
    
    % Find boundary of the mask
    boundary = bwperim(mask);
    [by, bx] = find(boundary);
    
    % Loop through quadrants
    for quadrant = 1:num_wedges
        % Create a blank mask for the current wedge
        wedge_mask = false(rows, cols);
        
        % Filter boundary points based on quadrant
        switch quadrant
            case 1 % Top-left
                valid = by < cy & bx < cx;
            case 2 % Top-right
                valid = by < cy & bx >= cx;
            case 3 % Bottom-left
                valid = by >= cy & bx < cx;
            case 4 % Bottom-right
                valid = by >= cy & bx >= cx;
        end
        bx_q = bx(valid);
        by_q = by(valid);
        
        % If boundary points exist in this quadrant
        if ~isempty(bx_q)
            % Fill wedge from center to boundary
            for k = 1:length(bx_q)
                wedge_mask = insertLine(wedge_mask, cx, cy, bx_q(k), by_q(k));
            end
            
            % Fill enclosed area
            wedge_mask = imfill(wedge_mask, 'holes');
        end
        
        % Store in output
        wedge_masks(:,:,quadrant) = wedge_mask;
    end
end

function img = insertLine(img, x1, y1, x2, y2)
    % Bresenham's line algorithm to draw a line on a binary image
    [x, y] = bresenham(x1, y1, x2, y2);
    for k = 1:length(x)
        if x(k) > 0 && y(k) > 0 && x(k) <= size(img,2) && y(k) <= size(img,1)
            img(y(k), x(k)) = 1;
        end
    end
end

function [x, y] = bresenham(x1, y1, x2, y2)
    % Calculate differences
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    
    % Determine step directions
    sx = sign(x2 - x1);
    sy = sign(y2 - y1);
    
    % Initialize starting points
    x = x1;
    y = y1;
    
    % Preallocate memory (worst-case scenario)
    n = max(dx, dy) + 1; 
    x_out = zeros(1, n);
    y_out = zeros(1, n);
    
    % Bresenham error term
    err = dx - dy;
    
    % Bresenham algorithm loop
    for i = 1:n
        x_out(i) = x;
        y_out(i) = y;
        
        if x == x2 && y == y2
            break;
        end
        
        e2 = 2 * err;
        if e2 > -dy
            err = err - dy;
            x = x + sx;
        end
        if e2 < dx
            err = err + dx;
            y = y + sy;
        end
    end
    
    % Trim outputs in case of preallocation overhead
    x = x_out(1:i);
    y = y_out(1:i);
end
