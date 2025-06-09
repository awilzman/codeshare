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
% Updated by Karen Troy 5/2025
% Example run:
% calibrate_slope = 0.00035619;
% calibrate_int = -0.00365584; 
%res in um, voxel edge length
% [bv,bmc,bmd,medial_left,angle_rot] = compare_dicoms(default_directory,res,LCV_name,
% mask1_name,mask2_name,calibrate_slope,calibrate_int,
% medial_left,angle_rot,first_full_slice)

% This version updated by Karen Troy to import the full grayscale mask of
% the segmented bone ([filename.AIM], to calulate BMD including marrow
% space.  Note that time 2 AIM files must be transformed prior to import.

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

    %printmasks option added to display slice 10  (line 204) - added by KLT
    printmasks=0;
    se = offsetstrel("ball",15,15); %shape for closing images later on

    %mask_LCV = get_mask(LCV_name,default_directory);
    mask_1 = get_mask(mask1_name,default_directory);
    mask_2 = get_mask(mask2_name,default_directory);
    
    %open full AIM files under a related variable name
    dir_path = strcat(chdir,'\',mask1_name,'_aim');
    mask_1fill=get_mask( '\', dir_path);
    dir_path = strcat(chdir,'\',mask2_name,'_aim');
    mask_2fill=get_mask( '\', dir_path);
    
    % un-Pad matrices
    mask_1 = pad_3dmat(mask_1);
    mask_2 = pad_3dmat(mask_2);

    %mask_2fill has an extra layer of zeros around all 4 edges compared to
    %mask_2. We need to remove them.
    mask_2fill = reduce_3dmat(mask_2fill);
    mask_2fill = padarray(mask_2fill,[0,0,1]); %add slice to top and bottom to match mask_2 dims

    tangle = pi/4;

    bone_threshold = 0.5/calibrate_slope;
    threshold_step = 0.01 / calibrate_slope;
    max_iter = 100;
    iter = 0;
    
        mask_1(mask_1 < 100) = 0;
        mask_2(mask_2 < 100) = 0;

    % fill in the holes of mask 1 and mask 2 and get the marrow gray values
    % from the aim files.
    bw=imclose(mask_1,se);
    bw(bw <100) = 0;
    bw(bw >=100) = 1;
    mask_1fill=bw.*mask_1fill;
    %figure; imshowpair(mask_1(:,:,40),mask_1fill(:,:,40),'montage');

    bw=imclose(mask_2,se);
    bw(bw <100) = 0;
    bw(bw >=100) = 1;
    mask_2fill=bw.*mask_2fill;
    %figure; imshowpair(mask_2(:,:,40),mask_2fill(:,:,40),'montage);
    clear bw;
    
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

        while true
            raw_image_b = raw_image > bone_threshold;
            raw_image_b = logical(raw_image_b);
        
            stats = regionprops(raw_image_b, 'Area', 'Centroid');
            stats = stats([stats.Area] > 5000);
        
            if numel(stats) == 2 || iter >= max_iter
                break;
            end
        
            bone_threshold = bone_threshold + threshold_step;
            iter = iter + 1;
        end
        disp(bone_threshold);
        figure; imshow(raw_image_b)

        if numel(stats) ~= 2
            error('Failed to isolate exactly two regions after threshold tuning.');
        end
        
        areas = [stats.Area];
        [~, sorted_idx] = sort(areas, 'descend');
        idx_tibia = sorted_idx(1);
        idx_fibula = sorted_idx(end);
        
        tibia_centroid = stats(idx_tibia).Centroid;
        fibula_centroid = stats(idx_fibula).Centroid;

        delta_y = fibula_centroid(1) - tibia_centroid(1);
        delta_x = fibula_centroid(2) - tibia_centroid(2);

        % Medial side calculation
        medial_left = tibia_centroid(2) >= fibula_centroid(2);

        % Angle calculation
        if medial_left
            angle_rot = atan2(delta_y, delta_x)-pi;
        else
            angle_rot = atan2(delta_y, delta_x)-tangle-pi;
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
    mask_1fill = rotate_mask(mask_1fill,angle_rot,ccenter);
    mask_2fill = rotate_mask(mask_2fill,angle_rot,ccenter);
    
    mask_1 = pad_3dmat(mask_1);  %crop to bone mask
    mask_2 = pad_3dmat(mask_2);
    mask_1fill = pad_3dmat(mask_1fill);
    mask_2fill = pad_3dmat(mask_2fill);

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
    
    % do the same for the filled mask. there's almost certainly a more
    % efficient way to do this, but I can't think of it now.
    mask_1f_med = zeros(size(mask_1fill));
    mask_1f_ant = zeros(size(mask_1fill));
    mask_1f_post = zeros(size(mask_1fill));
    mask_1f_lat = zeros(size(mask_1fill));
    
    mask_2f_med = zeros(size(mask_2fill));
    mask_2f_ant = zeros(size(mask_2fill));
    mask_2f_post = zeros(size(mask_2fill));
    mask_2f_lat = zeros(size(mask_2fill));
    
    % Iterate through slices

    med_s = 4;
    lat_s = 1;

    if medial_left
        med_s = 1;
        lat_s = 4;
    end

    %recalculate center after padding
    z_center = round(size(mask_1, 3) / 2); % Middle slice index
    center_slice = mask_1(:, :, z_center); % 2D slice along Z axis
    non_zero_voxels_2D = center_slice > 0;
    filled_slice = imfill(non_zero_voxels_2D, 'holes');
    [rows, cols] = find(filled_slice); 

    x = mean(cols);
    y = mean(rows);

    for z = 1:size(mask_1, 3)
        % Generate binary wedge masks for the current slice
        wedge_masks_1f= generateWedgeMasks(mask_1fill(:,:,z),x,y);
        wedge_masks_2f= generateWedgeMasks(mask_2fill(:,:,z),x,y);

        mask_1_slice = double(mask_1(:,:,z));
        mask_2_slice = double(mask_2(:,:,z));
    
        % Apply the binary wedge masks to extract the original values
        mask_1_med(:,:,z) = mask_1_slice.*double(wedge_masks_1f(:,:,med_s));
        mask_1_ant(:,:,z) = mask_1_slice.*double(wedge_masks_1f(:,:,2));
        mask_1_post(:,:,z) = mask_1_slice.*double(wedge_masks_1f(:,:,3));
        mask_1_lat(:,:,z) = mask_1_slice.*double(wedge_masks_1f(:,:,lat_s));

        mask_2_med(:,:,z) = mask_2_slice.*double(wedge_masks_2f(:,:,med_s));
        mask_2_ant(:,:,z) = mask_2_slice.*double(wedge_masks_2f(:,:,2));
        mask_2_post(:,:,z) = mask_2_slice.*double(wedge_masks_2f(:,:,3));
        mask_2_lat(:,:,z) = mask_2_slice.*double(wedge_masks_2f(:,:,lat_s));

        % Do the same for the filled masks so that we can calculate BMD
        mask_1_slice = double(mask_1fill(:,:,z));
        mask_2_slice = double(mask_2fill(:,:,z));
    
        % Apply the binary wedge masks to extract the original values
        mask_1f_med(:,:,z) = mask_1_slice.*double(wedge_masks_1f(:,:,med_s));
        mask_1f_ant(:,:,z) = mask_1_slice.*double(wedge_masks_1f(:,:,2));
        mask_1f_post(:,:,z) = mask_1_slice.*double(wedge_masks_1f(:,:,3));
        mask_1f_lat(:,:,z) = mask_1_slice.*double(wedge_masks_1f(:,:,lat_s));
    
        mask_2f_med(:,:,z) = mask_2_slice.*double(wedge_masks_2f(:,:,med_s));
        mask_2f_ant(:,:,z) = mask_2_slice.*double(wedge_masks_2f(:,:,2));
        mask_2f_post(:,:,z) = mask_2_slice.*double(wedge_masks_2f(:,:,3));
        mask_2f_lat(:,:,z) = mask_2_slice.*double(wedge_masks_2f(:,:,lat_s));
    end

    if printmasks==1
        slice_mask = mask_1_ant(:,:,10)+mask_1_post(:,:,10)+mask_1_med(:,:,10)+mask_1_lat(:,:,10);
        figure; imshow(slice_mask);
        imwrite(slice_mask,strcat(default_directory,mask1_name,"_maskslice10.png"));
        imwrite(mask_1(:,:,10),strcat(default_directory,mask1_name,"_slice10.png"));
    end

    [tv(1,1), bv(1,1), bmc(1,1), bmd(1,1)] = bv_bmc(mask_1_ant,mask_1f_ant,res,calibrate_slope,calibrate_int);
    [tv(1,2), bv(1,2), bmc(1,2), bmd(1,2)] = bv_bmc(mask_1_post,mask_1f_post,res,calibrate_slope,calibrate_int);
    [tv(1,3), bv(1,3), bmc(1,3), bmd(1,3)] = bv_bmc(mask_1_med,mask_1f_med,res,calibrate_slope,calibrate_int);
    [tv(1,4), bv(1,4), bmc(1,4), bmd(1,4)] = bv_bmc(mask_1_lat,mask_1f_lat,res,calibrate_slope,calibrate_int);
    [tv(2,1), bv(2,1), bmc(2,1), bmd(2,1)] = bv_bmc(mask_2_ant,mask_2f_ant,res,calibrate_slope,calibrate_int);
    [tv(2,2), bv(2,2), bmc(2,2), bmd(2,2)] = bv_bmc(mask_2_post,mask_2f_post,res,calibrate_slope,calibrate_int);
    [tv(2,3), bv(2,3), bmc(2,3), bmd(2,3)] = bv_bmc(mask_2_med,mask_2f_med,res,calibrate_slope,calibrate_int);
    [tv(2,4), bv(2,4), bmc(2,4), bmd(2,4)] = bv_bmc(mask_2_lat,mask_2f_lat,res,calibrate_slope,calibrate_int);
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

    slice_image = mask_1(:,:,round(size(mask_1,3)/2));
    [rows, cols] = size(slice_image);
    [xGrid, yGrid] = meshgrid(1:cols, 1:rows);
    ap_slope = tan(tangle);

    line_mask = abs(yGrid - (ap_slope * (xGrid - x) + y)) < 1;
    slice_image(line_mask) = max(slice_image(:));

    line_mask = abs(yGrid - ((1/(-ap_slope + 1e-6)) * (xGrid - x) + y)) < 1;
    slice_image(line_mask) = max(slice_image(:));

    imwrite(slice_image,strcat(default_directory,mask1_name,".png"))

    color_ant = [255, 0, 0];   % Red for anterior
    color_post = [0, 255, 0];  % Green for posterior
    color_med = [0, 0, 255];   % Blue for medial
    color_lat = [255, 255, 0]; % Yellow for lateral
    slice_image_rgb = zeros([size(slice_image), 3], 'uint8');
    
    for z = round(size(mask_1,3)/2) % 1:size(mask_1, 3)  -- this variable changed by KLT to only print the center slice for comparison
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
    se = strel('disk', 15);  % Adjust structuring element as needed
    closed_image = imclose(binary_image, se);
    closed_image(closed_image<20) = 0;

    % Calculate total filled area directly
    totalArea = sum(closed_image(:)); % Counts nonzero pixels (alternative to polyarea)
end

function [tv, bv, bmc, bmd] = bv_bmc(mask1,maskf, res, slope, int)
%mask1 = trabec segmented mask for bv, maskf = filled mask for other params    
    tv = 0;    
    bv = 0;
    bmc = 0;

    vox_ed = res / 10000.0; % Convert um to cm

    for z = 1:size(mask1, 3)
        slicef = double(maskf(:, :, z));
        slice = double(mask1(:,:,z));
        areaf = nnz(slicef);
        areab = nnz(slice);

        % Skip empty slices
        if areaf == 0
            continue;
        end

        bv_vol = areab * vox_ed^3;
        tv_vol = areaf * vox_ed^3;

        if bv_vol == 0
            continue;
        end
        %set non-bone values to -10000 so that we don't measure them
        slicef(slicef==0)=-10000;
        meanHU=mean(slicef(slicef > -9999), 'all')
        slice_density = mean(slicef(slicef > -9999), 'all') * slope + int;
        slice_content = slice_density * tv_vol;

        tv = tv + tv_vol; % Total Volume
        bv = bv + bv_vol; % Bone Volume
        bmc = bmc + slice_content;
    end

    if tv > 0
        bmd = bmc / tv;
    else
        bmd = 0;
    end
end

function new_mask = rotate_mask(mask, rotationAngle, centroid, interpolationMethod)
    assert(isnumeric(mask) && ndims(mask) >= 2, 'Input must be numeric 2D or 3D.');
    assert(isscalar(rotationAngle), 'Rotation angle must be scalar.');

    if nargin < 4
        interpolationMethod = 'linear';
    end

    % Compute padded size
    max_len = round(1.2 * hypot(size(mask,1), size(mask,2)));
    pad_TB = round((max_len - size(mask,1)) / 2);
    pad_LR = round((max_len - size(mask,2)) / 2);

    % Pad
    mask = padarray(mask, [pad_TB, pad_LR], 0, 'both');
    new_mask = zeros(size(mask));

    % Adjust centroid
    cy = centroid(1) + pad_TB;
    cx = centroid(2) + pad_LR;

    [rows, cols, slices] = size(mask);
    [X, Y] = meshgrid(1:cols, 1:rows);

    % Center relative to centroid
    x = X - cx;
    y = Y - cy;

    % Inverse rotation matrix
    x_r =  x * cos(rotationAngle) + y * sin(rotationAngle);
    y_r = -x * sin(rotationAngle) + y * cos(rotationAngle);

    % Map back to original coordinates
    xq = x_r + cx;
    yq = y_r + cy;

    for k = 1:slices
        new_mask(:,:,k) = interp2(double(mask(:,:,k)), xq, yq, interpolationMethod, 0);
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

function smaller_mask = reduce_3dmat(mask)
    % Trim a mask by 1 on each edge
    mask(isnan(mask)) = 0;
    [x, y, z] = ind2sub(size(mask), find(mask));
    
    % Compute bounds
    x1 = min(x); x2 = max(x);
    y1 = min(y); y2 = max(y);
    z1 = min(z); z2 = max(z);
    
    % Crop to bounding box
    smaller_mask = mask(x1+1:x2-1, y1+1:y2-1, z1:z2);
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

function wedge_masks = generateWedgeMasks(mask,cx,cy)
    % Ensure binary mask
    mask = logical(mask);
    
    [rows, cols] = size(mask); % Get image dimensions
    cx = round(cx); % X center
    cy = round(cy); % Y center

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
