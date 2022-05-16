classdef fMRIMethods
    %
    % Class of useful methods for fMRI data analysis and preprocessing
    % ---------------------------------------------------------------------
        
    properties
       
        u = Utils;
        
    end
    
    
    methods
        
        function [average_nii] = make_average_nii(self, file_paths, varargin)
            % Make an average nii from the first volume of a set of 4D niis 
            %(or 3d niis) (with Nifti toolbox load_nii).
            %
            % file_paths - either a cell array of full paths to the niis to
            %              average, or a pathable where varargin is the base_dir
            %
            % returns the nii object made with make_nii (in case 4D input
            % header is different to 3D output).
            % --------------------------------------------------------------
            
            if istable(file_paths)
                base_dir = varargin{1};
                file_paths = self.u.make_paths('dont_add_suffix', 'cell', [], base_dir, file_paths);
            end

            num_scans = length(file_paths);
            average_nii = [];
            for irun = 1:num_scans

                fprintf('Running scan %d / %d\n', irun, num_scans);
                nii = load_nii(...
                               char(file_paths{irun}));
                if irun == 1

                    average_nii = nii;
                    average_nii.img = nii.img(:, :, :, 1);
                    continue

                else
                    average_nii.img = average_nii.img + nii.img(:, :, :, 1); 
                end

            end
            average_nii.img = average_nii.img / num_scans;

            if size(nii.img, 4) ~= size(average_nii.img, 4)
                % give average nii header the new 4th dim size (if going from 4D to 3D) 
                average_nii.hdr.dime.dim(5) = size(average_nii.img, 4);
            end

        end



        function [xY_nii] = convert_SPM_xY_to_nii(~, template_nii, xY)
            %
            % www.slicer.org/wiki/Coordinate_systems
            % works for any timeseries output (e.g. can be brainmasked)
            % -------------------------------------------------------------

            template_dims = size(template_nii.img);
            vol_dims = template_dims(1:3);
            num_vols = template_dims(4);

            % convert real-world to matrix dims
            rot_mat = xY.spec.mat;
            pad_realworld_idx = [xY.XYZmm; 
                                 ones(1, size(xY.XYZmm, 2))];

            vox_idx = inv(rot_mat) * pad_realworld_idx;
            idx = sub2ind(vol_dims,  ...
                          vox_idx(1, :),  ...  x
                          vox_idx(2, :),  ...  y
                          vox_idx(3, :));   %  z

            % make a new volume filled with SPM xY.y
            output_4D = zeros([vol_dims,  ...
                               num_vols]);
            for ivol = 1:num_vols

                vol_3D = zeros(vol_dims);
                vol_3D(idx) = xY.y(ivol, :);
                output_4D(:, :, :, ivol) = vol_3D;

            end

            template_nii.img = output_4D;
            xY_nii = template_nii;

        end


        function [eigenvariate] = get_spm_eigenvariate(~, y)
            %
            % Taken directly from SPM's spm_regions fuction, written by
            % Karl Friston / Guillaume Flandin
            % y can be the first output of mask_nii_timeseries (transposed)
            % -------------------------------------------------------------

            [m,n]   = size(y);
            if m > n
                [v,s,v] = svd(y'*y);
                s       = diag(s);
                v       = v(:,1);
                u       = y*v/sqrt(s(1));
            else
                [u,s,u] = svd(y*y');
                s       = diag(s);
                u       = u(:,1);
                v       = y'*u/sqrt(s(1));
            end
            d       = sign(sum(v));
            u       = u*d;
            v       = v*d;
            eigenvariate = u*sqrt(s(1)/n);   

        end
        
        
        function [voxel_tSNR, voxel_mean, voxel_std] = calculate_tsnr(~, voxel_timeseries, varargin)
            % Calculates the tSNR per voxel for a voxel X time matrix.
            %
            % INPUTS: voxel_timeseries: voxel X volume array of voxel timeseries
            %
            %         varargin: Save volume of voxel tSNR - requires 3 inputs -
            %                   varargin{1} = full path (including filename) to saved for saved volume.
            %                   varargin{2} = binary 3D mask
            %
            % ---------------------------------------------------------------------------------------------------

            % calculate tSNR for each individual voxel
            voxel_mean = mean(voxel_timeseries, 2);
            voxel_std = std(double(voxel_timeseries), 0 , 2);
            voxel_tSNR = voxel_mean ./ voxel_std;

            % write to tSNR nifti
            if ~isempty(varargin)
                full_filepath_to_save = varargin{1};
                mask_voxels = varargin{2};
                mask_voxels = logical(mask_voxels);

                tSNR_volume = zeros(size(mask_voxels));
                tSNR_volume(mask_voxels) = voxel_tSNR;
                save_tSNR_nii = make_nii(tSNR_volume);
                save_nii(save_tSNR_nii, full_filepath_to_save);
            end

        end

        function [slice_range] = range_of_bold_across_slices(~, nifti, mask_voxels)
            % Find the range of BOLD values across the volume. Typically the
            % top slices will have large BOLD values but they should be fairly
            % consistent throughout the volume. If there are large spikes on
            % certain slices in the middle something has gone wrong.
            % ---------------------------------------------------------------------

            for slice = 1:size(nifti,3)
                for time = 1:size(nifti, 4)
                    slice_voxels = nifti(:, :, slice, time);
                    slice_timeseries(slice, time) = mean(slice_voxels(mask_voxels(:, :, slice)));
                end
            end
            slice_range = range(slice_timeseries, 2); 

        end


        function [coreg_correlation] = check_align(~, nii1, nii2, varargin)
            % perform a voxelwise correlation of two images.
            % Despite this very basic implimentation it
            % seems to perform similarly to the SPM function for MPMs and
            % fMRI and so can be used as a quick convenient check. 
            % ----------------------------------------------------------------------

            if ~isempty(varargin)
                mask_voxels = varargin{1};
                mask_voxels(isnan(mask_voxels)) = 0;
                nii1 = nii1 .* mask_voxels;
                nii2 = nii2 .* mask_voxels;
            end

            nii1 = reshape(nii1, [numel(nii1), 1]);
            nii2 = reshape(nii2, [numel(nii2), 1]);

            nii1 = double(nii1);
            nii2 = double(nii2);

            coreg_correlation = corr(nii1, nii2);

        end


        function [union_mask] = make_union_mask(~, masks)  
            %
            % Make a union mask nii structure from a cellstr of filenames
            % i.e. {{mask1.nii},{mask2.nii}, ...}
            %
            % ------------------------------------------
            
            first_mask = load_nii(fullfile(masks{1}));
            union_mask = first_mask;

            for irun = 2:length(masks)
                
                other_mask = load_nii(fullfile(masks{irun}));
                union_mask.img = union_mask.img | other_mask.img;
                
            end
            
        end
                
                   
        function [intersection_mask] = make_intersection_mask(self, mask1, mask2, thr1, thr2)  
            % Overlap two masks. If is a structure, assume load_nii input, otherwise a path.
            %
            % Output load_nii structure of interseection between all specified masks. 
            % 
            % Can handle binary masks or probability maps, if probability maps supply thresholds, otherwise []
            % --------------------------------------------------------

            if ~isstruct(mask1) && ~isstruct(mask2)
                mask1_nii = load_nii(mask1);
                mask2_nii = load_nii(mask2);

            elseif isstruct(mask1) && isstruct(mask2)
                mask1_nii = mask1;
                mask2_nii = mask2;

            else
                error('input types must match');
            end

            if ~isempty(thr1)
                mask1_nii = self.threshold_probability_map(mask1_nii, thr1);
            end
            
            if ~isempty(thr2)
                mask2_nii = self.threshold_probability_map(mask2_nii, thr2);
            end
            
            
            mask1_nii.img(isnan(mask1_nii.img)) = 0;  
            mask2_nii.img(isnan(mask2_nii.img)) = 0;  

            intersection_mask = mask1_nii;
            intersection_mask.img = mask1_nii.img & mask2_nii.img; 

        end


        function [thresholded_mask] = threshold_probability_map(~, map, thr) 
            %
            % Take a loaded probability map and threshold it into a binary mask.
            %
            % returns nii structure 
            % -------------------------------------------------------------
            
            map.img(map.img<thr) = 0;
            map.img(map.img>=thr) = 1;
                        
            thresholded_mask = map;

        end


        function [overlap] = get_num_voxel_overlap(self, nii1, nii2, thr1, thr2)
            % 
            % calculate number of voxels in overlap between two masks (dir input)
            % that are PMs (threshold) (use any threshold btween 0-1 if is binary mask
            %
            % seems to perfectly match for this use mask=false for proper interpolation
            % when resizing for proper interpolation (also it is not mask it is PM)
            % -------------------------------------------------------------------------
            if ~isempty(thr1)
                nii1 = self.threshold_probability_map(nii1, thr1);  % TODO: drive with above, own function
            end
            
            if ~isempty(thr2)
                nii2 = self.threshold_probability_map(nii2, thr2); 
            end
            
            intersection = logical(nii1.img) & logical(nii2.img);
            overlap = nnz(intersection);

        end


        function check_voxel_number(~, full_filepath)
            %
            % Convenience function to check and print the number of voxels in a binary mask
            %
            % -----------------------------------------------------------

            nii = load_nii(full_filepath);
            num_voxels = nnz(nii.img);
            
            [~, filename, ext] = fileparts(full_filepath);
            filename = [filename, ext];
            fprintf('%s NUM VOXELS: %d\n', filename, num_voxels);
                      
        end
        
        function resize_mutiple_niis(~, files, vox, bb, ismask)
            %
            % Use Ged Ridgeway's resize_img on multiple files. (see for arguments)
            % All files must be a mask or not a mask (ismask bool).
            % Output is r-prefixed niftis in same folder as original
            %
            % INPUT: currently the output of dir(). Update for pathtable. 
            % -----------------------------------------------

            for irun = 1:length(files)

                nii_path = fullfile(files(irun).folder, files(irun).name);
                resize_img(nii_path, vox, bb, ismask);

            end

        end


        function [dilated_nii] = dilate_nii(~, nii, voxels_to_add, smooth_)
            % Dilate a nii (used for basic brainmasking). 
            % Works poorly as strange speckling at edge but sufficient
            % for correcting segmentation error, can smooth to remove. 
            % ----------------------------------------------------------------------
            non_zero_idx = find(nii.img > 0);
            non_zero_idx_add = non_zero_idx + voxels_to_add;
            non_zerm_idx_minus = non_zero_idx - voxels_to_add;

            dilated_nii_idx = union(non_zero_idx_add,  ...
                                    non_zerm_idx_minus);

            dilated_nii = zeros(size(nii.img));
            dilated_nii(dilated_nii_idx) = 1;

            if smooth_
                dilated_nii = smooth3(dilated_nii);
                dilated_nii(dilated_nii > 0) = 1;
            end

            nii.img = dilated_nii;
            dilated_nii = nii;

        end


        function [centre_x, centre_y, centre_z] = convert_mrs_voxel_centre_point_to_mni_coords(~, file_path) 
            % Convert Centre of voxel coordinates to MNI
            % from https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;dae1f50e.0903
            % https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html
            %
            % Use median so the output is an actual voxel (i.e. multiple of the voxel dimension).
            % Can alternatively use mean, makes little difference, mean is probably better.
            % TODO: change when moving onto next project.
            %
            % An alternative (sanity check but actually nicer) implementation is:
            % vol = spm_vol(file_path);
            % [y, xyz] = spm_read_vols(x);
            %
            % x = median(xyz(1, find(y > 0));  where file_path is path to binary mask.
            % y = median(xyz(2, find(y > 0));
            % z = median(xyz(3, find(y > 0));
            % ----------------------------------------------------------------------------------------

            % load nii, find XYZ index of all mask voxels
            V = spm_vol(file_path);
            Vol = spm_read_vols(V);
            all_mask_voxels = find(logical(Vol)==1);
            [x, y, z] = ind2sub(size(Vol), all_mask_voxels);

            % find the centre of the voxel as half of each edge
            centre = [median(x), median(y), median(z), 1];
            assert(logical(Vol(centre(1), centre(2), centre(3))) == 1, 'centre voxel doesnt exist');

            % conver xyz to mni coordinates
            mni_coords = V.mat * centre';
            centre_x = mni_coords(1);
            centre_y = mni_coords(2);
            centre_z = mni_coords(3);

        end
 
        
        function [mask] = make_sphere_mask(~, ref_img, mni_center, radius)
            % Make a speherical mask from:
            %
            % ref_img - path to a reference image (used for dimensions),
            %           must be 3D with isotropic voxel size.
            % 
            % mni_center - [x,y, ] array of sphere center in mni
            %              coordinates
            %
            % radius - radius of sphere in mm
            %
            % -------------------------------------------------------------
            
            V = spm_vol(ref_img);
            warning('ensure voxels are isotropic (this warning is always shown');

            % convert mni centre and sizes from mni coordinate to idx 
            center_idx = inv(V.mat) * [mni_center 1]';
            r = int16(radius / abs(V.mat(1, 1)));  
         
            mask = load_nii(ref_img);
            mask.img(:, :, :, :) = 0;
            
            for x = -r:r
                for y = -r:r
                    for z = -r:r
                        
                        x_idx = int16(center_idx(1)) - x;
                        y_idx = int16(center_idx(2)) - y;
                        z_idx = int16(center_idx(3)) - z;
                        
                        within_threshold = r.^2 > abs(x).^2 + abs(y).^2 + abs(z).^2;
                        
                        if within_threshold
                            mask.img(x_idx, y_idx, z_idx) = 1;
                        end
                                          
                    end
                end
            end         
            
        end
        
        
        function [mask] = make_rectangle_mask(~, ref_img, mni_centre, x, y, z)
            %
            % See make sphere mask, however now supply length in mm for x,
            % y, z.
            %
            % -----------------------------------------------------------------
            
            V = spm_vol(ref_img);
            
            % convert mni centre and sizes from mni coordinate to idx 
            idx = inv(V.mat) * [mni_centre 1]';
            x = int16(x / V.mat(1, 1));
            y = int16(y / V.mat(2, 2));
            z = int16(z / V.mat(3, 3));
            
            % half x to move so it can be padded either side of centre
            x_size = int16(x / 2);
            y_size = int16(y / 2);
            z_size = int16(z / 2);
            
            % round centre idx and mask
            x_centre = int16(idx(1));
            y_centre = int16(idx(2));
            z_centre = int16(idx(3));
            
            mask = load_nii(ref_img);
            mask.img(:, :, :, :) = 0;
           
            mask.img(x_centre - x_size: x_centre + x_size,  ...
                    y_centre - y_size: y_centre + y_size,  ...
                    z_centre -  z_size: z_centre + z_size) = 1;

        end
        

        function [mask] = make_rough_mask(~, nii_to_mask, mask_opt, mask_setting, save_mask_opt)
            % Make a rough mask for fast checking of niftis. This is
            % not a substitue for proper masking in SPM. Only use to get a rough mask
            % for quick data checks during data collection
            %  
            % INPUTS: nii_to_mask:  must be 4D array from niftiread or load_nii (.img)
            %
            %         mask_opt / mask_settings: 
            %                     'cube': make a cube ROI with coordinates provided in
            %                           varargin. coordinates must be provided in a structure with mask.x,
            %                           mask.y, mask.z each containing an array of voxels to include in the mask
            %                           along that axis i.e. mask.x = 50:100
            %
            %                    'pseudobrain': takes the biggest contigious area in a
            %                                   binarised copy of the nii, with threshold set by
            %                                   varargin. This threshold is necessary because dark
            %                                   areas in an EPI are typically not 0. Look at nifti to
            %                                   determine a threshold to remove dark areas, 2000 is
            %                                   typical. 
            %                    'threshold': mask = 1 for anything above the given threshold
            %                                 (in BOLD units)
            %
            %         save_mask_opt: 'dont_save' for dont save, or full path + filename to save to.
            %         
            %  -------------------------------------------------------------------------------------------------------

            mask_voxels = zeros(size(nii_to_mask(:,:,:,1)));  

            switch mask_opt

                case 'cube'
                     alt_mask = mask_setting;
                     mask_voxels(alt_mask.x, alt_mask.y, alt_mask.z, 1) = 1;

                case 'pseudobrain'
                    B_mask_thr = mask_setting;
                    mask_voxels(nii_to_mask(:,:,:,1) > WB_mask_thr) = 1;

                    conncomp  = bwconncomp(mask_voxels, 6);  % find biggest contionus part of mask
                    [~, maxcell] = max(cellfun(@numel, conncomp.PixelIdxList));
                    mask_voxels = zeros(size(mask_voxels));      % Zero the image and assign to it the largest component.
                    mask_voxels(conncomp.PixelIdxList{1, maxcell}) = 1;

                case 'threshold'
                     threshold = mask_setting;
                     mask_voxels(nii_to_mask(:,:,:,1) > threshold) = 1;

                otherwise
                   error("mask option must be 'cube', 'pseudobrain', 'threshold'")
            end

            mask = logical(mask_voxels);

            % handle save option
            if ~strcmp(save_mask_opt, 'dont_save')
                save_mask_fullfilepath = save_mask_opt;
                mask_ROI = nii_to_mask(:,:,:,1);
                mask_ROI(mask) = 0;
                save_mask_nii = make_nii(mask_ROI);
                save_nii(save_mask_nii, save_mask_opt);
            end    

        end

        function [R] = get_motion_derivatives_squares_demean(~, R)
            % Take the derivative, squares and squared derivatives for
            % inclusion in GLM or use in denoising / for regressing out ICA
            % components before denoising. R - column vector of motion
            % -------------------------------------------------------------

            reg_deriv = [zeros(1, size(R,2)); 
                         diff(R, 1, 1)];

            reg_square = [R.^2 reg_deriv.^2];
            R = [R reg_deriv reg_square];

            % demean motion parameters. 
            for i = 1:size(R,2)
                R(:,i) = R(:,i) - mean(R(:,i));
            end

        end

        
        function [scans] = get_4D_nii_path_for_spm(~, nii_path)
            % If a nii is 4D, the path cannot be put in as-is to SPM. It
            % needs to be formatted with the vol number appended to the
            % path.
            % ---------------------------------------------------------------

            vols = load_nii_hdr(nii_path);
            vols = vols.dime.dim(5);    % read number of volumes from Nifti header
            for i = 1:vols
                scans{i,1} = strcat(nii_path, ',', num2str(i));
            end

        end
    
        function [filtered_mask_eigenvariate, filtered_masked_voxels,  ...
                    mask_eigenvariate, masked_voxels] = mask_eigenvariate_filter_data(self, nii, mask_nii, TR, low, high, filter_type, varargin)
            %
            % Main function for masking the brain, calculating eigenvariate
            % of the masked region and filtering it.
            %
            % INPUTS
            %
            % nii - 
            %
            %
            %
            
            % Voxel X Time input
            % Convenience function to extract the masked, bandpass-filetered eigenvariate
            % Have to transpose and re-transpose to use bp_filter 
            % -------------------------------------------------------------------------------------------

            if ~isequal(mask_nii, false)
                [masked_voxels] = self.mask_nii_timeseries(nii, mask_nii);
            else
                warning('no mask provided');
                [masked_voxels] = self.nii_to_vox_by_time(nii.img);
            end
       
            mask_eigenvariate = self.get_spm_eigenvariate(masked_voxels');
            filtered_mask_eigenvariate = self.bp_butter_filter(mask_eigenvariate, TR, low, high, filter_type);
            
            
            % Option for filtered vox by time matrix of voxel timeseries.
            % This cannot be done in bp_butter_filter as requires some
            % awkward transposes to fit into matlab function requirements.
            % Not run by default as takes quite a long time.
            if ~isempty(varargin) && varargin{1}
                
               [b, a] = self.make_bp_butter_filter(TR, low, high, filter_type);
               filtered_masked_voxels = filtfilt(b, a,  ... 
                                                 double(masked_voxels'));
               filtered_masked_voxels = filtered_masked_voxels';
               
            else
                filtered_masked_voxels = false;
            end
            
  
        end
        
        function [b, a] = make_bp_butter_filter(self, TR, low, high, filter_type)
           %
           % Convenience function to get filter coefficients for highpass
           % or bandpass
           %---------------------------------------------------------------------
           if strcmp(filter_type, 'highpass')
               [b, a] = self.make_bp_butter_filter_highpass(TR, low);         
               
           elseif strcmp(filter_type, 'bandpass')
               [b, a] = self.make_bp_butter_filter_bandpass(TR, low, high);
           end
          
            
        end
        
        
        function [filtered_data] = bp_butter_filter(self, data, TR, low, high, filter_type)
           %
           % zero-phase filter data with butterworth filter. 
           % Input (data) must be a single timeseries.
           %
           % INPUTS:
           %     data - 1 dimensional timeseries to filter
           %    
           %     TR - repetition time in seconds 
           % 
           %     low, high - filter cutoff in Hz (if highpass, low is
           %     ignored)
           %
           %     filter_type - 'bandpass' or 'highpass'
           %
           % NOTE: filter is 5th order but filtfilt passes throguh twice
           % for zero-phasing so effective order is 10. 
           % --------------------------------------------------------------

           [b, a] = self.make_bp_butter_filter(TR, low, high, filter_type);
          
           filtered_data = filtfilt(b, a,  ...
                                    double(data));  % filtfilt input must be double
                                
        end
        
        
        function [voxel_timeseries, voxel_timeseries_mean] = mask_nii_timeseries(~, nifti, mask_voxels)
            % extracts the timeseries per voxel of an 4D nifti (loaded with load_nii) within a binary 
            % nifti mask. Outputs to a voxel X time matrix and returns mean along time
            % dimension. 
            %
            % INPUTS: 
            %       1) a nifti loaded using niftiread (4d array), or load_nii (structure)            )
            %       2) a 3D binary mask in same space/dimensions as fMRI data (loaded with load_nii or nifti_read)      
            %
            %       NOTE: if input is structure will assume load_nii, else if nifti is
            %       array will assume niftiread load.
            %
            % OUTPUT: timeseries of all voxels within mask, average timeseries. voxel X volume
            % ---------------------------------------------------------------------------------------------------------------

            % Handle inputs        
            if isstruct(nifti)
                nifti = nifti.img;
            end
            if isstruct(mask_voxels)
                mask_voxels = mask_voxels.img;
            end
            mask_voxels(isnan(mask_voxels)) = 0;
            mask_voxels = logical(mask_voxels);        

            % mask out timeseries and save
            n_voxels_in_mask = sum(mask_voxels, 'all');
            n_time_steps = size(nifti, 4);
            voxel_timeseries = nan(n_voxels_in_mask, n_time_steps);
            
            for itime = 1:n_time_steps

                volume_voxels = nifti(:, :, :, itime);
                voxel_timeseries(:, itime) = volume_voxels(mask_voxels);

            end

            voxel_timeseries_mean = mean(voxel_timeseries); 

        end
        
        
        function [b, a] = make_bp_butter_filter_highpass(~, TR, low)
            
           fs = 1/TR;
           wn_max = low/(fs/2);
           
           [b, a] = butter(5, wn_max, 'high');
            
        end
        
        
        function [b, a] = make_bp_butter_filter_bandpass(~, TR, low, high)

           fs = 1/TR;
           wn_min = low/(fs/2);
           wn_max = high/(fs/2);

           [b, a] = butter(5, [wn_min, wn_max], 'bandpass');

        end
              
        function [coreg_cost_function] = check_coregistration_spm(self, nii1_path, nii2_path, mask_voxels, cost_fun)
            %
            % Quick convenience function for sanity checking alignment
            % -------------------------------------------------------------
            if isstruct(mask_voxels)
                mask_voxels = mask_voxels.img;
            end
            mask_voxels(isnan(mask_voxels)) = 0; 
            mask_voxels = logical(mask_voxels);  

            V1 = self.load_and_mask_spm_vol(nii1_path, mask_voxels);
            V2 = self.load_and_mask_spm_vol(nii2_path, mask_voxels);
            
            [coreg_cost_function] = jz_spm_checkalign(V1, V2, cost_fun);  % ncc of joint histogram

        end


        function [V] = load_and_mask_spm_vol(~, nii_path, mask_voxels)
            %
            %
            %--------------------------------------------------------------
            
            V = spm_vol(nii_path);
            if iscell(V)
                V = V{:};
            end

            V.uint8 = jz_orig_spm_loaduint8(V);

            if mask_voxels
                V.uint8 = V.uint8 .* uint8(mask_voxels);
            end

        end

        
        function [mask_niis] = load_all_masks_into_struct(~, base_dir)
            %
            % Load all masks into a structure (using load_nii). Useful when
            % using the same mask to mask multiple subjects / sessions as
            % avoids re-loading. 
            % 
            % INPUT: a directory containing masks to load. All nifti files 
            %        in the dir are assumed to be masks and will be loaded.
            %
            % OUTPUT:  mask_niis.(field_name) = load_nii() loaded mask
            %
            % -------------------------------------------------------------
                                                                                                     
            mask_paths = dir(fullfile(base_dir, '*.nii'));
            mask_names = extractBefore({mask_paths.name}, '.nii');

            mask_niis = struct();    
            for irun = 1:length(mask_names)

                name = mask_names{irun};
                mask_niis.(name) = load_nii(fullfile(mask_paths(irun).folder,  ...
                                                     mask_paths(irun).name));
            end

            
        end
        

        function [nii_4D] = quick_convert_3D_to_4D(self, base_dir, pathtable)
            % Quick script to take a pathtable of 3D nifti and load into a single 
            % 4D array. NOTE: not a substitute for SPM 3D to 4D conversion !
            % ----------------------------------------------------------------------

            first_nii = load_nii(self.u.make_paths('no_suffix', 'char', [], base_dir, pathtable(1, :))); 
            dims = size(first_nii.img);

            nii_4D = nan(dims(1), dims(2), dims(3), height(pathtable));
            nii_4D(:, :, :, 1) = first_nii.img;   
            for ivol = 2:height(pathtable)

                vol =  load_nii(self.u.make_paths('no_suffix', 'char', [], base_dir, pathtable(ivol, :)));   % most time here is spend loading file
                nii_4D(:, :, :, ivol) = vol.img;

            end
            
            
        end
   
        function [reshaped_nii, dims] = nii_to_vox_by_time(~, nii_4D)
            % convert a [x y z t] nifti (in 4D matrix form)
            % to 2D voxel x time matrix. 
            % ---------------------------------------------
                [x, y, z, num_vols] = size(nii_4D);
                num_voxels = x * y * z;
                reshaped_nii = reshape(nii_4D, num_voxels, num_vols);  
                dims = [x y z];
        end
        
        
        function [nii] = vox_by_time_to_nii(~, vox_by_t, dims)
            % Convert a 3D voxel x time matrix back into 
            % a 4D nii.
            %
            % vox_by_t : a 2D voxel by time matrix
            % dims : dimensions (3D) of the nii in the form [x y z]
            % ---------------------------------------------
            t = size(vox_by_t, 2);
            
            nii = reshape(vox_by_t,  ...
                          dims(1),  ...
                          dims(2), ...
                          dims(3),  ...
                          t);
        
        end
        
        function [summary_stat_ALFF, summary_stat_fALFF, ...
                    mean_voxelwise_ALFF, mean_voxelwise_fALFF] = get_alff_measures(self, summary_stat, voxels, low, high, fs)
            %
            % Voxels matrix must be t x voxel
            % return the alf from a summary statistic (i.e. eigenvariate and the voxels timeseries
            % only reason these are both input together is more convenience. 
            % ------------------------------------------------------------------------------------------
            assert(size(summary_stat, 1) == size(voxels, 1), 'inputs must be t x vox');
         
            [summary_stat_ALFF, summary_stat_fALFF] = self.compute_alff(summary_stat, low, high, fs);
            [voxelwise_ALFF, voxelwise_fALFF] = self.compute_alff(voxels, low, high, fs);
            
            mean_voxelwise_ALFF = mean(voxelwise_ALFF);
            mean_voxelwise_fALFF = mean(voxelwise_fALFF);                
                
        end
        
        
        function [ALFF, fALFF] = compute_alff(~, tc, f_low, f_high, fs)
            % 
            % For voxelwise, the matrix must be t x voxel
            %
            % coded according to the REST v1.22 toolbox    
            % adapted from vmk25@cam.ac.uk script  'extract_roi_alff.m' 
            % ---------------------------------------------------------
            N = size(tc,1);
            N_padded = 2^nextpow2(N);
            idx_low = ceil(f_low * N_padded / fs + 1);      % because first entry corresponds to DC signal
            idx_high = fix(f_high * N_padded / fs + 1);
            tc = [tc; zeros(N_padded - N, size(tc,2))];	% padded with zeros

            % run ALFF calculation on input image
            new_img = 2*abs(fft(tc))/N;
            ALFF = single(mean(new_img(idx_low:idx_high,:)));
            fALFF = single(sum(new_img(idx_low:idx_high,:)) ./ sum(new_img(2:(N_padded/2 + 1),:)));
            fALFF(~isfinite(fALFF)) = 0;
            
        end
 
    end
    
end    
    
    
    
    
    