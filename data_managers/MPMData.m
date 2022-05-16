%
% Class for storing and retrieving MPM data for the glass-patterns-three-session-mri 
% project. 
%
% ----------------------------------------------------------------------------------
%
% This requires the path to a .mat file (self.mpm_cluster) containing a structure with
% MPM data in the form:
%
%     results_mean.(smooth_type).(map_type).(mask).(tissue_type).[subj x scan table of MPM values]
%     e.g. results_mean.mni_gm_6mmsmoothed.MTsat.GM_DLPFC_R.GM
%
% e.g. see extract_map_values_master_clusters.m
%
% The class methods return data subtracted across sessions within subjects 
% (e.g. S1 scan 3 - S1 scan 2) for all subjects. There are keyword arguments for
% performing controls e.g. regress session 1. 
% 
% NOTE: there are two very similar functions here, one for 'ROI' and for 'cluster'. 
% For future use 'ROI' can be depreciated, this is a legacy of
% how the data structure was organised. Use the cluster method.
% 
% example useage:
%
% Initialise the class - 
%   mpm_data = MPMData(config). If not using mpm_roi you can delete self.mpm_roi fields and config from
%                               the constructure and just call as MPMData()
%              Make sure to fill in the correct path to the mpm data .mat file
%
% Get subtracted session data - 
%   mpm_data.get_subtract_mpm_cluster(["mni_gm_6mmsmoothed"], ...
%                                     ["MTsat"],  ...
%                                     ['all'],  ...
%                                     {'subject_01', 'subject_01'},  ...
%                                     'regress_s1',  ...
%                                     'Scan_3',  ...
%                                     'Scan_2');
%   see get_subtract_mpm_cluster() for more details on arguments.
%
%   This returns a structure with fields data (Subject X data) array and
%                                         labels data x 1 cell of labels containing info on data 
%
%   This output is compatable with the correlate_data() function of the class BivariateMethods.
% 
% ----------------------------------------------------------------------------------------------

classdef MPMData
    
    properties
        
        config = struct();
        mpm_roi = struct();
        mpm_cluster = struct();
        bvm = BivariateMethods;
        
    end

    
    methods
        
        function self = MPMData(config)
            
            self.config = config;  
            self.mpm_roi = load(self.config.results_mpm_roi);
            self.mpm_roi = self.mpm_roi.results_mean;
            self.mpm_cluster = load('D:\fMRIData\glass-patterns-three-session-mri\analysis\masks\analysis_masks\cluster_masks\final\mpm\cluster_extraction_values.mat'); 
            self.mpm_cluster = self.mpm_cluster.results_mean;

        end
        
        
        function [mpm_cluster_to_cor] = get_subtract_mpm_cluster(self, ...
                                                             analysis_fields, map_fields,  ...
                                                             cluster_fields,  ...
                                                             mask,  ...
                                                             control,  ...
                                                             varargin)
            %                                             
            % Return the subtracted mpm data for the field for the SPM clusters.
            %
            % see self.get_subtract_mpm_roi for most inputs. Only difference is 'cluster fields'
            % can be specified as an int containig the index of the field to return, or 'all' to return
            % all clusters
            %
            % ------------------------------------------------------------------------------------------
            
            [scan_to_subtract_from, scan_to_subtract] = handle_scan_to_subtract_varargin(varargin);
                        
            mpm_cluster_to_cor = struct('data', [NaN(length(mask), 1)],  ...  % cannot pre-allocate exact size as to difficult to detemrine
                                        'labels', {cell(1, 1)});              % number fo analyses that will be run

            cnt = 0;    
            for analysis = analysis_fields
                for map = map_fields
                        
                    % show specfic cluster fields or all 
                    all_cluster_fields = fields(self.mpm_cluster.(analysis).(map));

                    if isnumeric(cluster_fields)
                        cluster_fields = string(all_cluster_fields(cluster_fields))';
                    elseif strcmp(cluster_fields, 'all')
                        cluster_fields = string(all_cluster_fields)';
                    end

                    for cluster = cluster_fields        
                        label = strcat(analysis, '__', map, '__', cluster);
                        
                        subfield = string(fields(self.mpm_cluster.(analysis).(map).(cluster)));  
                        assert(length(subfield) == 1, 'too many subfields in MPM cluster');      

                        subtract_data = subtract_sessions(self.mpm_cluster.(analysis).(map).(cluster).(subfield)(mask, :),  ...
                                                          scan_to_subtract_from,  ...
                                                          scan_to_subtract);

                        s1 = self.mpm_cluster.(analysis).(map).(cluster).(subfield){mask, 'Scan_1'};
                        switch control
                            case 'divide_s1'
                                subtract_data = subtract_data ./ s1;
                            case 'regress_s1'
                                subtract_data = self.bvm.regress_data(s1, subtract_data, 'not_robust');
                            case 'no_control'
                        
                            otherwise
                                error('control not specified correcty');
                        end

                        cnt = cnt + 1;
                        mpm_cluster_to_cor.('data')(:, cnt) = subtract_data;
                        mpm_cluster_to_cor.('labels'){cnt} = label;

                    end
                end
            end
        end

        function [mpm_roi_to_cor] = get_subtract_mpm_roi(self, ... TODO: very similar to cluster - DRY
                                                         analysis_fields, map_fields,  ...
                                                         region_fields, seg_mask_fields,  ...
                                                         mask,  ...
                                                         control,  ...
                                                         varargin)
            %
            % Return specified MPM data subtracted across scans with/without controls applied
            %
            % INPUTS:
            %       MPM data structures are in the format:
            %       analysis_fields (e.g. mni/smoothed) > map_fields (e.g. MTsat, R2) > region_fields (e.g. LOC) > seg_mask_fields (e.g white_matter).
            %
            %       See self.mpm_roi for all options within the structre.
            %
            %       mask: subject_mask (e.g. {'subject_01', 'subject_02'}
            %       control: 'divide_s1' or 'regress_s1' to divide or regress normalise to Scan_1
            %       varargin: by default this will subtract Scan_3 - Scan_2. Otherwise input scans (see handle_scan_to_subtract_varargin)
            %
            % --------------------------------------------------------------------------------
            
            mpm_roi_to_cor = struct('data', [NaN(length(mask), 1)],  ...  % cannot pre-allocate exact size as to difficult to detemrine
                                    'labels', {cell(1, 1)});              % number fo analyses that will be run
                                
            [scan_to_subtract_from, scan_to_subtract] = handle_scan_to_subtract_varargin(varargin);

            if ismember('Scan_1', {scan_to_subtract, scan_to_subtract_from}) &&  ...
                ismember(control, {'clean_data', 'regress_s1', 'divide_s1'})
                error('Due to regression of S1, cleaning data is only valid for S3-S2 currently');                          
            end
            
            cnt = 0;
            for analysis = analysis_fields
                for map = map_fields
                    for region = region_fields
                        for seg_mask = seg_mask_fields
                            
                            if ~ismember(seg_mask,  ...
                                         fields(self.mpm_roi.(analysis).(map).(region)))
                                continue
                            end

                            label = strcat(analysis, '__', map, '__', region, '__', seg_mask);
                            subtract_data = subtract_sessions(self.mpm_roi.(analysis).(map).(region).(seg_mask)(mask, :),  ...
                                                              scan_to_subtract_from,  ...
                                                              scan_to_subtract);

                            % controls    
                            s1 = cell2mat(self.mpm_roi.(analysis).(map).(region).(seg_mask){mask, 'Scan_1'});
                            
                            switch control
                                case 'divide_s1'
                                    subtract_data = subtract_data ./ s1;
                                case 'regress_s1'
                                    subtract_data = self.bvm.regress_data(s1, subtract_data, 'not_robust');
                                case 'no_control'
                        
                                otherwise
                                    error(sprintf('control not specified correcty'));
                            end

                            cnt = cnt + 1;
                            mpm_roi_to_cor.('data')(:, cnt) = subtract_data;
                            mpm_roi_to_cor.('labels'){cnt} = label;
                            
                        end
                    end
                end
            end
        end
    
    end    
end