% Class for the analysis of MRS Data including methods to get data subtracted 
% across scans with / without controls. See glass-patterns-three-session-mri
% repo for example use.
% 
% Requires:
%     self.mrs = a struct containing MRS data where each key describes the data type and
%                the data is a subject x session table of MRS data.
%                 
%                The struct keys should be in the format mrstype_region_normalisation e.g. 'gaba_loc_naa' or 'glu_ips_water'
%
%                In the glass-patterns project this is generated by the script mrs_data_msgp() from .csv data
%
%     self.bvm = the BivariateMethods class (see glass-patterns project)
%
%     self.config = a configuration structure, here the only required fields is:
%                   config.self.config.mrs_csv_path and self.config.mrs_voxel_tissue_proportion
%                   which are paths to a .csv containig GABA data and a .mat file containing tissue proportion data.
%   
%
% -------------------------------------------------------------------------------------------------------------------------


classdef MRSData
    
    properties
        
        config = struct();
        mrs = struct();
        tissue_proportion = struct();
        bvm = BivariateMethods;   
        
    end

    
    methods
        
        function self = MRSData(config)
            
            self.config = config;                                      
            self.mrs = mrs_data_msgp();
            self.tissue_proportion = load(config.mrs_voxel_tissue_proportion);  
            self.tissue_proportion = self.tissue_proportion.tissue_proportion;
            
        end
    
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%       Get subtracted data processed after / not after controls 
% --------------------------------------------------------------------------------------------------------------------------------------------------------------

        function [mrs_to_cor] = get_subtract_mrs(self, fields_to_get, mask, control, varargin)
            % 
            % Get MRS data subtracted between sessions and with / without controls applied. 
            %
            % INPUTS: 
            %       fields_to_get: str array of keys key for self.mrs struct e.g. ["gaba_loc_naa", "glu_loc_naa"].
            %
            %       mask: cell str of subject IDS to include in the analysis. The rownames for the self.mrs subject x data
            %             tables are used for this.
            %
            %       control: type of control to apply (see swtich statement below for possible inputs)
            %
            %       varargin: by default Scan_2 is subtracted from Scan_3, but input values to specifiy. The VariableNames
            %                 from the self.mrs subject x data tables are used as Scan names.
            %
            % OUTPUT:
            %       an output structure with fields data - a subject x num-data-requested table of subtracted values
            %                                       labels - labels describing the data for future reference. This
            %                                       output is compatible with BiviarateMethods class correlate_data() function
            %
            % EXAMMPLE:
            %       mrs_data = MRSData(config);
            %       mrs_to_cor = mrs_data.get_subtract_mrs(["gaba_loc_water", "glu_loc_naa"],  {'subject_01', 'subject_02', subject_03'}, ...
            %                                             'regress_csf', 'Scan_2', 'Scan_1');
            %
            %       This will return a structure where data is a 3 x 2 array with rows as subject and columns are session 2 - session 1
            %       data with csf regressed, column 1 is GABA LOC over Water and column 2 is GLU LOC over NAA.
            % -----------------------------------------------------------------------------------------
            
            mrs_to_cor = struct('data', [NaN(length(mask),  length(fields_to_get))],  ...
                               'labels', {cell(1, length(fields_to_get))}); 
                            
            [scan_to_subtract_from, scan_to_subtract] = handle_scan_to_subtract_varargin(varargin);

            if ismember('Scan_1', {scan_to_subtract, scan_to_subtract_from}) &&  ...
                ismember(control, {'clean_data', 'regress_s1', 'divide_s1'})
                error('Due to regression of S1, cleaning data is only valid for S3-S2 currently');                          
            end
            
            cnt = 0;
            for analysis = fields_to_get

                subtract_data = subtract_sessions(self.mrs.(analysis)(mask, :),  ...
                                                  scan_to_subtract_from, scan_to_subtract);
                                              
                % handle controls - only 1 control may be entered per call
                s1 = self.mrs.(analysis){mask, 'Scan_1'};
                switch control
                    case 'divide_s1'
                         subtract_data = subtract_data ./ s1;
                         
                    case 'regress_s1'
                         subtract_data = self.bvm.regress_data(s1, subtract_data, 'not_robust');
                    
                    case 'regress_glu_loc'
                        glu = subtract_sessions(self.mrs.('glu_loc_naa')(mask, :),  ...
                                                scan_to_subtract_from,  ...
                                                scan_to_subtract);
                        subtract_data = self.bvm.regress_data(glu, subtract_data, 'not_robust'); 
                    
                    case 'regress_linewidth'
                        subtract_data = self.regress_out_linewidth(analysis, mask, scan_to_subtract_from, scan_to_subtract);
                        
                    case 'regress_csf' 
                        subtract_data = self.regress_or_divide_csf_mrs(analysis, mask, 'regress',  ...
                                                                       scan_to_subtract_from, scan_to_subtract); 
                    case 'divide_csf'
                        subtract_data = self.regress_or_divide_csf_mrs(analysis, mask, 'divide',  ...
                                                                       scan_to_subtract_from, scan_to_subtract); 
                    case 'alpha'  
                        subtract_data = self.regress_or_divide_csf_mrs(analysis, mask, 'alpha',  ...
                                                                       scan_to_subtract_from, scan_to_subtract);                                                               
                    case 'divide_crlb' 
                        subtract_data = self.regress_or_divide_crlb(analysis, mask, 'divide',  ...
                                                                    scan_to_subtract_from, scan_to_subtract);
                    case 'regress_crlb'
                         subtract_data = self.regress_or_divide_crlb(analysis, mask, 'regress',   ...
                                                                     scan_to_subtract_from, scan_to_subtract);
                    case 'lcmodel_csf'
                        subtract_data = self.regress_or_divide_csf_mrs(analysis, mask, 'lcmodel_csf',  ...
                                                                       scan_to_subtract_from, scan_to_subtract);
                        
                    case 'no_control'

                    otherwise
                        error('control not specified correcty');
                end

                cnt = cnt + 1;
                mrs_to_cor.('data')(:, cnt) = subtract_data; 
                mrs_to_cor.('labels'){cnt} = analysis;
            end
        end
       
        
        function [subtract_data] = regress_out_linewidth(self, analysis, mask, scan_to_subtract_from, scan_to_subtract)
            
            loc_or_ips = strsplit(analysis, '_');
            loc_or_ips = loc_or_ips{2};
                        
            s_from_LW = self.get_linewidth_data(loc_or_ips, mask, scan_to_subtract_from);
            s_sub_LW = self.get_linewidth_data(loc_or_ips, mask, scan_to_subtract);
            
            s_from_MRS = self.mrs.(analysis){mask, scan_to_subtract_from};
            s_sub_MRS = self.mrs.(analysis){mask, scan_to_subtract};
            
            subtract_LW = s_from_LW - s_sub_LW;
            subtract_MRS = s_from_MRS - s_sub_MRS;
            
            subtract_data = self.bvm.regress_data(subtract_LW, subtract_MRS, 'not_robust');  
            
        end

        
        function [subtract_data] = regress_or_divide_crlb(self, analysis, mask, regress_or_divide,  ...
                                                          scan_to_subtract_from, scan_to_subtract)
            %
            % Regress or divide out the Cramer-Rao Lower Bounds to control for data quality before
            % subtracting the scans.
            %
            % It makes most sense to regress out the change in data quality from the change in GABA
            % rather than regress out of a specific scan
            % ---------------------------------------------------          
            
            s_from_CRLB = self.get_crlb(analysis, mask, scan_to_subtract_from);
            s_sub_CRLB = self.get_crlb(analysis, mask, scan_to_subtract);
            
            s_from_MRS = self.mrs.(analysis){mask, scan_to_subtract_from};
            s_sub_MRS = self.mrs.(analysis){mask, scan_to_subtract};

            
            subtract_MRS = s_from_MRS - s_sub_MRS;
            subtract_CRLB = s_from_CRLB - s_sub_CRLB;
            
            if strcmp(regress_or_divide, 'divide')
                subtract_data = subtract_MRS ./ subtract_CRLB;

            elseif strcmp(regress_or_divide, 'regress')    
                subtract_data = self.bvm.regress_data(subtract_CRLB, subtract_MRS, 'not_robust');  

            end

            subtract_data = s_from_MRS - s_sub_MRS;

        end
        
        
    function [subtract_data] = regress_or_divide_csf_mrs(self, analysis, mask, divide_or_regress,  ... 
                                                             scan_to_subtract_from, scan_to_subtract)         
            %
            % Handle CSF controls from subtracted session data. First divide, regress or use alpha method to get
            % tissue corrected GABA values from sessions of interest. Then subtract these sessions to get the difference.
            %
            % -------------------------------------------------------------------------------------------------

            % get tissue and mrs data
            [s_from_GM, s_from_WM, s_from_CSF] = self.get_tissue_percentage(analysis, mask, scan_to_subtract_from);
            [s_sub_GM, s_sub_WM, s_sub_CSF] = self.get_tissue_percentage(analysis, mask, scan_to_subtract);
                                                         
            s_from_MRS = self.mrs.(analysis){mask, scan_to_subtract_from};
            s_sub_MRS = self.mrs.(analysis){mask, scan_to_subtract};
                      
            switch divide_or_regress
                
                case 'divide'  % TODO: remove, this is identical to lcmodel_csf
                    s_from_MRS = s_from_MRS ./ (s_from_GM + s_from_WM);
                    s_sub_MRS = s_sub_MRS ./ (s_sub_GM + s_sub_WM);

                case 'regress'                         
                    s_from_MRS = self.bvm.regress_data(s_from_CSF, s_from_MRS, 'not_robust');  
                    s_sub_MRS = self.bvm.regress_data(s_sub_CSF, s_sub_MRS, 'not_robust');       
                                  
                case 'alpha'
                    s_from_MRS = self.calculate_alpha(s_from_MRS, s_from_GM, s_from_WM); 
                    s_sub_MRS = self.calculate_alpha(s_sub_MRS, s_sub_GM, s_sub_WM);
                   
                case 'lcmodel_csf'
                    s_from_MRS = s_from_MRS ./ (1 - s_from_CSF); 
                    s_sub_MRS = s_sub_MRS ./ (1 - s_sub_CSF);
                    
                otherwise
                    error(sprintf('wrong divide, regress or alpha input: %s, must be divide or regress or alpha\n', divide_or_regress));
            end

            subtract_data = s_from_MRS - s_sub_MRS;

    end
           
    function [corrected_GABA] = calculate_alpha(~, MRS, GM, WM)
        %
        % Calculate alpha as specified in Harris et al., et al., 2016,
        % Tissue correction for GABA-edited MRS: considerations of voxel composition, tissue segmentation and tissue relaxations
        %
        % MRS: concentration of metabolite (e.g. GABA)
        % GM: tissue proportion of grey matter in voxel
        % WM: tissue proportion of white matter in voxel
        % -----------------------------------------------------------------
        a = 0.5;
        u_GM = mean(GM);
        u_WM = mean(WM);

        corrected_GABA = MRS .* ( 1 ./ (GM + a .* WM) ) .* ( (u_GM + a .* u_WM) ./ (u_GM + u_WM));
                
    end
        
% --------------------------------------------------------------------------------------------------------------------------------------------------------------
%       Helper functions
% --------------------------------------------------------------------------------------------------------------------------------------------------------------

        function [crlb] = get_crlb(self, analysis, mask, scan)
            %
            % Get CRLB for the specified analysis type, masked subjects and scan. 
            %
            % This function is a workaround to compensate bt some poorly named header
            % named in the old raw GABA MRS file where the mrs data is saved.
            %
            %
            %
            % First the analysis string is reformatted to batch the analysis format
            % on the raw gaba data structure. Then the Scan number is reformatted
            % to match the scan name on the table.
            %
            % TODO: go back and refactor the original tables
            % -----------------------------------------------------------------------

            raw_mrs_data = load(fullfile(self.config.mrs_csv_path, 'raw_mrs_data.mat'));
            split_fields = strsplit(analysis, '_');                            % ignore gaba and water suffix
            out_field = strcat('out_', split_fields{2}, '_', split_fields{1}); % CRLB for both gaba and water are in single table

            switch scan

                case 'Scan_1'
                    crlb_key = 'S1_CRLB';
                case 'Scan_2'
                    crlb_key = 'S2_CRLB';
                case 'Scan_3'
                    crlb_key = 'S3_CRLB';
                otherwise
                    error('wrong scan number specified')
            end

             crlb = raw_mrs_data.mrs_data.(out_field){mask, {crlb_key}};

        end
        
        function [linewidth] = get_linewidth_data(self, loc_or_ips, mask, scan)
            %
            %
            % ----------------------------------------------------------------
            
            linewidth_data = load(fullfile(self.config.mrs_csv_path, 'linewidth_data.mat'));
            
            linewidth = linewidth_data.linewidth_data.(loc_or_ips){mask, scan};
            
            
        end
        
        
        function [gm, wm, csf] = get_tissue_percentage(self, analysis, mask, scan)
            %
            % Get tissue percentage of voxel from saved structure.
            % Some awkward formatting for the different headers. 
            % TODO: Next time keep naming and organisation consistent!
            % ------------------------------------------------------------
            scan_num = extractAfter(scan, 'Scan');
            
            if contains(analysis, 'loc')
                region = 'LOC';
            elseif contains(analysis, 'ips')
                region = 'IPS';
            end

            gm = self.tissue_proportion.(region){strcat(mask, scan_num),  ...
                                                 'gray_matter'};                                                    
            wm = self.tissue_proportion.(region){strcat(mask, scan_num),  ...
                                                 'white_matter'};
            csf = self.tissue_proportion.(region){strcat(mask, scan_num),  ...
                                                  'CSF'}; 
        end
              
        
     end    
end
