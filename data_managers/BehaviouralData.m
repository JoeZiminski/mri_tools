classdef BehaviouralData
    
   properties
       
       config = struct();
       behav = struct();
       bvm = BivariateMethods;
        
   end
   
   
   methods
       
        % constructer
        function self = BehaviouralData(config)

            self.config = config;   
            self.behav = behaviour_data_prepro_master;

        end     


        function [behav_to_cor] = get_subtract_behav(self, fields_to_get, mask, control, varargin)
        %
        % Get difference in behavioural data across sessions. 
        %
        % INPUT:
        %       fields_to_get: behavioural data to retrieve (naming based on structure output from 
        %                      'behaviour_data_prepro_master'
        %                      More than one can be given e.g. ["avg_percent_correct", "learning_rate_1_6"]
        %       mask: cell array of subjects to include (e.g. {'subject_01', 'subject_02'}
        %       control: perform control on the data, 'divide_s1' or 'regress_s1'
        %       varargin: by default Scan_3 - Scan_2 is performed. Enter scan to subtract from and scan to subtract
        %                 if want to subtract other (see handle_scan_to_subtract_varargin)
        % 
        %       e.g. BehaviouralData.get_subtract_behav(["avg_percent_correct"],   ...
        %                          config.masks.learners, 'divide_s1', 'Scan_2', 'Scan_1');
        %
        % -------------------------------------------------------------------------------------

        behav_to_cor = struct('data',   [NaN(length(mask),  length(fields_to_get))],  ...
                              'labels', {cell(1, length(fields_to_get))}); 

        [scan_to_subtract_from, scan_to_subtract] = handle_scan_to_subtract_varargin(varargin);

        if ismember('Scan_1', {scan_to_subtract, scan_to_subtract_from}) &&  ...
                ismember(control, {'clean_data', 'regress_s1', 'divide_s1'})
                error('Due to regression of S1, cleaning data is only valid for S3-S2 currently');                          
        end
        
        
        cnt = 0;
        for analysis = string(fields_to_get)

            subtract_data = subtract_sessions(self.behav.(analysis)(mask, :),  ...
                                              scan_to_subtract_from, scan_to_subtract);

            % handle controls
            if ~contains(analysis, 'learning_rate')

                s1 = self.behav.(analysis){mask, 'Scan_1'};
                switch control
                    case 'divide_s1'
                        subtract_data = subtract_data ./ s1;
                    case 'regress_s1'
                        subtract_data = self.bvm.regress_data(s1, subtract_data, 'not_robust');
                end
            end
            
            cnt = cnt + 1;
            behav_to_cor.('data')(:, cnt) = subtract_data; 
            behav_to_cor.('labels'){cnt} = analysis;

        end
    end
    
   end
end


   