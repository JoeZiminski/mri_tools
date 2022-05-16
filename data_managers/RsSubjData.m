classdef RsSubjData 
    %
    % Class for storing data for connectivity analysis. See
    % correlate_region_timeseries.m for generation of inputs.
    %
    % ---------------------------------------------------------------------
    
    properties
            
       plotter = Plotter; 
       
       num_regions = []; 
       labels = [];
       plot_labels = [];
       corr_values = [];
       correlation_matrices = struct();
       subj_id = [];
       ts_discard_idx = [];
       log = {};
       
    end
    
    methods
        
        
        function [self] = RsSubjData(num_regions, corr_values, plot_labels, subj_id, ts_discard_idx, log)
            
             self.num_regions = num_regions; 
             self.plot_labels = plot_labels;
             self.corr_values = corr_values;    
             self.subj_id = subj_id;
             self.correlation_matrices = self.make_correlation_matrices();
             self.ts_discard_idx = ts_discard_idx;
             self.log = log;

             
        end
        
                
        function [corr_mat] = make_correlation_matrices(self)
            %
            % Reshape the subjects data into a region x region correlation
            % matrix
            %
            % -------------------------------------------------------------
            corr_mat =  struct('Scan_1', [],  ...
                               'Scan_2', [],  ...
                               'Scan_3', []);
            
            
            for scan = ["Scan_1", "Scan_2", "Scan_3"]
                
               data = self.corr_values{:, scan};
               corr_matrix = reshape(data, self.num_regions, self.num_regions);
               
               corr_mat.(scan) = corr_matrix;
                
            end
           
            
        end
        
                
        function plot_correlation_matrix(self, scans, varargin)
            
            if ~isempty(varargin)
                caxis_ = varargin{1};

            else
                caxis_ = [-1 1];
            end
            
            if strcmp(scans, 'all')
                scans = ["Scan_1", "Scan_2", "Scan_3"];
            end
            
            for scan = scans
                
                 self.plotter.make_heatmap(self.correlation_matrices.(scan),  ...
                                           self.plot_labels,  ...
                                           strcat(self.subj_id, '_', scan),  ...
                                           caxis_)
                  
            end
            
        end
       
    end
      
end
