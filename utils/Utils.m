classdef Utils
     
    methods
              
       
        
        function [subjects, sessions] = get_subjects_and_sessions(self, config, num_subjects, num_sessions)
            % Convenience function to return subjects and a list of
            % sessions from configs.
            % OUTPUT: subjects: a cell array of subject names (e.g. {subject_01, subject_02}
            %         sessions: a cell array of session names (e.g. {subject_01_1, subject_01_2, subject_01_3, 
            %                                                       {subject_02_1, subject_02_2, subject_02_3};
            % --------------------------------------------------------------------------------------------
            subjects = config.subjects.Properties.RowNames;
            subjects = subjects(num_subjects);
            sessions = self.make_list_of_sessions(subjects, num_sessions);
            
        end
        
        
        function [list_of_sessions] = make_list_of_sessions(~, subjects, num_sessions)
        % make a list of sessions for every subject. If a session is missing for a
        % subject this script will need to be refactored. 
        % -------------------------------------------------------------------------

        cnt = 0;
        for isubj = 1:length(subjects)
            for isesh = 1:length(num_sessions)

                cnt = cnt + 1;
                list_of_sessions{cnt, 1} =  strcat(subjects{isubj}, '_', num2str(isesh));

            end
        end
            
        end
        
        
        function [average_upper_tri] = get_average_upper_triangular(~, matrix)
                    
               average_upper_tri =  sum(sum(triu(matrix, 1))) / sum(1:size(matrix, 1)-1); 
            
        end
                
        function [upper_triangular_no_diag] = get_upper_triangular_no_diag(~, matrix)            

             matrix(matrix == 0 ) = eps;
             matrix = triu(matrix, 1);
             matrix(matrix == 0) = [];
             matrix(matrix == eps) = 0;

             upper_triangular_no_diag = matrix;
            
        end
        
        
    end
    
        
    methods (Static)
        
        function [transposed_table] = transpose_table(table)
            
            cell_table = table2cell(table);
            trans_cell_table = cell_table';
            
            transposed_table = cell2table(trans_cell_table, 'RowNames', table.Properties.VariableNames,  ...
                                                            'VariableNames', table.Properties.RowNames);
        end
        
        function copy_files(files_to_copy, new_base_dir, varargin) 
            %
            % Copy files from an old base directory to a new base directory
            % 
            % if no varargin, files_to_copy should be the output of dir()
            %
            % if an original base dir is provided in varargin, files_to_copy
            % should be a pathtable
            %
            % ------------------------------------------
            
            if ~isempty(varargin)
                base_dir = varargin{1};
                
                for irun = 1:height(files_to_copy)

                    copyfile(cjoin([{base_dir}, files_to_copy{irun, :}], filesep),   ...
                             cjoin([{new_base_dir}, files_to_copy{irun, :}], filesep));
                         
                end
                                    
            else
            
                for irun = 1:length(files_to_copy)
                    
                    copyfile(fullfile(files_to_copy(irun).folder, files_to_copy(irun).name),  ...
                             fullfile(new_base_dir, files_to_copy(irun).name));
                end
            end
        end
        
        
        function delete_files(cell_array_of_paths, varargin)  
            % Take a cell array of paths and delete.
            % if structure in form returned from dir is passed, will turn
            % into cell array of paths to delete.
            % if varargin is supplied, it is assumed input is a pathtable
            % and cell array of paths made from the pathtable and base dir
            % passed in the varargin
            %  -----------------------------------------------------------------------------------------------
            
            if isstruct(cell_array_of_paths)
                cell_array_of_paths = fullfile({cell_array_of_paths.folder},  ...
                                               {cell_array_of_paths.name});
            end
            
            if ~isempty(varargin)
                assert(istable(cell_array_of_paths), sprintf('pathtable must be first argument if varargin supplied'));
                cell_array_of_paths = fullfile(varargin{1}, join(cell_array_of_paths{:, :}, filesep));
            end
            
            for ipath = 1:length(cell_array_of_paths)
                
                path_ = cell_array_of_paths(ipath);
                delete(char(path_));
          
            end
            
        end
        
% ----------------------------------------------------------------------------------------------------------
% Pathtable Methods 
% ----------------------------------------------------------------------------------------------------------

  
        function varargout = parallel_copy(pathtable_to_copy, orig_base_dir, new_base_dir, varargin)  
        % Copy files parrallised file-by-file. Requires input files to be
        % at the same level of the directory tree. ~30% faster than a straight copy.
        %
        %
        % INPUT: pathtable_to_copy: a table of filepaths where each row is a
        %                           seperate path and each column contains one level of the directory tree
        %                           above the base_dir. The final column should be the file name.
        %                           i.e. 'subject'    'session'  'filename'
        %                            'S01'        'S01_1'      'epi.nii'
        %       orig_base_dir: base directory (which pathtable_to_copy dirs are on) to
        %                  be copeied from.
        %       new_base_dir: directory to be copied to
        %
        %       varargin:  optional argument to rename files. All files copied will
        %                 be prefixed with char found in vargargin e.g. 'coreg_'.
        %                 Varargout returns a new pathtable with the new filenames
        %                 in place of the old. 
        %  -----------------------------------------------------------------------------------------------


            if isempty(gcp('nocreate'))
                error('No active parpool detected')
            end

            % rename destination files if varargin exists
            destination_pathtable = pathtable_to_copy;
            if ~isempty(varargin)  % test for varargin
                prefix = char(varargin{1});

                for irun = 1:height(destination_pathtable)
                    filename = destination_pathtable{irun, end};
                    new_name = strcat(prefix, filename);
                    destination_pathtable{irun, end} = new_name;
                end
                varargout{1} = destination_pathtable;
            end

            % copy files
            parfor irun = 1:height(pathtable_to_copy) 
                new_path  =  char(join([{new_base_dir}, pathtable_to_copy{irun, 1:end-1}], ...   % make new base dir
                                        filesep));
                if ~exist(new_path, 'dir')
                    mkdir(new_path)
                end

                copyfile(char(join([{orig_base_dir}, pathtable_to_copy{irun, :}], ... 
                                     filesep)),  ...
                         char(join([{new_path}, destination_pathtable{irun, end}], ...
                                    filesep)));
            end
            
        end
              
        
        function move_files(pathtable_to_move, orig_base_dir, new_base_dir)
            
            for irun = 1:height(pathtable_to_move)
            
                new_path = char(join([new_base_dir, pathtable_to_move{irun, 1:end-1}], filesep));
                
                if ~exist(new_path)
                    mkdir(new_path);
                end
                
                movefile(char(join([{orig_base_dir}, pathtable_to_move{irun, :}], filesep)),  ...
                         char(join([{new_path}, pathtable_to_move{irun, end}], filesep)));
            end
            
        end
        
        
        function [filepaths] = make_paths(add_suffix, type, rename, base_path, varargin)
            % Make cell array inputs for SPM from a pathtable and base directory(see make_path_table).
            %
            % INPUT: add_suffix: if 'add_suffix' ',1' will be added to the
            %                    end of each filename for SPM
            %        type: 'cell' or 'char' to return as cell or char array
            %       
            %        rename: if not empty, filename will be prefixed with
            %                supplied string/char (requires pathtable input
            %                with 'filename' field. 
            %        base_path: base directory for the pathtable
            %       
            %        varargin: pathtable or cell array of paths to make filepaths from (these
            %                  will be vertically concatenated)
            %  -----------------------------------------------------------------------------------------------
            
            if ~isempty(rename)
                for irn = 1:length(varargin)
                    varargin{irn}{:, 'filename'} = strcat(rename, varargin{irn}{:, 'filename'});
                end
            end

            for iconv = 1:length(varargin)
                if istable(varargin{iconv})
                    varargin{iconv} = varargin{iconv}{:, :};  
                end
            end

            filepaths = fullfile(base_path,  ...
                                 join(vertcat(varargin{:}), filesep));

            if strcmp(add_suffix, 'add_suffix')
                for iname = 1:length(filepaths)
                    filepaths(iname) = strcat(filepaths(iname), ',1');
                end
            end

            if strcmp(type, 'char')
                filepaths = char(filepaths);
            end
        end

        
        function [path] = quick_path(base_path, pathtable)
            % Quick function to make a single path (char) from base path and
            % pathtable
            % -------------------------------------------------------
            
            % check inputs
            assert(height(pathtable) == 1, sprintf('Pathtable must only be 1 row in height\n'));
            if ~iscell(base_path)
                base_path = {base_path};
            end
            
            path = join([base_path, pathtable{1, :}],  ...
                   filesep);
               
            path = char(path);
            
        end
                
        
        function [pathtable] = slice_pathtable(pathtable, column, logic, string_, varargin)
            % Slice into a pathtable to return only cells that match
            % string. For multiple inputs, iteratively slice with each
            % input from first to last, whittling down the pathtable.
            %
            % e.g. slice_pathtable(pathtable, {'subj', 'sesh'}, {'exact', 'ends_with'}, {subject', '_1'});
            %
            %     will first slice out all the subjects matching subject, then
            %     all the sessions ending in 1.
            %
            % INPUT: pathtable: tabletable to slice
            %        
            %       column: column to match with in e.g. 'subj' to find
            %               cells that match 'subject_01'
            %
            %       logic: type of matching to perform (see switch
            %              statement below, includes regexp)
            %       string_: string to match column contents with logic
            %                conditional
            %       varargin: 'to_cell' or 'to_char' to return as cell/char 
            %                 array of paths rather than pathtable. 
            % ----------------------------------------------------------------------------------------------
           
            for run = 1:length(column)
                
                col = column{run};
                log = logic{run};
                str = string_{run};

                switch log
                    case 'starts_with'
                        pathtable = pathtable(startsWith(pathtable{:, col}, str), :);
                    case 'ends_with'
                        pathtable = pathtable(endsWith(pathtable{:, col}, str), :);
                    case 'contains' 
                        pathtable = pathtable(contains(pathtable{:, col}, str), :);
                    case '~starts_with'
                        pathtable = pathtable(~startsWith(pathtable{:, col}, str), :);
                    case '~ends_with'
                        pathtable = pathtable(~endsWith(pathtable{:, col}, str), :);
                    case '~contains' 
                        pathtable = pathtable(~contains(pathtable{:, col}, str), :);
                    case 'exact'
                        pathtable = pathtable(strcmp(pathtable{:, col}, str), :);
                    case 'regexp'
                        pathtable = pathtable(regexpcmp(pathtable{:, col}, str), :);
                    otherwise
                        error('logic not valid, options are starts_with, ends_with, exact, contains or regexp');
                end
            end    
            
            if strcmp(varargin, 'to_cell')
                pathtable = join(table2cell(pathtable), filesep);
            elseif strcmp(varargin, 'to_char')
                pathtable = cjoin(table2cell(pathtable), filesep);
            end

        end
        
        
    end
end
