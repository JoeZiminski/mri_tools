classdef SPMMethods
    %
    % Concept for class to run SPMMethods in. Seems to work well 
    % (better than running in script) but it would be even better to have
    % each spm function in a standalone class that contains its batch
    % and run method. See 7T_TUNING project. This minimises code in the 
    % main running scripts
    % ---------------------------------------------------------------------

    properties  
        
        batches = [];
    
    end
    
    
    methods
        
        function self = SPMMethods()
            
           self.batches = SPMBatches;                   
            
        end
        
        function init_fmri_defaults(~)
            
           spm('defaults', 'fmri');
           spm_jobman('initcfg'); % initialise SPM for fMRI analysis
           
        end
        
                  
        function dcm_2_nii_file(self, input_filepath, output_filepath)
           % inputs must be char/ str
           % --------------------------------------------------------------
            
           dcms = dir(fullfile(input_filepath, '*.dcm'));
           dcm_paths = fullfile({dcms.folder}, {dcms.name})';
           
           if ~isfolder(output_filepath) 
               mkdir(output_filepath);
           end
           
           
           SPMbatch = self.batches.dicom_to_nifti_hmri();
           
           SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.data = dcm_paths;
           SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.convopts.metaopts.mformat = 'ext';
           SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.outdir = {output_filepath};
           
           spm_jobman('run', SPMbatch.matlabbatch); 
            
        end
        
        function normalise_and_write(self, input_filepath, vox)
            
            SPMBatch = self.batches.normalise_and_write();
            
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = vox;
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.subj(1).vol = {input_filepath};
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.subj(1).resample = {input_filepath};
            
            spm_jobman('run', SPMbatch.matlabbatch); 
            
        end
        
        function normalise_write(self, def_filepath, anat_filepath, interp_method, vox)
            
            SPMbatch = self.batches.normalise_write();
            
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.write.subj = struct('def', {{def_filepath}}, 'resample', {{anat_filepath}});
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vox;
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = interp_method;
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
            
            spm_jobman('run', SPMbatch.matlabbatch); 
            
        end
        
        function hmri_segmentation(self, hmri_path, nii_filepath, output_path, vox_dims, bounding_box)
            
           SPMbatch = self.batches.hmri_segmentation(hmri_path, output_path);
           
           SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.vox = vox_dims;
           SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.bb = bounding_box;
           SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.channel.vols = {nii_filepath};
           SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.vols_pm = {};
            
           spm_jobman('run', SPMbatch.matlabbatch); 
            
        end
        
        
        function run_cat12_segmentation(self, spm_path, input_volume_path)
        
            SPMbatch = self.batches.cat12_segment(spm_path);
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.data = input_volume_path;
            
            spm_jobman('run', SPMbatch.matlabbatch); 
            
        end
            
        function run_spm_one_way_repeat_measures_ANOVA(self, output_path, subject_nii_paths, conditions, varargin)
            %
            % subject_nii_paths is a 1 x n cell array of subject filepaths e.g. subjects{1} = 3x1 char array of paths for three scans in order (e.g. session 1, 2 3)
            % conditions is the condition numbers for the scans, assumed to be the same across all subjects (e.g. [1 2 3]
            % -----------------------------------------------------------------------------------------------------------------------------------------------
            
            mask_path = varargin;
            
            % Generate and run batch
            SPMbatch = self.batches.one_way_repeat_measures_ANOVA();
            SPMbatch.matlabbatch{1,  1}.spm.stats.factorial_design.dir = {output_path}; 
            
            for s = 1:length(subject_nii_paths)
                SPMbatch.matlabbatch{1, 1}.spm.stats.factorial_design.des.anovaw.fsubject(s).scans  = subject_nii_paths{s};          % Session X filename  array per subj         
                SPMbatch.matlabbatch{1, 1}.spm.stats.factorial_design.des.anovaw.fsubject(s).conds = conditions;
            end
            
            if isempty(mask_path)
                SPMbatch.matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
                SPMbatch.matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            else
                SPMbatch.matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
                SPMbatch.matlabbatch{1}.spm.stats.factorial_design.masking.em = mask_path;
            end
            spm_jobman('run', SPMbatch.matlabbatch(1));
            
            % Estimate model
            SPMbatch.matlabbatch{1,  2}.spm.stats.fmri_est.spmmat = {fullfile(output_path, 'SPM.mat')};
            spm_jobman('run', SPMbatch.matlabbatch(2)); 
        
        end
        
        
        function run_spm_second_level_one_sample_w_covariate(self, output_path, subject_nii_paths, covariates, covariate_name, varargin)
            %
            % -----------------------------------------------------------------------------------------------------------------------
            
            mask_path = varargin;
            
            % run second level analysis 
            SPMbatch = self.batches.second_level_one_sample_w_covariate();
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.dir = {output_path};
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = subject_nii_paths;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.cov.c = covariates;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.cov.cname = covariate_name;
            
            if isempty(mask_path)
                SPMbatch.matlabbatch{1,  1}.spm.stats.factorial_design.masking.im = 1;
                SPMbatch.matlabbatch{1, 1}.spm.stats.factorial_design.masking.em = {''};
            else
                SPMbatch.matlabbatch{1,  1}.spm.stats.factorial_design.masking.im = 0;
                SPMbatch.matlabbatch{1,  1}.spm.stats.factorial_design.masking.em = mask_path;
            end
            spm_jobman('run', SPMbatch.matlabbatch);

            % estimate
            SPMbatch = struct();
            SPMbatch.matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(output_path, 'SPM.mat')};
            spm_jobman('run', SPMbatch.matlabbatch);

            % contrast
            SPMbatch = struct();
            SPMbatch.matlabbatch{1}.spm.stats.con.spmmat = {fullfile(output_path, 'SPM.mat')};
            SPMbatch.matlabbatch{1}.spm.stats.con.consess{1,  1}.tcon.name = 'positive';
            SPMbatch.matlabbatch{1}.spm.stats.con.consess{1,  1}.tcon.weights = [0 1];
            SPMbatch.matlabbatch{1}.spm.stats.con.consess{1,  2}.tcon.name = 'negative';
            SPMbatch.matlabbatch{1}.spm.stats.con.consess{1,  2}.tcon.weights = [0 -1];
            spm_jobman('run', SPMbatch.matlabbatch);

        end
        
        
    end
    
end



