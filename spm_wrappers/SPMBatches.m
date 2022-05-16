
classdef SPMBatches
    %
    % To run, ensure fMRIBatches is on path.
    % initialise e.g.:  batches = fMRIBatches;. 
    % Functions may then be accessed e.g:
    %   [SPMbatch] = batches.realign_and_unwarp();
    %
    %
    % Note: to save a batch to add to the class, load into Workspace and click 'Save as' > 'MATLAB Script'
    % -----------------------------------------------------------------------------------------------------
    
    methods (Static)
        
% ----------------------------------------------------------------------------------------------------------
% General Batches (Used in both fMRI and MPM)
% ----------------------------------------------------------------------------------------------------------
        
        
        function [SPMbatch] = convert_3D_to_4D()
            
            SPMbatch.matlabbatch{1}.spm.util.cat.vols = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.util.cat.name = '4D.nii';
            SPMbatch.matlabbatch{1}.spm.util.cat.dtype = 0;  % keep same data format
            
        end
        
        function [SPMbatch] = convert_4D_to_3D()
            
            SPMbatch.matlabbatch{1}.spm.util.split.vol = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.util.split.outdir = {''};           
            
        end
        
        function [SPMbatch] = dicom_to_nifti_hmri()
            
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.data = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.root = 'flat';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.outdir{1} = '';   % current working directory
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.protfilter = '.*';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.convopts.format = 'nii';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.convopts.metaopts.mformat = 'sep';  % 'ext' for header 
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.convopts.icedims = 0;
            
        end
        
        
        function [SPMbatch] = cat12_longitundinal_alignment()
            
            SPMbatch.matlabbatch{1}.spm.tools.cat.tools.series.data = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.cat.tools.series.noise = NaN;
            SPMbatch.matlabbatch{1}.spm.tools.cat.tools.series.bparam = 1.0E+6;
            SPMbatch.matlabbatch{1}.spm.tools.cat.tools.series.use_brainmask = 0;
            SPMbatch.matlabbatch{1}.spm.tools.cat.tools.series.reduce = 0;
            SPMbatch.matlabbatch{1}.spm.tools.cat.tools.series.reg.rigid = 1;
            SPMbatch.matlabbatch{1}.spm.tools.cat.tools.series.write_rimg = 0;
            SPMbatch.matlabbatch{1}.spm.tools.cat.tools.series.write_avg = 1;
        
        end
    

        function [SPMbatch] = spm_normalise_from_deformation()                                

            SPMbatch.matlabbatch{1}.spm.spatial.normalise.write.subj.def = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.write.subj.resample = '<UNDEFINED>';  % file to resampel
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70;
                                                                                78 76 85];
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;            % trilinear interpolation for brain images
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'mni_';
                        
        end
        
     
% ----------------------------------------------------------------------------------------------------------
% fMRI Batches
% ----------------------------------------------------------------------------------------------------------        

        function [SPMbatch] = realign_estimate()

            SPMbatch.matlabbatch{1}.spm.spatial.realign.estimate.data = {'<UNDEFINED>'};
            SPMbatch.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
            SPMbatch.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
            SPMbatch.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
            SPMbatch.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 1;
            SPMbatch.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
            SPMbatch.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
            
        end


        function [SPMbatch] = realign_and_unwarp()
                                           
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.data.scans = '<UNDEFINED>';    % all scans (put SBRef first)
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan = '';              % vdm5 for unwarping
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 0];
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 3;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
            SPMbatch.matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'ur_';
            
        end
        
        
        function [SPMbatch] = coregistration_estimate_and_reslice_epi(varargin)
            
            if ~isempty(varargin)
                cost_fun = varargin{1};
            else
                cost_fun = 'nmi';
            end
            
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.ref = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.source = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.other{1} = '';
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r_';
            
        end
        
        
        function [SPMbatch] = realign_reslice()
            
            SPMbatch.matlabbatch{1}.spm.spatial.realign.write.data = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 1];
            SPMbatch.matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
            SPMbatch.matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 0;
            SPMbatch.matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';

        end
        
        
        function [SPMbatch] = coregistration_reslice()
            
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.write.ref = '<UNDEFINED>';    % imaging defining space
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.write.source = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

        end
        
        
        function [SPMbatch] = coregistration_estimate_only_epi(varargin)
            
            if ~isempty(varargin)
                cost_fun = varargin{1};
            else
                cost_fun = 'nmi';
            end

            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estimate.ref = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estimate.source = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estimate.other{1} = '';
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = cost_fun;
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            SPMbatch.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

        end
        
        
        function [SPMbatch] = normalise_and_write()
            
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.subj(1).vol = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.subj(1).resample = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'D:\spm12\tpm\TPM.nii'};
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                                                   78 76 85];
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
            SPMbatch.matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

        end
        
        
        function [SPMBatch] = normalise_write()
            % 
            % see SPM GUI for argumetns
            % 0 - nearest neighbour
            % 4 - 4th degree B-Spline
            %
            
            SPMBatch.matlabbatch{1}.spm.spatial.normalise.write.subj = struct('def', {}, 'resample', {});
            SPMBatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                      78 76 85];
            SPMBatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
            SPMBatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
            SPMBatch.matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
            
        end
        
        
        function [SPMbatch] = spm_segment(spm_path)
        
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.channel.vols = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm_path,'tpm', 'TPM.nii,1')};
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm_path,'tpm', 'TPM.nii,2')};
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm_path,'tpm', 'TPM.nii,3')};
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm_path,'tpm', 'TPM.nii,4')};
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm_path,'tpm', 'TPM.nii,5')};
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm_path,'tpm', 'TPM.nii,6')};
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
            SPMbatch.matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                                                   NaN NaN NaN];

        end
                       
        
        function [SPMbatch] = cat12_segment(spm_path)
            % see spm GUI for option descriptions
            % --------------------------------------
            
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.data = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;                            % option to run parallel
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {fullfile(spm_path, '/tpm/TPM.nii')};
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.dartel.darteltpm = {fullfile(spm_path, '/toolbox/cat12/templates_1.50mm/Template_1_IXI555_MNI152.nii')};
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.best = [0.5 0.1];  % use the best native approach for 7T. Change this to use defaults
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
            % SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;  % don't write ROI atlas out
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;  % save the native space
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
            SPMbatch.matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0];                      
            
        end
        
        
        function [SPMbatch] = skull_strip_anatomical()
            % i1 - 3 are tissue segmentation results (c1-3)
            % i4 is EPI to skullstrip
            SPMbatch.matlabbatch{1}.spm.util.imcalc.input = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.util.imcalc.output = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.util.imcalc.outdir{1} = '';
            SPMbatch.matlabbatch{1}.spm.util.imcalc.expression = '((i1+i2+i3)>.5).*i4';  
            SPMbatch.matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            SPMbatch.matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            SPMbatch.matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            SPMbatch.matlabbatch{1}.spm.util.imcalc.options.interp = -2;  % second order sinc interpolation
            SPMbatch.matlabbatch{1}.spm.util.imcalc.options.dtype = 16;

        end
        
        
        function [SPMbatch] = skull_strip_epi()
            
            SPMbatch.matlabbatch{1}.spm.util.imcalc.input = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.util.imcalc.output = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.util.imcalc.outdir{1} = '';
            SPMbatch.matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
            SPMbatch.matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            SPMbatch.matlabbatch{1}.spm.util.imcalc.options = struct;
            SPMbatch.matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            SPMbatch.matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            SPMbatch.matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            SPMbatch.matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            
        end

        
        function [SPMbatch] = smooth_epi()
           
            SPMbatch.matlabbatch{1}.spm.spatial.smooth.data = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.spatial.smooth.fwhm = [4 4 4];
            SPMbatch.matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            SPMbatch.matlabbatch{1}.spm.spatial.smooth.im = 0;
            SPMbatch.matlabbatch{1}.spm.spatial.smooth.prefix = 'sm_';
            
        end
        

        function [SPMbatch] = first_level_glm()
            
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.dir = '<UNDEFINED>';                                % where to write SPM.mat file
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';                             % TR units (num scans 'scans' vs. seconds 'secs')
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = '<UNDEFINED>';                           % TR in seconds
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;                                 % The microtime resolution. See https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;70bf5c7f.0904
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;                                 % no slice-scan timing correction was performed so leave the default values - this and above value make no difference for RS data
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.sess.scans = '<UNDEFINED>';                         % path to scans for this session
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, ...     % conditions, ignore for RS
                                                                           'duration', {}, ...
                                                                           'tmod', {}, 'pmod', {}, ...
                                                                           'orth', {});                                                
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};                               % multiple conditions, ignore for RS       
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});    % regressors
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''} ;                          % multiple regressors (option to include mat file of regressors)
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;                                  % high-pass filter
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});            % factorial design    
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];                           % include basis function derivatives (e.g. for cannonical HRF)
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;                                           % Volterra interactions (1 for off, 2 for on)
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';                                    % global normalisation
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.mthresh = '<UNDEFINED>';                         % explicit mask threshold. Set to zero for RS and 0.8 (SPM default) for other? see akiraoconnor.org/2010/04/07/masking-in-spm
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.mask = {''};                                       % explicing mask - brain mask suggested
            SPMbatch.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';                                   % autocorrelation correction method
            
        end
        
        function [SPMbatch] = rs_estimate()
            % DOC
            SPMbatch.matlabbatch{1}.spm.stats.fmri_est.spmmat = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            SPMbatch.matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                        
        end         

        
        function [SPMbatch] = rs_contrast()   
            % DOC
            SPMbatch.matlabbatch{1}.spm.stats.con.spmmat = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'test_con';
            SPMbatch.matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 0;
            SPMbatch.matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            SPMbatch.matlabbatch{1}.spm.stats.con.delete = 0;
        end
        
        
% ----------------------------------------------------------------------------------------------------------
% MPM Batches
% ----------------------------------------------------------------------------------------------------------        
        
        function [SPMbatch] = hmri_dicom_to_nii_with_metadata()
            
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.data = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.root = 'flat';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.outdir{1} = '';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.protfilter = '.*';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.convopts.format = 'nii';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.convopts.metaopts.mformat = 'sep';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.dicom.convopts.icedims = 0;
            
        end

        
        function [SPMbatch] = hmri_auto_reorient()
            
            SPMbatch.matlabbatch{1}.spm.tools.hmri.autoreor.reference = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.autoreor.template = {'D:\spm12\canonical\avg152T1.nii'};
            SPMbatch.matlabbatch{1}.spm.tools.hmri.autoreor.other = {''};
            SPMbatch.matlabbatch{1}.spm.tools.hmri.autoreor.output.indir = 'yes';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.autoreor.dep = 'individual';

        end
        

        function [SPMbatch] = hmri_normalisation(num_maps, varargin)
            
            % number of image types to normalise (one channel per class e.g. MT, R2 etc. )
            for map = 1:num_maps
                SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.mni_norm.multsdata.vols_pm{map} = '<UNDEFINED>';
            end
            
            % in dir vs. out dir
            if ~isempty(varargin)
                SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.mni_norm.output.outdir = varargin(1);
            else
                SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.mni_norm.output.indir = 1;
            end
            
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.mni_norm.template = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.mni_norm.multsdata.vols_tc{1} = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.mni_norm.multsdata.vols_tc{2} = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.mni_norm.multsdata.vols_field = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.mni_norm.vox = [NaN NaN NaN];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.mni_norm.bb = [NaN NaN NaN;
                                                                                            NaN NaN NaN];
           
        end
        
        function [SPMbatch] = hmri_segmentation(hmri_path, varargin)
            
            % in dir vs. out dir
            if ~isempty(varargin{1})
                SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.output.outdir = varargin(1);
            else
                SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.output.indir = 1;
            end
            
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.channel.vols = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.vols_pm{1} = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.vox = [1 1 1];                % Voxel Dimenions and  
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.bb = [-78 -112 -70;           % bounding box
                                                                                              78 76 85];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.channel.biasreg = 0;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.channel.biasfwhm = Inf;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.many_sdatas.channel.write = [0 0];         % Bias field corrected images
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(1).tpm{1} = fullfile(hmri_path, 'etpm',' eTPM.nii,1');
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(1).ngaus = 2;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(1).native = [1 1];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(1).warped = [1 1];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(2).tpm{1} = fullfile(hmri_path, 'etpm', 'eTPM.nii,2');
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(2).ngaus = 2;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(2).native = [1 1];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(2).warped = [1 1];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(3).tpm{1} = fullfile(hmri_path, 'etpm', 'eTPM.nii,3');
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(3).ngaus = 2;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(3).native = [1 1];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(3).warped = [1 1];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(4).tpm{1} = fullfile(hmri_path, 'etpm', 'eTPM.nii,4');
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(4).ngaus = 3;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(4).native = [0 0];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(4).warped = [0 0];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(5).tpm{1} = fullfile(hmri_path, 'etpm', 'eTPM.nii,5');
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(5).ngaus = 4;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(5).native = [0 0];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(5).warped = [0 0];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(6).tpm{1} = fullfile(hmri_path, 'etpm', 'eTPM.nii,6');
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(6).ngaus = 2;
            SPMbatch. matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(6).native = [0 0];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.tissue(6).warped = [0 0];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.warp.mrf = 1;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.warp.cleanup = 1;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.warp.reg = [0 0.001 0.5 0.05 0.2];
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.warp.affreg = 'mni';
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.warp.fwhm = 0;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.warp.samp = 3;
            SPMbatch.matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_us.warp.write = [1 1];  % forward and backward deformations

        end

% ----------------------------------------------------------------------------------------------------------
% DARTEL Batches
% ----------------------------------------------------------------------------------------------------------  

        function [SPMbatch] = make_dartel_template()

            % in dir vs. out dir
            if ~isempty(varargin{1})
                matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.output.outdir =  varargin(1);
            else
                matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.output.indir = 1;
            end
            
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.images{1} = '<UNDEFINED>';
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.images{2} = '<UNDEFINED>';

            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.template = 'Template';
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.rform = 0;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(1).its = 3;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(1).rparam = [4 2 1E-6];
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(1).K = 0;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(1).slam = 16;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(2).its = 3;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(2).rparam = [2 1 1E-6];
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(2).K = 0;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(2).slam = 8;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(3).its = 3;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(3).rparam = [1 0.5 1E-6];
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(3).K = 1;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(3).slam = 4;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(4).its = 3;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(4).rparam = [0.5 0.25 1E-6];
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(4).K = 2;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(4).slam = 2;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(5).its = 3;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(5).rparam = [0.25 0.125 1E-6];
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(5).K = 4;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(5).slam = 1;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(6).its = 3;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(6).rparam = [0.25 0.125 1E-6];
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(6).K = 6;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.param(6).slam = 0.5;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.optim.lmreg = 0.01;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.optim.cyc = 3;
            matlabbatch{1}.spm.tools.hmri.proc.proc_modul.proc_dart.warp.settings.optim.its = 3;

        end


        function [SPMbatch] = spm_flow_field_to_deformation(varargin) 
            
            if ~isempty(varargin)
                SPMbatch.matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = varargin(1);
            else
                SPMbatch.matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.savepwd = 1;  % save on current directory
            end
            
            SPMbatch.matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
            SPMbatch.matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
            SPMbatch.matlabbatch{1}.spm.util.defs.comp{1}.dartel.template{1} = '';
            SPMbatch.matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = '<UNDEFINED>';
            
        end
        
        
% -------------------------------------------------------------------------------------------------------------------------------------------------
% General Analyis Batches
% -------------------------------------------------------------------------------------------------------------------------------------------------


        function [SPMbatch] = one_way_repeat_measures_ANOVA()

            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.dir = '<UNDEFINED>';                                             % Ouput dir
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(1).scans = '<UNDEFINED>';                    % Path to Scans 
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(1).conds = '<UNDEFINED>';                    % vector of condition numbers
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.des.anovaw.dept = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.des.anovaw.variance = 1;  % 1 = unequal variance, 0 = equal variance
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.des.anovaw.ancova = 0;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

        end
    
        
        function [SPMbatch] = second_level_one_sample_w_covariate()
            
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.dir = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.cov.c = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.cov.cname = '<UNDEFINED>';
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            SPMbatch.matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
        end
        
    end
end




