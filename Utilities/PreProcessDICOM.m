function PreProcessDICOM(cfg)
% function PreProcessDICOM(cfg)
%
% reads in DICOM files and performs Realignment, Coregistration and
% Segmentation in that order. 'cfg' is a configuration structure with the
% following fields:
% cfg.
%     dicom_dir = directory of the dicom files
%     nifti_dir = where the niftis should be written
%     run_nr = optional, which runs to include

[~,subjectname] = fileparts(fileparts(fileparts(cfg.nifti_dir)));

if ~exist(cfg.nifti_dir,'dir') % if the map doesn't exist, make it
    mkdir(cfg.nifti_dir);
end

if ~exist(fullfile(cfg.nifti_dir,'dicom_conversion_done'),'file') % check if converted dicoms do not already exist
    %% Dicom import
    dicoms = str2fullfile(cfg.dicom_dir,'*.IMA');
    
    if isempty(dicoms)
        folders = dir(cfg.dicom_dir); folders = {folders.name};
        
        if iscell(cfg.dicom_dir)
            for s = 1:size(cfg.dicom_dir,2)
                if isfield(cfg,'run_nr') && ~isempty(cfg.run_nr{s})
                    for r = 1:size(cfg.run_nr{s},2)
                        dicoms = [dicoms,str2fullfile([cfg.dicom_dir{s} '/' folders{strncmp(folders,sprintf('%03d',1),3)}],'*.IMA')];%sprintf('*SKYRA.%04d*.IMA',cfg.run_nr{s}(r)))];
                    end
                else
                    dicoms = [dicoms,str2fullfile(cfg.dicom_dir{s},'*.IMA')];
                end
                
            end
        else
            if isfield(cfg,'run_nr') && ~isempty(cfg.run_nr)
                for r = 1:size(cfg.run_nr,2)
                    dicoms = [dicoms,str2fullfile([cfg.dicom_dir '/' folders{strncmp(folders,sprintf('%03d',cfg.run_nr(r)),3)}],'*.IMA')];%[dicoms,str2fullfile(cfg.dicom_dir,sprintf('*SKYRA.%04d*.IMA',cfg.run_nr(r)))];
                end
            else
                dicoms = [dicoms;str2fullfile(cfg.dicom_dir,sprintf('*%s*.IMA',subjectname))];
            end
        end
    end
    
    dicom{1}.spm.util.dicom.data = dicoms';
    dicom{1}.spm.util.dicom.root = 'flat';
    dicom{1}.spm.util.dicom.outdir = {cfg.nifti_dir};
    dicom{1}.spm.util.dicom.convopts.format = 'nii';
    dicom{1}.spm.util.dicom.convopts.icedims = 0;
    
    fprintf('Running dicom import for subject %s\n',subjectname)
    
    spm_jobman('run',dicom)
    unix(sprintf('touch %s', fullfile(cfg.nifti_dir,'dicom_conversion_done')))
    clear dicoms
else
    fprintf('dicom conversion for subject %s already done\n',subjectname)
end

%% Realignment

if ~exist(fullfile(cfg.nifti_dir,'realignment_done'),'file')
    
    niftis = str2fullfile(cfg.nifti_dir,'f*.nii'); % functional nifti's
    
    realign{1}.spm.spatial.realign.estwrite.data = {niftis'};
    realign{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    realign{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    realign{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    realign{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    realign{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    realign{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    realign{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    realign{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    realign{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    realign{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    realign{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    realign{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    fprintf('running realignment for subject %s\n',subjectname)
    
    spm_jobman('run',realign)
    clear niftis
    
    unix(sprintf('touch %s', fullfile(cfg.nifti_dir,'realignment_done')))
    
else
    fprintf('realignment for subject %s already done \n',subjectname)
end

%% Coregistration
if ~exist(fullfile(cfg.nifti_dir,'coregistration_done'),'file')
    
    mean_file = str2fullfile(cfg.nifti_dir,'meanf*.nii');
    mean_coreg_file = fullfile(cfg.nifti_dir,'mean_coreg.nii');
    
    % copy mean file to be used for coregistration estimation
    if ~exist(mean_coreg_file, 'file')
        copyfile(mean_file,mean_coreg_file);
    end
    
    real_scans = str2fullfile(cfg.nifti_dir,'rf*.nii');
    
    if contains(cfg.nifti_dir,'ses1') % structural is in the first session
        structural = str2fullfile(cfg.nifti_dir, 's*128-01*.nii'); %
    else
        structural = str2fullfile(fullfile(fileparts(fileparts(...
            cfg.nifti_dir)),'ses1','Niftis'),'str*.nii');
    end
    
    coreg{1}.spm.spatial.coreg.estimate.ref = {structural};
    coreg{1}.spm.spatial.coreg.estimate.source = {mean_coreg_file};
    coreg{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    coreg{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    coreg{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    coreg{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    fprintf('running coregistration for subject %s\n',subjectname)
    
    spm_jobman('run',coreg)
    
    % SPM coregistration parameters
    flags.estimate.cost_fun = 'nmi';
    flags.estimate.sep = [4 2];
    flags.estimate.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    flags.estimate.fwhm = [7 7];
    
    % create transformation matrix
    coregParameters = spm_coreg(spm_vol(mean_coreg_file), ...
        spm_vol(mean_file), ...
        flags.estimate);
    EPItransformationMatrix = spm_matrix(coregParameters);
    
    % apply transformation matrix to functionals
    
    % add mean functional to the list
    real_scans{end+1} = mean_file;
    
    for iFunc = 1:size(real_scans,2)
        if mod(iFunc,100)==0
            fprintf('\tworking on file %d\n',iFunc)
        end
        iFuncSpace = spm_get_space(char(real_scans(iFunc)));
        spm_get_space(char(real_scans(iFunc)), EPItransformationMatrix\iFuncSpace);
    end
    
    unix(sprintf('touch %s',fullfile(cfg.nifti_dir,'coregistration_done')))
    clear real_scans
else
    fprintf('coregistration for subjects %s already done \n',subjectname)
    
end

%% Segmentation
if ~exist(fullfile(cfg.nifti_dir,'segmentation_done'),'file')...
        && contains(cfg.nifti_dir,'ses1')
    
    segm{1}.spm.spatial.preproc.data = {structural};
    segm{1}.spm.spatial.preproc.output.GM = [0 0 1];
    segm{1}.spm.spatial.preproc.output.WM = [0 0 1];
    segm{1}.spm.spatial.preproc.output.CSF = [0 0 0];
    segm{1}.spm.spatial.preproc.output.biascor = 1;
    segm{1}.spm.spatial.preproc.output.cleanup = 0;
    segm{1}.spm.spatial.preproc.opts.tpm = {
        '/vol/optdcc/spm12/toolbox/OldSeg/grey.nii'
        '/vol/optdcc/spm12/toolbox/OldSeg/white.nii'
        '/vol/optdcc/spm12/toolbox/OldSeg/csf.nii'
        };
    segm{1}.spm.spatial.preproc.opts.ngaus = [2
        2
        2
        4];
    segm{1}.spm.spatial.preproc.opts.regtype = 'mni';
    segm{1}.spm.spatial.preproc.opts.warpreg = 1;
    segm{1}.spm.spatial.preproc.opts.warpco = 25;
    segm{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
    segm{1}.spm.spatial.preproc.opts.biasfwhm = 60;
    segm{1}.spm.spatial.preproc.opts.samp = 3;
    segm{1}.spm.spatial.preproc.opts.msk = {''};
    
    spm_jobman('run',segm)
    
    fprintf('running segmentation for subject %s\n',subjectname)
    
    grey_matter_mask = str2fullfile(cfg.nifti_dir, 'c1s*128-01*.nii');
    white_matter_mask = str2fullfile(cfg.nifti_dir, 'c2s*128-01*.nii');
    
    % reslice to functional space (note that mean_func is already coregistered to the space of the segmentations
    
    % SPM reslice parameters
    resliceParameters = struct(...
        'prefix', 'r',...
        'mask', 1,...
        'interp', 4, ...
        'wrap', [0 0 0], ...
        'which', [2 1]);
    
    resFlags = struct(...
        'interp', resliceParameters.interp,... % interpolation type
        'wrap', resliceParameters.wrap,...     % wrapping info (ignore...)
        'mask', resliceParameters.mask,...     % masking (see spm_reslice)
        'which', 1,...                  % what images to reslice
        'mean', 0);                     % write mean image
    
    spm_reslice({mean_file, grey_matter_mask}, resFlags);
    spm_reslice({mean_file, white_matter_mask}, resFlags);
    spm_reslice({mean_file, structural}, resFlags);
    
    unix(sprintf('touch %s', fullfile(cfg.nifti_dir,'segmentation_done')))
else
    fprintf('segmentation of subject %s already done \n',subjectname)
    
end

%% Rename key files for ease of use later
if contains(cfg.nifti_dir,'ses1')

new_struct_name = fullfile(cfg.nifti_dir, 'struct.nii');
new_resliced_struct = fullfile(cfg.nifti_dir, 'resliced_struct.nii');
new_resliced_gm = fullfile(cfg.nifti_dir, 'resliced_gm.nii');
new_resliced_wm = fullfile(cfg.nifti_dir, 'resliced_wm.nii');
new_mean_func = fullfile(cfg.nifti_dir, 'mean_func.nii');

if ~exist(new_mean_func,'file')
    
    mean_file = str2fullfile(cfg.nifti_dir,'meanf*.nii');
    mean_coreg_file = fullfile(cfg.nifti_dir,'mean_coreg.nii');
    structural = str2fullfile(cfg.nifti_dir, 's*128-01*.nii'); % even checken!
    resliced_structural = str2fullfile(cfg.nifti_dir, 'rs*128-01*.nii'); % even checken!
    resliced_grey_matter_mask = str2fullfile(cfg.nifti_dir, 'rc1*128-01*.nii');
    resliced_white_matter_mask = str2fullfile(cfg.nifti_dir, 'rc2*128-01*.nii');
    
    fprintf('renaming key files for subject %s\n',subjectname)
    movefile(resliced_structural,new_resliced_struct);
    movefile(resliced_grey_matter_mask,new_resliced_gm);
    movefile(resliced_white_matter_mask,new_resliced_wm);
    movefile(structural,new_struct_name);
    movefile(mean_file,new_mean_func);
    delete(mean_coreg_file);
end

end
%% Functional masking
if ~exist(fullfile(cfg.nifti_dir,'funcmasking_done'),'file')
    fprintf('running funcmasking for subject %s \n',subjectname)
    
    % create functional mask out of gm and wm
    hdr = spm_vol(fullfile(fullfile(fileparts(fileparts(...
            cfg.nifti_dir)),'ses1','Niftis'),'resliced_gm.nii'));
    gm = spm_read_vols(hdr);
    hdr = spm_vol(fullfile(fullfile(fileparts(fileparts(...
            cfg.nifti_dir)),'ses1','Niftis'),'resliced_wm.nii'));
    wm = spm_read_vols(hdr);
    
    func_brain_mask = gm > 0.5 | wm > 0.5;
    
    write_nii(hdr,func_brain_mask,fullfile(cfg.nifti_dir,'func_brain_mask.nii'))
    clear wm gm func_brain_mask
    
    % use the functional mask to mask functionals
    real_scans = str2fullfile(cfg.nifti_dir,'rf*.nii');
    
    [~, maskY] = read_nii(fullfile(cfg.nifti_dir,'func_brain_mask.nii'));
    
    nrscans = size(real_scans,2);
    
    for iFunc = 1:nrscans
        
        if mod(iFunc,100)==0
            fprintf('\tworking on file %d\n',iFunc)
        end
        
        [funcV, funcY] = read_nii(real_scans{iFunc});
        newY = maskY .* funcY;
        
        write_nii(funcV,newY,funcV.fname);
    end
    
    unix(sprintf('touch %s', fullfile(cfg.nifti_dir,'funcmasking_done')))
else
    fprintf('funcmasking for subject %s already done \n',subjectname)
    
end

% change access settings
unix(sprintf('chmod 770 -R -f %s',cfg.nifti_dir))