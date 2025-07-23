function SPMNiftitoFreesurfer(cfg)
% function SPMNiftitoFreesurfer(cfg)
%
% cfg.root    = root directory where subjects directory spring from
% cfg.subject = subjectID (e.g. 'S01'), which is the name of the subject
% directory in the root and also in freesurfer's output
% cfg.nifti   = name of directory where nifti's are, amongst which
% resliced_struct.nii 
% cfg.map     = name of the map to transform - from cfg.root/cfg.subject
% cfg.space   = 'subject' or 'mni'
%
% creates a SPM_to_Freesurfer.dat transformation matrix that allows
% transformation between the two spaces. 

fsreconDir = fullfile(cfg.root,'Freesurfer');
if strcmp(cfg.space,'subject')
    subjectDir = fullfile(cfg.root,cfg.subject);
end

% create transformation func to FS space
if strcmp(cfg.space,'subject')
    source_image = fullfile(subjectDir,cfg.nifti,'resliced_struct.nii'); % resliced structural
    target_image = fullfile(fsreconDir,cfg.subject,'mri/orig.mgz');
elseif strcmp(cfg.space,'mni')
    cfg.subject = 'Group';  
    
    source_image = '/vol/ccnlab1/naddij/Templates/MNI/mni_t1.nii';    
    target_image = fullfile(fsreconDir,cfg.subject,'mri/orig.mgz');
end

funct_to_FS = fullfile(fsreconDir,cfg.subject,'SPM_to_Freesurfer.dat');

% set subject dir
unix(['export SUBJECTS_DIR=' fsreconDir]);

unix(['tkregister2 --mov ' source_image ' --targ ' target_image...
    ' --noedit --s ' cfg.subject     ' --sd ' fsreconDir ' --regheader --reg ' funct_to_FS])


% transform the map to surfaces
if strcmp(cfg.space,'mni')
    map = fullfile(cfg.root,cfg.map);
elseif strcmp(cfg.space,'subject')
    map = fullfile(subjectDir,cfg.map);
end
[~,name] = fileparts(cfg.map);

outputFileLH = fullfile(fsreconDir,cfg.subject,['lh_' name '.w']);

unix(['mri_vol2surf --src ' map ' --srcreg ' funct_to_FS ...
    ' --out ' outputFileLH ' --out_type paint  --hemi lh'])

outputFileRH = fullfile(fsreconDir,cfg.subject,['rh_' name '.w']);

unix(['mri_vol2surf --src ' map ' --srcreg ' funct_to_FS ...
    ' --out ' outputFileRH ' --out_type paint  --hemi rh'])
