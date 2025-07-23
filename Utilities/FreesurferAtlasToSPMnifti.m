function FreesurferAtlasToSPMnifti(cfg)
% function FreesurferAtlasToSPMnifti(cfg)
%
% cfg.root    = root directory where subjects directory spring from
% cfg.subject = subjectID (e.g. 'S01'), which is the name of the subject
% directory in the root and also in freesurfer's output
% cfg.nifti   = name of directory where nifti's are, amongst which
% resliced_struct.nii and mean_func.nii
% cfg.atlas   = the name of the atlas to convert, e.g. 'aparc.DKTatlas40'
%
% creates a SPM_to_Freesurfer.dat transformation matrix that allows
% transformation between the two spaces. 

subjectDir = fullfile(cfg.root,cfg.subject);
fsreconDir = fullfile(cfg.root,'Freesurfer');

labelDir = fullfile(subjectDir,'labels'); % where do you want the labels
if ~exist(labelDir,'dir'); mkdir(labelDir); end
    
% convert annotation to labels
annotation = cfg.atlas;
hemi = 'lh'; % for left hemisphere
unix(['mri_annotation2label --subject ' cfg.subject ' --sd ' fsreconDir ...
    ' --annotation ' annotation ' --hemi ' hemi ' --labelbase '...
    fullfile(labelDir,hemi)])
hemi = 'rh'; % for right hemisphere
unix(['mri_annotation2label --subject ' cfg.subject ' --sd ' fsreconDir ...
    ' --annotation ' annotation ' --hemi ' hemi ' --labelbase '...
    fullfile(labelDir,hemi)])

% create transformation func to FS space
source_image = fullfile(subjectDir,cfg.nifti,'mean_func.nii'); % mean functional
target_image = fullfile(fsreconDir,cfg.subject,'mri/orig.mgz');
funct_to_FS = fullfile(fsreconDir,cfg.subject,'SPM_to_Freesurfer.dat');
unix(['tkregister2 --mov ' source_image ' --targ ' target_image...
    ' --noedit --s ' cfg.subject     ' --sd ' fsreconDir ' --regheader --reg ' funct_to_FS])

% transform the labels to volumes
labels = str2fullfile(labelDir,'*.label');
nlabels = numel(labels);

for l = 1:nlabels
    labelfile = labels{l};
    unix(['mri_label2vol --label ' labelfile ' --temp ' source_image ' --reg ' funct_to_FS ...
        ' --o ' [labelfile '.nii.gz']])    
end

% make one volume containing all rois
labelvols = str2fullfile(labelDir,'*.nii.gz');
V = read_nii(labelvols{1});
mask = nan(V.dim);

for l = 1:nlabels
    [~,Y] = read_nii(labelvols{l});
    mask(Y == 1) = l;        
end

write_nii(V,mask,fullfile(labelDir,'ROImask.nii'))
save(fullfile(labelDir,'ROImask'),'mask','labels')
    