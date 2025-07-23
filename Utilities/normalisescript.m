%% Set up path names etc.

addpath(genpath('/vol/ccnlab1/naddij/ImageryPerceptionfMRI/Analyses'))
addpath(genpath('/vol/ccnlab-scratch1/naddij/spm8'))
addpath /vol/optdcc/fieldtrip-latest/qsub
root = '/vol/ccnlab-scratch1/naddij/ImageryPerceptionfMRI/';

subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S16','S17','S18','S19','S20','S21','S23','S25','S26','S27','S28','S29'};


%% Do the normalisation
maps = 'spmT_0*.img';
NormaliseVolumes(subjects, root, '/DCM/Univariate',maps)

%% Load the normalised pictures and average
dim = [79 95 68];
nsubjects = length(subjects);
norm_brain_corr = zeros(nsubjects,dim(1),dim(2),dim(3));


for sub = 1:nsubjects
    
    subjectname = subjects{sub};
    
    searchlight_dir = [root subjectname '/Searchlight'];
    
    norm_brain = str2fullfile(searchlight_dir,'wfaces_corr.nii');
    hdr = spm_vol(norm_brain);
    norm_brain_corr(sub,:,:,:) = spm_read_vols(hdr);
end




cd([root '/GroupResults/Searchlight'])
norm_braincorr = squeeze(nanmean(norm_brain_corr,1));
write_nii(hdr,norm_braincorr,'faces_braincorr.nii')


% check some things
subnans = zeros(size(norm_brain_corr));
nannormbrain = isnan(norm_brain_corr);

renormbrain = reshape(nannormbrain,nsubjects,prod(dim));
nvoxels = size(renormbrain,2);
x = zeros(nvoxels,1);

for v = 1:nvoxels
    x(v,1) = sum(double(renormbrain(:,v)));
end


% do a one sided t-test
rmpath(genpath('/vol/ccnlab-scratch1/naddij/spm8'))

renormbrain = reshape(norm_brain_corr,nsubjects,prod(dim));
nvoxels = size(renormbrain,2);

h = zeros(nvoxels,1);
p = zeros(nvoxels,1);

for v = 1:nvoxels
    if  isempty(find(isnan(renormbrain(:,v)),1)) % no nans
        [h(v,1),p(v,1)] = ttest(renormbrain(:,v));
    end
end

% write the results
addpath(genpath('/vol/ccnlab-scratch1/naddij/spm8'))

normbrain_sig = reshape(h,dim(1),dim(2),dim(3));

write_nii(hdr,normbrain_sig,'norm_brainsig.nii')
