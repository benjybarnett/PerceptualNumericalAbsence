% Delete files from folders when you have to rerun analyses

root = '/vol/ccnlab-scratch1/naddij/DCM/';

subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S16','S17','S18','S19','S20','S21','S23','S25','S26','S27','S28','S29'};
nsubjects = length(subjects);

for sub = 1:nsubjects
    
    fprintf('Deleting files for subject %s \n',subjects{sub})
    
    %VOIs = str2fullfile(fullfile(root,subjects{sub},'GLM/CollCat'),'VOI*');
    %for v = 1:length(VOIs)
    %    delete(VOIs{v});
    %end
    
    DCMs = str2fullfile(fullfile(root,subjects{sub},'GLM/CollCat'),'DCM_PercDirInd.mat');
    delete(DCMs)
    %rmdir(fullfile(root,subjects{sub},'PPI/V1_perception'),'s')
    

end