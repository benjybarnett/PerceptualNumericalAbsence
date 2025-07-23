function clusterCorrection(cfg)

% function clusterCorrection(cfg)
% takes in uncorrected pvalues and transforms these into an image of
% significant clusters based on the set q-threshold and required cluster
% size
%          cfg.pvalImage = a map with uncorrected p-values
%          cfg.qvalue = the sifnificance cut-of, will use FDR correction
%          at this level
%          cfg.clustersize = required cluster size for cut-off
%          cfg.df = degrees of freedom (N-1)

cfg.mask = fullfile('/vol/ccnlab1/naddij/Templates','wmni_icbm152_t1_tal_nlin_sym_09a_mask.nii');

[dir,~,~] = fileparts(cfg.pvalImage); % get the directory to put the files in

[V,Y] = read_nii(cfg.pvalImage); % read in the p-values

% get the mask
[~,brainmask] = read_nii(cfg.mask);

% mask the pvalues and write to file
pVal = brainmask.*Y; pMap = fullfile(dir,'masked_pval.nii');
write_nii(V,pVal,pMap);

% do the fdr correction

cfg.inputfile = pMap; % use the masked image for further processing
FDR = bb_fdr(cfg); % get the FDR value

% get the cluster table

cfg.threshold       = 1-FDR;
cfg.atlas           = '/vol/ccnlab-scratch1/naddij/fieldtrip-20140614/template/atlas/aal/ROI_MNI_V4.nii';
cfg.searchRadius    = 5;
cfg.outputfile      = fullfile(dir,'cluster_table.txt');
table = bb_cluster_table(cfg);

small_clusters = table(str2double(table(:,2)) < cfg.clustersize,1);

% seperate the clusters

bb_separate_clusters(cfg);
for c = 1:length(small_clusters) % delete maps small clusters
    delete(fullfile(dir,['masked_pval_cluster' small_clusters{c} '.nii']))
end

% put all together in one significance mask
clustermaps = str2fullfile(dir,'masked_pval_cluster*.nii');
clusters = zeros([length(clustermaps), V.dim]);

for c = 1:length(clustermaps)
    [~,clusters(c,:,:,:)] = read_nii(clustermaps{c});
end

sigmap = any(clusters,1);
write_nii(V,sigmap,fullfile(dir,sprintf('sig_FDR_%d_map.nii',round(cfg.qvalue),2)))



