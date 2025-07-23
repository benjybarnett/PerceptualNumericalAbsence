function clusterTable = bb_cluster_table(cfg)
% clusterTable = bb_cluster_table(cfg)
%
% cfg.inputfile       = statisical image (e.g. t-map)
% cfg.pImage          = image containing (un)corrected p-values, usually in
%                           the 1-p format for easy thresholding
% cfg.threshold       = t-value or (1 minus) p-value (in this case also 
%                           specify cfg.pImage)
% cfg.df              = integer, degrees of freedom (e.g. N-1) for conversion
%                          T-to-P-toZ conversion (only for parametric
%                          stats!)
% cfg.atlas           = string, can be 'aal', 'afni' or a filename of a
%                           atlas file supported by ft_read_atlas
% cfg.searchRadius    = radius of sphere in which anatomical label may be
%                           searched
% cfg.outputfile      = tab delim txt file where the table is stored

% load variables from cfg
get_vars_from_struct(cfg)


if ~isfield(cfg,'atlas') || isempty(cfg.atlas)
    findLabels = false;
else
    findLabels = true;
end 

if ~isfield(cfg,'searchRadius') || isempty(cfg.searchRadius)
    searchRadius = 0;
end

if isfield(cfg,'outputfile') && ~isempty(cfg.outputfile) 
    writeOutput = true;
else
    writeOutput = false;
end

if isfield(cfg,'pImage') && ~isempty(cfg.pImage) 
    
else

end



%%

fslCmd = sprintf('cluster -i %s --thresh=%f --mm',cfg.inputfile,cfg.threshold);

[~, clusterOutput] = unix(fslCmd);

disp(clusterOutput)

cellTable        = textscan(clusterOutput,'%s %s %s %s %s %s %s %s %s');
cellTable        = horzcat(cellTable{:});
cellTable(1:3,:) = [];
cellTable(:,7:9) = [];

tVals  = cellfun(@str2num,cellTable(:,3));

%[pVals zVals] = t2p2z(tVals,df);
pVals = tVals; zvals = tVals;

mniCoord = cellfun(@str2num,cellTable(:,4:6));


%% reorder table into following format:
% 1: cluster number
% 2: number of voxels
% 3: X coord mni
% 4: Y coord mni
% 5: Z coord mni
% 6: peakZ
% 7: peakT
% 8: peakP
% 9: anatomical labels

clusterTable = cellTable(:,[1 2 4 5 6 3]);

nRois = size(clusterTable,1);

if findLabels
    
    switch atlas
        case 'aal'
            atlasFn = 'D:\bbarnett\Documents\ecobrain\fieldtrip-master-MVPA\template\atlas\aal\ROI_MNI_V4.nii';
        case 'afni'
            atlasFn = '/vol/ccnlab-scratch1/naddij/fieldtrip-20140614/template/atlas/afni/TTatlas+tlrc.BRIK';
        otherwise
            atlasFn = atlas;
    end
    
    addpath('/vol/ccnlab-scratch1/naddij/fieldtrip-20140614')
    addpath('/vol/ccnlab-scratch1/naddij/fieldtrip-20140614')
    atlasStruct = ft_read_atlas(atlasFn);
    
    labels = atlas_label_from_coord(mniCoord,atlasStruct.tissue,atlasStruct.tissuelabel,atlasStruct.transform,searchRadius);
else
    labels = cell(nRois,1); % empty cells
end

clusterTable(:,9:18) = {[]};

for iRoi = 1:nRois
    
    clusterTable{iRoi,6} = sprintf('%.2f',zVals(iRoi));
    clusterTable{iRoi,7} = sprintf('%.2f',tVals(iRoi));
    clusterTable{iRoi,8} = sprintf('%.2e',pVals(iRoi));
    
    if ~isempty(labels{iRoi})
        
        nRegions = length(labels{iRoi}{1});
        clusterTable(iRoi,9:9+nRegions-1) = labels{iRoi}{1};
    else
        clusterTable{iRoi,9} = '';
    end
end

display(clusterTable)


%% write output to csvfile

if writeOutput
    
    fid = fopen(outputfile,'w');
    
    fprintf(fid,'clusterId\tnVox\tmniX\tmniY\tmniZ\tpeakZ\tpeakT\tpeakP\tlabel\n');
    
    for iRoi = 1:nRois
        fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',clusterTable{iRoi,:});
    end
    fclose(fid);
    
end

