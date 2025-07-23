function bb_fsl_smooth(cfg)
% bb_fsl_smooth(cfg)
%
% Uses FSL machinary to smooth an image. Full-Width Half Maximum is
% coverted to sigma of gaussian kernel.
%
% cfg.inputfile     = filename of image to be smoothed
% cfg.outputfile    = output filename of smoothed image
% cfg.FWHM          = Full-Width Half Maximum, which is coverted to sigma 
%                       of gaussian kernel by this function

% load variables from cfg
get_vars_from_struct(cfg)


%%
sigmaGauss = FWHM / (2*sqrt(2*log(2)));

fslCommand = ['fslmaths ' inputfile ' -s ' num2str(sigmaGauss) ' ' outputfile];

unix(fslCommand);
