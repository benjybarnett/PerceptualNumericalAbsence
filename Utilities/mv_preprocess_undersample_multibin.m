function [pparam, X, clabel] = mv_preprocess_undersample_multibin(pparam, X, clabel)
% For use when decoding a multiclass problem using a series 
% of One Vs. All classifiers. First, undersamples the majority class(es) from the non-target 
% classes. Then, binarises the non-target and target classes and undersamples the majority class
% again. Works out as 
%
%Usage:
% [pparam, X, clabel] = mv_preprocess_undersample(pparam, X, clabel)
%
%Parameters:
% X              - [samples x ... x ...] data matrix
% clabel         - [samples x 1] vector of class labels
%
% pparam         - [struct] with preprocessing parameters
% .sample_dimension - which dimension(s) of the data matrix represent the samples
%                     (default 1)
% .undersample_test_set - by default, if undersampling is used during
%                         cross-validation, only the train set is
%                         undersampled (default 0). If the test data has to
%                         be undersampled, too, set to 1.
%
% Undersampling can safely be performed globally (on the full dataset)
% since it does not introduce any dependencies between samples.

if pparam.is_train_set || pparam.undersample_test_set

    sd = sort(pparam.sample_dimension(:))';
    nclasses = max(clabel);
    ixs_rm = []; % to hold all the indexes of X to remove

    % Sample count for each class
    N = arrayfun( @(c) sum(clabel==c) , 1:nclasses);
    % there can be multiple sample dimensions. Therefore, we build a colon
    % operator to extract the train/test samples irrespective of the
    % position and number of sample dimensions
    s = repmat({':'},[1, ndims(X)]);

    % undersample the majority class(es) in the non-target classes
    non_T_classes = N(unique(clabel) ~= pparam.target_class);
    rm_nonT_samples = abs(non_T_classes - min(non_T_classes));
    classes = unique(clabel);
    nonT_labels = classes(classes ~= pparam.target_class);
    for cc = 1:length(nonT_labels)
        class = nonT_labels(cc);
        if rm_nonT_samples(cc)>0
            ix_this_class = find(clabel == class);
            ix_rm = randperm( numel(ix_this_class), rm_nonT_samples(cc));
            %ixs_rm = [ixs_rm ix_rm];
           
             % Remove samples from all sample dimensions
            for rm_dim=sd
                s_dim = s;
                s_dim(rm_dim) = {ix_this_class(ix_rm)};
                X(s_dim{:})= [];
            end
            clabel(ix_this_class(ix_rm))= [];

        end
    end
    
    %Create a binary list of target vs. non target classes
    bin_labels = (clabel ~= pparam.target_class) + 1; 
    %Sample count for binary classes
    bin_N = arrayfun(@(c) sum(bin_labels == c), 1:2);
    %Samples to remove to balance binary target/non-target classes
    rm_bin_samples = abs(bin_N - min(bin_N));


    % undersample the majority class(es) in binary target vs non-target setup
    if bin_N(1) > bin_N(2) %if target class is the majority
        ix_this_class = find(bin_labels == 1);
        ix_rm = ix_this_class(randperm(rm_bin_samples(1)))';
        ixs_rm = [ixs_rm ix_rm];
       
    elseif bin_N(2) > bin_N(1) %if non-target classes in majority
        rm_per_class = floor(rm_bin_samples(2) / (nclasses-1));
        for cc = 1:length(nonT_labels)
            class = nonT_labels(cc);
            
            ix_this_class = find(clabel == class);
            ix_rm = ix_this_class(randperm(rm_per_class))';
            ixs_rm = [ixs_rm ix_rm];
        
        end
    end
  


    bin_labels(ixs_rm)= []; %remove undersampled trials
    clabel = bin_labels; %return binarised labels

    % Remove samples from all sample dimensions
    for rm_dim=sd
        s_dim = s;
        s_dim(rm_dim) = {ixs_rm};
        X(s_dim{:})= [];
    end

    
elseif ~ pparam.is_train_set
    bin_labels = (clabel ~= pparam.target_class) + 1;
    clabel = bin_labels;
end

