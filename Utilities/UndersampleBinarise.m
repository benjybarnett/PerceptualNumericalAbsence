function [data,labels] = UndersampleBinarise(cfg0,data,labels,tClass)
    % Function to undersmpale classes before binarisation. 
    % Target class indicated by tClass argument.
    % Results in non-target classes all having same number of samples
    % when binarised into target vs non-target.
    nclasses = max(labels); 
    
    % Sample count for each class
    N = arrayfun( @(c) sum(labels==c) , 1:nclasses);
    
    
    % undersample the majority class(es)
    rm_samples = abs(N - min(N));
    for cc=1:nclasses
        if rm_samples(cc)>0
            ix_this_class = find(labels == cc);
            ix_rm = randperm( numel(ix_this_class), rm_samples(cc));
    
            % Remove samples from data and labels
            
            data(ix_this_class(ix_rm),:,:)= [];
            labels(ix_this_class(ix_rm))= [];
        end
    end
    
    N = arrayfun( @(c) sum(labels==c) , 1:nclasses);

    % Now cut down the non-target classes to 1/5th of number of target trials
    nClassNzero = round(N(tClass)/cfg0.numNTClass); 
    rm_samples = N - nClassNzero;
    for cc = 1:nclasses
        if cc == tClass
            continue
        end
        ix_this_class = find(labels == cc);
        ix_rm = randperm( numel(ix_this_class), rm_samples(cc));
        data(ix_this_class(ix_rm),:,:)= [];
        labels(ix_this_class(ix_rm))= [];
        disp(arrayfun( @(c) sum(labels==c) , 1:nclasses));
        
    end

end