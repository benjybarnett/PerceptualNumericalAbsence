function trialMatrix = pivotTrialsDots(trialInfo)
    
    trialMatrix = zeros(size(trialInfo,1),8);
    trialMatrix(:,1) = trialInfo(:,1);
    trialMatrix(:,2) = trialInfo(:,2);
    trialMatrix(:,3) = trialInfo(:,3);
    trialMatrix(:,4) = repmat([1;2],size(trialInfo,1)/2,1); %sample or test stim
    trialMatrix(:,6) = trialInfo(:,27);
    trialMatrix(:,7) = trialInfo(:,30);
    trialMatrix(:,8) = trialInfo(:,31);
    trialMatrix(:,9) = trialInfo(:,32);
    trialMatrix(:,10) = trialInfo(:,33);
    count = 1;
    for row = 1:size(trialInfo,1)
        if count ==1
            trialMatrix(row,5) = trialInfo(row,28);
            count = 2;
        else
            trialMatrix(row,5) = trialInfo(row,29);
            count = 1;
        end
       

    end

end