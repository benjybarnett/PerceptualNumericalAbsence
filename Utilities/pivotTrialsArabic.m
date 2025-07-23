function trialMatrix = pivotTrialsArabic(trialInfo)
    
    trialMatrix = zeros(size(trialInfo,1),8);
    trialMatrix(:,1) = trialInfo(:,1);
    trialMatrix(:,2) = trialInfo(:,2);
    trialMatrix(:,3) = trialInfo(:,3);
    trialMatrix(:,6) = trialInfo(:,24);
    trialMatrix(:,7) = trialInfo(:,25);
    trialMatrix(:,8) = trialInfo(:,26);
    count = 1;
    for row = 1:size(trialInfo,1)
        if count <= 10
            trialMatrix(row,4) = trialInfo(row,3+count);
            trialMatrix(row,5) = trialInfo(row,13+count);
            count = count+1;
        else
            trialMatrix(row,4) = trialInfo(row,4);
            trialMatrix(row,5) = trialInfo(row,14);
            count = 2;
        end
       

    end

end