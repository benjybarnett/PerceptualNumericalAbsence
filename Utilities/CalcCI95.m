function CIs = CalcCI95(data)
%% Function to calculate 95% confidence intervals of group data
% Input: Data = [Subjects x Time Points] data matrix
% Output: CIs = [1 x Time Points] array of 95% confidence intervals
    
    std_dev = std(data,1,'omitnan');
    CIs = zeros(1,size(data,2));
    for i =1:size(data,2)
        sd = std_dev(i);
        n = size(data,1);
        
        CIs(i) = 1.96*(sd/sqrt(n));
    end
    

end