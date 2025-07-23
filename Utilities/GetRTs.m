function [means,sd] = GetRTs(data)

    %Calculates the average RT for each distance between sample and test numerosity.
    %Only looks at correct trials
    
    %when sample and test are the same
    zeroMean = mean(data((data(:,28) - (data(:,29)) == 0)& (data(:,33) ==1),32));
    zeroSD =  std(data((data(:,28) - (data(:,29)) == 0)& (data(:,33) ==1),32));            
      
    %when the distance between sample and test is 1
    oneMean = mean(data((abs(data(:,28) - (data(:,29))) == 1)& (data(:,33) ==1),32));
    oneSD = std(data((abs(data(:,28) - (data(:,29))) == 1)& (data(:,33) ==1),32));

    twoMean = mean(data((abs(data(:,28) - (data(:,29))) == 2)& (data(:,33) ==1),32));
    twoSD = std(data((abs(data(:,28) - (data(:,29))) == 2)& (data(:,33) ==1),32));

    threeMean = mean(data((abs(data(:,28) - (data(:,29))) == 3)& (data(:,33) ==1),32));
    threeSD = std(data((abs(data(:,28) - (data(:,29))) == 3)& (data(:,33) ==1),32));

    fourMean = mean(data((abs(data(:,28) - (data(:,29))) == 4)& (data(:,33) ==1),32));
    fourSD = std(data((abs(data(:,28) - (data(:,29))) ==4)& (data(:,33) ==1),32));
    
    fiveMean = mean(data((abs(data(:,28) - (data(:,29))) == 5)& (data(:,33) ==1),32));
    fiveSD = std(data((abs(data(:,28) - (data(:,29))) ==5)& (data(:,33) ==1),32));

          
    means = [zeroMean,oneMean,twoMean, threeMean,fourMean,fiveMean];
    sd = [zeroSD,oneSD,twoSD,threeSD,fourSD,fiveSD];
end