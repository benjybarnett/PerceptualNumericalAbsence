function [zero_correct, one_correct, two_correct, three_correct,four_correct,five_correct] = GetCorrectRTs(data)

    %Calculates the average RT for each kind of match trial
    
    zero_correct = mean(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 0,32));
                
    one_correct = mean(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 1,32));
    
    two_correct = mean(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 2,32));

    three_correct = mean(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 3,32));

    four_correct = mean(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 4,32));
    
    five_correct = mean(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 5,32));

    
end