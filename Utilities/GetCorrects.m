function [zero_correct, one_correct, two_correct, three_correct,four_correct,five_correct] = GetCorrects(data)

    %Calculates the correct rate for match trials
    
    zero_correct = size(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 0,:),1)/...
                    size(data(data(:,28) == 0 & data(:,30) == 1),1) * 100;
                
    one_correct = size(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 1,:),1)/...
        size(data(data(:,28) == 1 & data(:,29) == 1),1) * 100;

    two_correct = size(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 2,:),1)/...
        size(data(data(:,28) == 2 & data(:,29) == 2),1) * 100;
       
    three_correct = size(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 3,:),1)/...
        size(data(data(:,28) == 3 & data(:,29) == 3),1) * 100;

    four_correct = size(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 4,:),1)/...
        size(data(data(:,28) == 4 & data(:,29) == 4),1) * 100;

    five_correct = size(data(data(:,30) == 1 & data(:,33) == 1 & data(:,28) == 5,:),1)/...
        size(data(data(:,28) == 5 & data(:,29) == 5),1) * 100;
end