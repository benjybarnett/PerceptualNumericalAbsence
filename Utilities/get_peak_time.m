function [peak_sample,peak_time] = get_peak_time(cfg0,data)
    %Returns the sample and time where peak decoding occurs.
    %Input: [subj x timepoints] matrix of accuracy values
    %Output: peak sample and peak time
    group_acc = mean(data,1);
    [peak_acc,peak_sample] = max(group_acc);
    time = cfg0.time;
    peak_time = time(peak_sample);
   
    fprintf('The peak group accuracy was %0.3f at %0.3f ms after stimulus presentation for %s task \n',peak_acc,peak_time, cfg0.name)
    if cfg0.plot
        figure;
        plot(time,group_acc,'LineWidth',1); 
        hold on; 
        plot(xlim,[0.5 0.5],'k--','LineWidth',2)
        ylim([0.3,0.7])
        xlabel('Time (s)'); ylabel('Accuracy');
        xlim([time(1) time(end)]); 
        title(strcat(['Mean Decoding Accuracy For All ',' ',cfg0.name,' ', 'Classifiers']));
        xline(peak_time,'--r')
        hold off;
    end
end