function dataMain = AlignPDiode(cfg0,data,subject,block)

    cfg = [];
    cfg.channel = cfg0.trialdef.pdiodetype ;
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 0.5;
    lightDiodeSignal = ft_preprocessing(cfg, data);
   
    % determine the onset of the visual stimulus
    visOnset = [];
    for iTrial = 1:length(lightDiodeSignal.trial)
       
        PD_on = lightDiodeSignal.trial{iTrial} < (mean(lightDiodeSignal.trial{iTrial})-std(lightDiodeSignal.trial{iTrial})*1); %get when the PD is on (when is less than 1 SD below mean)
        PD_on_idx = find(PD_on);   
        visOnset(iTrial) = lightDiodeSignal.time{iTrial}(PD_on_idx(1));
       
    end
    
    %figure;
    %scatter(1:length(visOnset),visOnset)
    if cfg0.plot
        figure;
        for i = 1:length(lightDiodeSignal.time)
            plot(lightDiodeSignal.time{i},lightDiodeSignal.trial{i})
            %plot(data.time{i},data.trial{i}(314,:),'Color','cyan')
            hold on
            
            
        end
        xline(0,'r')
        title(sprintf('Unaligned Trials: %s, Block %d',subject,block))
    end
    
    % realign the trials to this onset
    cfg = [];
    cfg.offset = -visOnset * data.fsample;
    dataMain = ft_redefinetrial(cfg, data);

    %cut segments of interest
    cfg =[];
    cfg.toilim = [-cfg0.prestim cfg0.poststim];
    dataMain = ft_redefinetrial(cfg,dataMain);

    
    %plot New redefined trials
    cfg = [];
    cfg.channel = cfg0.trialdef.pdiodetype ;
    lightDiodeSignal = ft_preprocessing(cfg, dataMain);
    
    if cfg0.plot
        figure;
        for i =1:length(lightDiodeSignal.time)
            plot(lightDiodeSignal.time{i},lightDiodeSignal.trial{i})
            %plot(dataMain.time{i},dataMain.trial{i}(314,:),'Color','cyan')
            hold on
        end
        xline(0,'r')
        title(sprintf('Aligned Trials: %s, Block %d',subject,block))
        drawnow;
    end
    
    clear lightDiodeSignal visOnset
    
end