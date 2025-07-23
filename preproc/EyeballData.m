function EyeballData(cfg,subject)
    %Function to eyeball data 
    %Produces a topographic and temporal plot of the trail-averaged ERFs

    %load data
    data = load(fullfile(cfg.datadir,cfg.datatype,subject,[cfg.file_name,'.mat']));
    data = struct2cell(data);
    data = data{1};
    
    % Eyeball data
    tl_data = ft_timelockanalysis(cfg,data);
    figure; imagesc(tl_data.avg); %xticklabels(tl_data.time(1):0.2:tl_data.time(end)); xticks(linspace(1,size(tl_data.time,2),numel(xticklabels)));
    title(strcat(subject, ' ',cfg.file_name));
    hold on;
    figure; ft_topoplotER(cfg,tl_data); %title(strcat(subject, ' ',cfg.file_name));
    drawnow;

end