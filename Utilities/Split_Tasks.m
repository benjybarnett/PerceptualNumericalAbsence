function Split_Tasks(cfg0,subject)
    

    datadir = cfg0.datadir;
    load(fullfile(datadir,'CleanData',subject,'data.mat'));  

  
    %% Save By Task

    %Arabic task
    root = 'D:\bbarnett\Documents\Zero\data\CleanData';
    %cfg.channel = 'MEG';
    cfg.trials = data.trialinfo(:,3) == 2; %select detection task
    arabic_data = ft_selectdata(cfg,data);
    save(fullfile(root,subject,'arabic_data.mat'),'arabic_data');
    
    %numerical task
    cfg = [];
    root = 'D:\bbarnett\Documents\Zero\data\CleanData';
    %cfg.channel = 'MEG';
    cfg.trials = data.trialinfo(:,3) == 1;
    dot_data = ft_selectdata(cfg,data);
    save(fullfile(root,subject,'dot_data.mat'),'dot_data');




    clear data dot_data arabic_data 


end