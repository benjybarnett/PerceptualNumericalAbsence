function data = appendTasks(cfg0,subject)
    
    exp_data = load(fullfile('D:\bbarnett\Documents\Zero\data\VARData\',subject,'data_VAR.mat'));
    loc_data = load(fullfile('D:\bbarnett\Documents\Zero\data\VARData\',subject,'loc_data_VAR.mat'));
    
    exp_data = exp_data.data;
    dataS{1} = exp_data;

    %make trial info arrays matching dimensions
    loc_data = loc_data.data;
    full_loc_trial_info = zeros(size(loc_data.trialinfo,1),size(exp_data.trialinfo,2));
    full_loc_trial_info(:,1) = loc_data.trialinfo;
    loc_data.trialinfo = full_loc_trial_info;

    dataS{2} = loc_data;

    %append
    cfgA = []; 
    cfgA.keepsampleinfo='yes';
    data = ft_appenddata(cfgA,dataS{:});


    % save and clean up
    save(fullfile(cfg0.saveDir,subject,cfg0.saveName),'data','-v7.3')
    clear cfg data
        

    
end