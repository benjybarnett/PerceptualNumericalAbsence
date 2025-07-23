function appendLocalizer(cfg0,subject)
    
    VARData             =  fullfile(cfg0.datadir,'VARData',subject);
    
    % load the VA removed data
    data = load(fullfile(VARData,cfg0.inputData));
    data = data.data;
    
end