function PreprocessingICA_BOB(cfg0,subject)
% function PreprocessingICA_BOB(subjectdata)

% some settings
outputComp              = fullfile(cfg0.datadir,'ICAData',subject);
outputData             = fullfile(cfg0.datadir,'CleanData',subject);
if ~exist(outputComp,'dir'); mkdir(outputComp); end
if ~exist(outputData,'dir'); mkdir(outputData); end
VARData             =  fullfile(cfg0.datadir,'VARData',subject); 

%% ICA

% load the VA removed data
data = load(fullfile(VARData,cfg0.inputData));
data = struct2cell(data);data = data{1};

comp_output = strcat(cfg0.compOutput,'.mat');

% check if ICA already done
if ~exist(fullfile(outputComp,comp_output),'file')
    
    % perform the independent component analysis (i.e., decompose the data)
    cfg                 = [];
    cfg.channel         = 'MEG';
    cfg.method          = 'runica';
    cfg.demean          = 'no';
    comp                = ft_componentanalysis(cfg,data);
    
    % save the components
    save(fullfile(outputComp,comp_output),'comp','-v7.3')
    
else 
   load(fullfile(outputComp,comp_output))
end

% identify EOG and ECG components    

% correlate to eye tracking
ET = cell2mat(data.trial);
ET = ET(ismember(data.label, {'UADC001','UADC002'}), :);

disp(size(cell2mat(comp.trial)'))
disp(size(ET'))

r = corr(cell2mat(comp.trial)', ET');
[ro, i] = sort(abs(r),'descend');

fprintf('Highest correlations: \n \t X: comp %d [%.4f] \n \t Y: comp %d [%.4f] \n \t ',i(1,1),r(i(1,1),1),i(1,2),r(i(1,2),2))

if ~isempty(find(ro(1,:)<0.3,1)) % if some are below 0.3
    % manually check the components
    warning('low correlations, manually check components!')
    %return
end

% inspect these components
figure;
cfg                = [];
cfg.component      = i(1:6,:);
cfg.layout         = 'CTF275.lay';
cfg.commment       = 'no';
ft_topoplotIC(cfg,comp)




% plot the time course
tmp_comp = cell2mat(comp.trial);
figure;
flat_idx = cfg.component(:);
nComps = length(flat_idx);
for c = 1:nComps/2
    subplot(nComps/2,1,c);
    plot(tmp_comp(flat_idx(c),1:2000))
    title(sprintf('Component %d \n',flat_idx(c)))
end
drawnow

%print top 5 correlations of components with EOG  
for j = 1:5
    disp(sprintf('The %d th top correlations are %.2f and %.2f, for components %d and %d', j, ro(j,1),ro(j,2),i(j,1),i(j,2)))
end

% decide which components to remove and save decision
remove_comp = input('Enter the components you would like to remove in the format [a b]');
comp_removed        = remove_comp;


%view plots of all components to find heart artefacts
%time course of the independent components
n_total_comp = length(ro);
pos = 1;
figure('Units','normalized','Position',[0 0 1 1])
for c = 1:n_total_comp
    if mod(c,40) == 0
        figure('Units','normalized','Position',[0 0 1 1])
        pos=1;
    end
    subplot(20,2,pos)
    plot(tmp_comp(c,1:2000))
    xlabel("   "+newline+"   ")
    title(sprintf('Component %d \n',c))
    pos= pos+1;
end
drawnow

ecg_comp = input('Enter the component(s) with ECG artefacts');
%plot topography to check
figure('Units','normalized','Position',[0 0 1 1])
cfg = [];
cfg.component = ecg_comp;
disp(ecg_comp)
cfg.layout = 'CTF275.lay';
cfg.commment       = 'no';
ft_topoplotIC(cfg,comp)
drawnow

remove = input('Which of these components would you like to remove? Press enter if none');
comp_removed = [comp_removed remove];
disp('Removing the following components:')
disp(comp_removed);

% remove them from the data
cfg = [];
cfg.component       = comp_removed;
cfg.demean ='no';
data             = ft_rejectcomponent(cfg, comp, data);
close all
save(fullfile(outputComp,comp_output),'comp_removed','-append')



save(fullfile(outputData,strcat(cfg0.outputData,'.mat')),'data','-v7.3'); clear data


clear comp comp_removed


clear data

end


