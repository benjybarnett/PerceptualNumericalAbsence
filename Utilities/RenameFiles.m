
root = '/vol/ccnlab1/naddij/UPCPIM/P02';

oldPart = 'P01'; newPart = 'P02';
files = str2fullfile(root,['*' oldPart '*.IMA']);

% loop through the files and rename them
for file = 2:length(files)
    
    if mod(file,50) == 0 
        fprintf('Renaming file %d out of %d \n',file,length(files))
    end
    
   oldname = files{file}; 
   idx = strfind(oldname,oldPart);
   idx = idx:idx+length(oldPart)-1;
   
   newname = oldname; 
   newname(idx) = newPart;
    
   movefile(oldname,newname); 
    
end



