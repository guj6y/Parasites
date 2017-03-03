
web = 'BahiaFalsa';



model_files = cell(7,1);

model_files{1} = 'Model1';
model_files{2} = 'Model2';
model_files{3} = 'Model3';
model_files{4} = 'Model4';


web_names = cell(7,1);
web_names{1} = 'model1Links_';
web_names{2} = 'model2Links_';
web_names{3} = 'model3Links_';
web_names{4} = 'model4Links_';


for ii = 1:4
    webs_dir = sprintf('Generated_webs/%s/%s/%s'...
        ,web,model_files{ii},web_names{ii});
    for jj = 1:1000
        
        web_jj_dir = sprintf('%s%04u.csv',webs_dir,jj);
        
        [M] = csvread(web_jj_dir);
        
        res = M(:,1);
        cons = M(:,2);
        
        web_jj_dir_new = sprintf('Generated_webs/%s/%s/a%s%04u.txt'...
        ,web,model_files{ii},web_names{ii},jj);
        
        fid = fopen(web_jj_dir_new,'w');
        fprintf(fid,'%u,%u\n',[cons res]');
        fclose(fid);
        
        
    end
    
    
end