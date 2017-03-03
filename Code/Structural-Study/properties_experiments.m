
%compare the 7 levels of niche model with parasites to the 7 (4) empirical
%webs.

%Four levels of niche model:
%Model 0: Increase S in ONM
%Model 1: Parastes eat anywhere, match 2 connectances (Cff,C)
%Model 2: Parasites eat above,  match 2 connectances (Cff,C)
%Model 3: Parasites eat anywhere, are restricted, match 2 connectances
%(Cff,C)
%Model 4: Parasites eat above, are restricted, match 2 connectances (Cff,C)
%Model 5: Parasites eat anywhere, are restricted, match 2 connectances
%(Cff,C)
%Model 6: Parasites eat above, are restricted, match 2 connectances (Cff,C)
%Model 7: Parasites eat anywhere, match all connectances (not restricted)
%Model 8: Parasites eat above, match all connectances (not restricted)
%%%*4b,5b have different ranges.

%Empirical Webs: (one of these is bad.)
%Bahia Falsa
%Carpinteria
%Flensberg Fjord
%Ythan Estuary (this one is bad.)
%Sylt Tidal Basin
%Otago Harbor
%Punta de la Banda (something like that)

%Try all of these for now.
web_names = cell(6,1);
web_names{5} = 'BahiaFalsa';
web_names{6} = 'Carpinteria';
web_names{4} = 'Flensburg'; %Doesn't work even with new stuff.
web_names{2} = 'Otago';
web_names{3} = 'Punta'; %My methods fail here
web_names{1} = 'Sylt';
%I don't trust ythan since it doesn't include p->f links or p->p links.
%web_name = 'Ythan';
np_max = 0.3;
%Where to save data: change to match your system
Data_dir0 = '/Volumes/Buster/Research/Parasitism/InverseNicheModels';
Code_dir = pwd;
try
    cd(Data_dir0);
    cd(Code_dir);
catch
    error(sprintf('Buster is neither seen nor heard..\n(Can''t cd to data directory)'))
end
Trials = 100;
num_webs = 0;

%11 models right now.
%[np_min np_max para_above match_fp
n_models = 6;
parms = [0 1 0 0; 
        0 1 0 1; 
        0 1 1 0; 
        0 1 1 1;
        %0.7 0.9 0 0; %These still aren't doing what I want; can't match
        %overall without knowing the expected overlap.
        0.7 0.9 0 1; 
        %0.1 0.3 1 0; %These still aren't doing what I want.  ditto above
        0.1 0.3 1 1];
for jj = 4  %Have more webs, but these are the only ones that it seems practical to do this for.
    
    %TODO: organize directory and make what paths I can relative.
    %First, load the empirical web
    web_name = web_names{jj};
    web_dir = sprintf('/Users/Nick/Documents/Research/FoodWebData/Matlab_web_data2/%s.mat',web_name);
    
    web = load(web_dir);
    web=web.web;
    
    Np = web.para.taxa.metrics.Np;
    N = web.para.taxa.metrics.N;
    Nf = N-Np;
    
    Cff = web.para.taxa.metrics.Cff;
    Cpf = web.para.taxa.metrics.Cpf;
    Cfp = web.para.taxa.metrics.Cfp;
    Cpp = web.para.taxa.metrics.Cpp;
    C = web.para.taxa.metrics.C;
    C5 = [C,Cff,Cpf,Cfp,Cpp];
    
    Data_dir = sprintf('%s/%s',Data_dir0,web_name);
    try
        cd(Data_dir);
        cd(Code_dir);
    catch
        mkdir(Data_dir);
    end
    
    for ii = 0:11
        Model_dir = sprintf('%s/Model%u',Data_dir,ii);
        try
            cd(Model_dir)
            cd(Code_dir)
        catch
            mkdir Model_dir
            mkdir(sprintf('%s/webs',Model_dir))
            mkdir(sprintf('%s/niches',Model_dir))
        end
        fid = fopen(sprintf('%s/properties.csv',Model_dir),'w');
        fprintf(fid,'C,Cff,Cpf,Cfp,Cpp,svd,clustMN,gen1,gen2,gen3,gen4,gen5,vul1,vul2,vul3,vul4,vul5,link1,link2,link3,link4,link5,nsccs,floop,fbasal,fherb,fomn,ftop,fint,fcann,meansp,fsp,diam,meanswTL,stdswTL,skewswTL,maxswTL,mxsim\n');
        %              1, 2,  3,  4,  5,  6,  7,       8,    9,  10,  11,  12, 13,  14,  15,   16,   17,  18,  19,  20,    21,    22, 23,    24,   25,    26,    27, 28, 29, 30,    31,   32,  33,  34,       35,     36,      37      38      ,
        %          
        %              
        fclose(fid);
        
    end
    
    
    while num_webs < Trials
        num_webs = num_webs + 1;

            
           

                %-----------------------Model 0
                Model_dir = sprintf('%s/Model0',Data_dir);
                [res, con, n, c, r] = NicheModel_nk(N,C);
                props = calculateProperties(res,con);
                
                C_web = length(res)/N^2;
                %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n ',[C_web,C_web,C_web,C_web,C_web,props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx,n,c,r\n');
                fprintf(fid,'%u,%.5f,%.5f,%.5f\n',[(1:N)' n c r ]');
                fclose(fid);

                

        fprintf('Model 1\n')   

                %-----------------------Model 1
                para_above = 0;
                Model_dir = sprintf('%s/Model1',Data_dir);
                [res,con,n,c,r,C_all,para]= NicheModelPara_test(N, C, Np, para_above);
                props = calculateProperties(res,con);
                
                %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n ',[C_all,props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%u,%u,%.9f,%.9f,%.9f\n',[(1:N)' para n c r ]');
                fclose(fid);
                

        fprintf('model2\n')

                %-----------------------Model 2
                para_above = 1;
                Model_dir = sprintf('%s/Model2',Data_dir);
                [res,con,n,c,r,C_all,parasites]= NicheModelPara_test(N, C, Np, para_above);
                props = calculateProperties(res,con);
                
                 %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n ',[C_all,props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%u,%u,%.9f,%.9f,%.9f\n',[(1:N)' para n c r ]');
                fclose(fid);
                
        for kk = 1:n_models
            fprintf('model%u\n',kk+2)
            model_num = 2+kk;
            Model_dir = sprintf('%s/Model%u',Data_dir,model_num);
            %-----------------------Model kk
            [~,~,n_new,c_new,r_new]= NicheModel_nk(Nf, C) ;
            [Res,Cons,n,c,r,C_all,para] = ...
                AddParasites(n_new,c_new,r_new,Np,C5,...
                parms(kk,1),parms(kk,2),parms(kk,3),parms(kk,4));
            props = calculateProperties(res,con);
                
                 %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n ',[C_all,props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%u,%u,%.9f,%.9f,%.9f\n',[(1:N)' para n c r ]');
                fclose(fid);

        end
                
        clc
        percents01 = num_webs;
        fprintf('%s is %u%% complete\n',web_name,percents01)

    end
    num_webs = 0;
end

%}