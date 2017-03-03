
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
web_names{4} = 'Flensburg'; %My methods fail here; too few links.
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
Trials = 1000;
num_webs = 0;

%11 models right now.
parms = [0 1 0 0; 
        0 1 0 1; 
        0 1 1 0; 
        0 1 1 1;
        0.1 0.3 0 0; 
        0.1 0.3 0 1; 
        0.1 0.3 1 0; 
        0.1 0.3 1 1];
for jj = 1  %Have more webs, but these are the only ones that it seems practical to do this for.
    
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
   
    
    Data_dir = sprintf('%s/%s',Data_dir0,web_name);
    try
        cd(Data_dir);
        cd(Code_dir);
    catch
        mkdir(Data_dir);
    end
    
    for ii = 0:8
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
        %fprintf(fid,'C,Cff,Cpf,Cfp,Cpp,svd,clustMN,gen1,gen2,gen3,gen4,gen5,vul1,vul2,vul3,vul4,vul5,link1,link2,link3,link4,link5,nsccs,floop,fbasal,fherb,fomn,ftop,fint,fcann,meansp,fsp,diam,meanswTL,stdswTL,skewswTL,maxswTL,mxsim\n');
        %              1, 2,  3,  4,  5,  6,  7,       8,    9,  10,  11,  12, 13,  14,  15,   16,   17,  18,  19,  20,    21,    22, 23,    24,   25,    26,    27, 28, 29, 30,    31,   32,  33,  34,       35,     36,      37      38      ,
        %          
        %              
        fclose(fid);
        
    end
    
    
    while num_webs < Trials
        num_webs = num_webs + 1;
        if saving(1)
            Model_dir = sprintf('%s/Model0',Data_dir);
            if resaving(1)
                LL = csvread(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs));
                res = LL(:,1);
                con = LL(:,2);
                
                props = calculateProperties(res,con);
                
                C_ = length(LL)/N^2;
                indxPara = datasample(1:N,Np);
                para = zeros(N,1)&false;
                para(indxPara) = true;
                
                free = ~para;
                
                nicheWeb = sparse(res,con,1,N,N);
                Cff_ = sum(sum(nicheWeb(free,free)))/Nf^2;
                Cpf_ = sum(sum(nicheWeb(para,free)))/Nf/Np;
                Cfp_ = sum(sum(nicheWeb(free,para)))/Nf/Np;
                Cpp_ = sum(sum(nicheWeb(para,para)))/Np^2;
                
                C_all = full([C_ Cff_ Cpf_ Cfp_ Cpp_]);
                
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f, \n ',[C_all,props]);
                fclose(fid);
                
                nicheSpace = csvread(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),1,0);
                nicheSpace = reshape(nicheSpace',[5,N])';
                nicheSpace(:,2) = para;
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%u,%u,%.5f,%.5f,%.5f\n',nicheSpace');
                fclose(fid);
                
            else
                %-----------------------Model 0
                [res, con, n, c, r,C_all,para] = NicheModelPara0(N,C,Np);
                props = calculateProperties(res,con);
                
                %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f, \n ',[C_all,props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%.5f\t%.5f\t%.5f\n',[(1:N)' para n c r ]');
                fclose(fid);
            end
        end
        
        if saving(2)
            Model_dir = sprintf('%s/Model1',Data_dir);
            if resaving(2)
                
                LL = csvread(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs));
                res = LL(:,1);
                con = LL(:,2);
                
                props = calculateProperties(res,con);
                
                C_ = length(LL)/N^2;
                indxPara = datasample(1:N,Np);
                para = zeros(N,1)&false;
                para(indxPara) = true;
                
                free = ~para;
                
                nicheWeb = sparse(res,con,1,N,N);
                Cff_ = sum(sum(nicheWeb(free,free)))/Nf^2;
                Cpf_ = sum(sum(nicheWeb(para,free)))/Nf/Np;
                Cfp_ = sum(sum(nicheWeb(free,para)))/Nf/Np;
                Cpp_ = sum(sum(nicheWeb(para,para)))/Np^2;
                
                C_all = full([C_ Cff_ Cpf_ Cfp_ Cpp_]);
                
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f, \n ',[C_all,props]);
                fclose(fid);
                
            else
                %-----------------------Model 1
                [~, ~, n0, c0, r0,~] = NicheModelPara0(Nf,C,Np);
                [res, con, n, c, r,C_all] = ...
                    AddParasites1(n0,c0,r0,Np,[C,Cfp],0);
                props = calculateProperties(res,con);
                
                %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[C_all,props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%u,%.5f\t%.5f\t%.5f\n',[(1:N)' n c r ]');
                fclose(fid);
            end
        end
        
        if saving(3)
            Model_dir = sprintf('%s/Model2',Data_dir);
            if resaving(3)
                LL = csvread(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs));
                res = LL(:,1);
                con = LL(:,2);
                
                props = calculateProperties(res,con);
                
                C_ = length(LL)/N^2;
                indxPara = datasample(1:N,Np);
                para = zeros(N,1)&false;
                para(indxPara) = true;
                
                free = ~para;
                
                nicheWeb = sparse(res,con,1,N,N);
                Cff_ = sum(sum(nicheWeb(free,free)))/Nf^2;
                Cpf_ = sum(sum(nicheWeb(para,free)))/Nf/Np;
                Cfp_ = sum(sum(nicheWeb(free,para)))/Nf/Np;
                Cpp_ = sum(sum(nicheWeb(para,para)))/Np^2;
                
                C_all = full([C_ Cff_ Cpf_ Cfp_ Cpp_]);
                
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f, \n ',[C_all,props]);
                fclose(fid);
                
            else
                
                
                %-----------------------Model 2
                [~,~,n0,c0,r0,~]= NicheModelPara(Nf, C, Np);
                [res, con, n, c, r, C_all] = ...
                    AddParasites1(n0,c0,r0,Np,[C,Cfp],1);
                props = calculateProperties(res,con);
                
                %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[C_all,props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%u\t%u\t%.5f\t%.5f\t%.5f\n',[(1:N)' [zeros(Nf,1);ones(Np,1)] n c r ]');
                fclose(fid);
            end
            
        end
        
        if saving(4)
            Model_dir = sprintf('%s/Model3',Data_dir);
            if resaving(4)
                LL = csvread(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs));
                res = LL(:,1);
                con = LL(:,2);
                
                props = calculateProperties(res,con);
                
                C_ = length(LL)/N^2;
                indxPara = datasample(1:N,Np);
                para = zeros(N,1)&false;
                para(indxPara) = true;
                
                free = ~para;
                
                nicheWeb = sparse(res,con,1,N,N);
                Cff_ = sum(sum(nicheWeb(free,free)))/Nf^2;
                Cpf_ = sum(sum(nicheWeb(para,free)))/Nf/Np;
                Cfp_ = sum(sum(nicheWeb(free,para)))/Nf/Np;
                Cpp_ = sum(sum(nicheWeb(para,para)))/Np^2;
                
                C_all = full([C_ Cff_ Cpf_ Cfp_ Cpp_]);
                
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f, \n ',[C_all,props]);
                fclose(fid);
                
            else
                
            end
            %-----------------------Model 3
            np_min = 0;
            [~,~,n0,c0,r0,~]= NicheModelPara(Nf, C, Np);
            [res, con, n, c, r,C_all] = ...
                AddParasites2(n0,c0,r0,Np,[C,Cfp],np_min,np_max,0);
            props = calculateProperties(res,con);
            
            %-------------------%%%Saving output%%-------------------%
            fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
            fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[C_all,props]);
            fclose(fid);
            
            fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
            fprintf(fid,'%u,%u\n',[res con]');
            fclose(fid);
            
            fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
            fprintf(fid,'indx\tpara\tn\tc\tr\n');
            fprintf(fid,'%.5f\t%.5f\t%.5f\n',[(1:N)' [zeros(Nf,1);ones(Np,1)] n c r ]');
            fclose(fid);
        end
        
        if saving(5)
            Model_dir = sprintf('%s/Model4',Data_dir);
            if resaving(5)
                LL = csvread(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs));
                res = LL(:,1);
                con = LL(:,2);
                
                props = calculateProperties(res,con);
                
                C_ = length(LL)/N^2;
                indxPara = datasample(1:N,Np);
                para = zeros(N,1)&false;
                para(indxPara) = true;
                
                free = ~para;
                
                nicheWeb = sparse(res,con,1,N,N);
                Cff_ = sum(sum(nicheWeb(free,free)))/Nf^2;
                Cpf_ = sum(sum(nicheWeb(para,free)))/Nf/Np;
                Cfp_ = sum(sum(nicheWeb(free,para)))/Nf/Np;
                Cpp_ = sum(sum(nicheWeb(para,para)))/Np^2;
                
                C_all = full([C_ Cff_ Cpf_ Cfp_ Cpp_]);
                
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f, \n ',[C_all,props]);
                fclose(fid);
                
            else
                %-----------------------Model 4
                [~,~,n0,c0,r0,~]= NicheModelPara(Nf, C, Np);
                [res, con, n, c, r, C_all] = ...
                    AddParasites2(n0,c0,r0,Np,[C,Cfp],np_min,np_max,1);
                props = calculateProperties(res,con);
                
                %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[C_all,props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%.5f\t%.5f\t%.5f\n',[(1:N)' [zeros(Nf,1);ones(Np,1)] n c r ]');
                fclose(fid);
            end
        end
        if saving(6)
            Model_dir = sprintf('%s/Model5',Data_dir);
            if resaving(6)
                LL = csvread(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs));
                res = LL(:,1);
                con = LL(:,2);
                
                props = calculateProperties(res,con);
                
                C_ = length(LL)/N^2;
                indxPara = datasample(1:N,Np);
                para = zeros(N,1)&false;
                para(indxPara) = true;
                
                free = ~para;
                
                nicheWeb = sparse(res,con,1,N,N);
                Cff_ = sum(sum(nicheWeb(free,free)))/Nf^2;
                Cpf_ = sum(sum(nicheWeb(para,free)))/Nf/Np;
                Cfp_ = sum(sum(nicheWeb(free,para)))/Nf/Np;
                Cpp_ = sum(sum(nicheWeb(para,para)))/Np^2;
                
                C_all = full([C_ Cff_ Cpf_ Cfp_ Cpp_]);
                
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f, \n ',[C_all,props]);
                fclose(fid);
                
            else
                %-----------------------Model 5
                np_min = 0.2;
                np_max = 0.5;
                [~,~,n0,c0,r0,~,]= NicheModelPara(Nf, C, Np);
                [res, con, n, c, r] = ...
                    AddParasites2(n0,c0,r0,Np,[C,Cfp],np_min,np_max,0);
                props = calculateProperties(res,con);
                
                %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[C_all,props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%.5f\t%.5f\t%.5f\n',[(1:N)' [zeros(Nf,1);ones(Np,1)] n c r ]');
                fclose(fid);
            end
        end
        
        if saving(7)
            Model_dir = sprintf('%s/Model6',Data_dir);
            if resaving(7)
                LL = csvread(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs));
                res = LL(:,1);
                con = LL(:,2);
                
                props = calculateProperties(res,con);
                
                C_ = length(LL)/N^2;
                indxPara = datasample(1:N,Np);
                para = zeros(N,1)&false;
                para(indxPara) = true;
                
                free = ~para;
                
                nicheWeb = sparse(res,con,1,N,N);
                Cff_ = sum(sum(nicheWeb(free,free)))/Nf^2;
                Cpf_ = sum(sum(nicheWeb(para,free)))/Nf/Np;
                Cfp_ = sum(sum(nicheWeb(free,para)))/Nf/Np;
                Cpp_ = sum(sum(nicheWeb(para,para)))/Np^2;
                
                C_all = full([C_ Cff_ Cpf_ Cfp_ Cpp_]);
                
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f, \n ',[C_all,props]);
                fclose(fid);
                
            else
                %-----------------------Model 6
                [~,~,n0,c0,r0,~]= NicheModelPara(Nf, C, Np);
                [res, con, n, c, r, C_all] = ...
                    AddParasites2(n0,c0,r0,Np,[C,Cfp],np_min,np_max,1);
                props = calculateProperties(res,con);
                
                %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[C_all,props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%.5f\t%.5f\t%.5f\n',[(1:N)' [zeros(Nf,1);ones(Np,1)] n c r ]');
                fclose(fid);
            end
        end
        if saving(8)
            Model_dir = sprintf('%s/Model7',Data_dir);
            if resaving(8)
                LL = csvread(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs));
                res = LL(:,1);
                con = LL(:,2);
                
                props = calculateProperties(res,con);
                
                C_ = length(LL)/N^2;
                indxPara = datasample(1:N,Np);
                para = zeros(N,1)&false;
                para(indxPara) = true;
                
                free = ~para;
                
                nicheWeb = sparse(res,con,1,N,N);
                Cff_ = sum(sum(nicheWeb(free,free)))/Nf^2;
                Cpf_ = sum(sum(nicheWeb(para,free)))/Nf/Np;
                Cfp_ = sum(sum(nicheWeb(free,para)))/Nf/Np;
                Cpp_ = sum(sum(nicheWeb(para,para)))/Np^2;
                
                C_all = full([C_ Cff_ Cpf_ Cfp_ Cpp_]);
                
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f, \n ',[C_all,props]);
                fclose(fid);
                
            else
                %-----------------------Model 7
                [res,con,n,c,r, C_all] = two_axisNicheModel(Nf,Np,[Cff Cpf Cfp Cpp],0);
                props = calculateProperties(res,con);
                
                %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[C_all',props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%.5f\t%.5f\t%.5f\n',[(1:N)' [zeros(Nf,1);ones(Np,1)] n c r ]');
                fclose(fid);
            end
        end
        if saving(9)
            Model_dir = sprintf('%s/Model8',Data_dir);
            if resaving(9)
                LL = csvread(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs));
                res = LL(:,1);
                con = LL(:,2);
                
                props = calculateProperties(res,con);
                
                C_ = length(LL)/N^2;
                indxPara = datasample(1:N,Np);
                para = zeros(N,1)&false;
                para(indxPara) = true;
                
                free = ~para;
                
                nicheWeb = sparse(res,con,1,N,N);
                Cff_ = sum(sum(nicheWeb(free,free)))/Nf^2;
                Cpf_ = sum(sum(nicheWeb(para,free)))/Nf/Np;
                Cfp_ = sum(sum(nicheWeb(free,para)))/Nf/Np;
                Cpp_ = sum(sum(nicheWeb(para,para)))/Np^2;
                
                C_all = full([C_ Cff_ Cpf_ Cfp_ Cpp_]);
                
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f, \n ',[C_all,props]);
                fclose(fid);
                
            else
                %-----------------------Model 8
                [res,con,n,c,r, C_all] = two_axisNicheModel(Nf,Np,[Cff Cpf Cfp Cpp],1);
                props = calculateProperties(res,con);
                
                %-------------------%%%Saving output%%-------------------%
                fid = fopen(sprintf('%s/properties.csv',Model_dir),'a');
                fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[C_all',props]);
                fclose(fid);
                
                fid = fopen(sprintf('%s/webs/web_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'%u,%u\n',[res con]');
                fclose(fid);
                
                fid = fopen(sprintf('%s/niches/nichespace_%04u.csv',Model_dir,num_webs),'w');
                fprintf(fid,'indx\tpara\tn\tc\tr\n');
                fprintf(fid,'%.5f\t%.5f\t%.5f\n',[(1:N)' [zeros(Nf,1);ones(Np,1)] n c r ]');
                fclose(fid);
            end
        end
        clc
        percents01 = num_webs/10;
        fprintf('%s is %.1f%% complete\n',web_name,percents01)
        
        
        
        
    end
    num_webs = 0;
end

%}