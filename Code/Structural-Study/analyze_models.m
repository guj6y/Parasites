% This script analyzes the models made in properties_experiments.
close all
web_names = cell(4,1);
web_names{1} = 'BahiaFalsa';
web_names{2} = 'Carpinteria';
web_names{3} = 'Otago';
web_names{4} = 'Sylt';
web_names{5} = 'Punta';
property_names = {'C','Cf','Cp','Cfp','Cpp','svd','clustMN',...
    'gen1','gen2','gen3','gen4','gen5',...
    'vul1','vul2','vul3','vul4','vul5',...
    'link1','link2','link3','link4','link5',...
    'nsccs','floop','fbasal','fherb','fomn','ftop',...
    'fint','fcann','meansp','fsp','diam','meanswTL',...
    'stdswTL','skewswTL','maxswTL','mxsim'};

Data_dir0 = '/Volumes/Buster/Research/Parasitism/InverseNicheModels';
Code_dir = pwd;
%property_names([9:12,14:17,19:22]) = [];
models = 0:8;

num_hit = zeros(4,7);
normalized_error = zeros(4,7);
bin_hits = zeros(length(property_names),9,6);
for ii = 1:5
    
    %First, load the empirical web
    web_name = web_names{ii};
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
    
    rese = web.para.taxa.links.Resources;
    cone = web.para.taxa.links.Consumers;
    
    nicheweb = sparse(rese,cone,1,N,N);
    
    
    Data_dir = sprintf('%s/%s',Data_dir0,web_name);
    look_at = [1,2,3,24,32,34];
    emp_props = calculateProperties(rese,cone);
    emp_propsC = [C Cff Cpf Cfp Cpp emp_props];
    
    emp_propsC(2) = (emp_propsC(2)*Nf*Nf + emp_propsC(3)*Np*Nf)...
        /((Nf-(emp_propsC(25)*(Nf+Np)))*(Nf+Np));
    emp_propsC(3) = (emp_propsC(4)*Np*Nf + emp_propsC(5)*Np^2)...
        /(Np*(Nf+Np));
    
    
    
    %emp_propsC = emp_propsC(:,look_at);
    %emp_propsC([9:12,14:17,19:22]) = [];
    
    
    %for each model
    for jj = 0:8
        %fprintf(fid,'C,Cff,Cpf,Cfp,Cpp,svd,clustMN,gen1,gen2,gen3,gen4,gen5,vul1,vul2,vul3,vul4,vul5,link1,link2,link3,link4,link5,nsccs,floop,fbasal,fherb,fomn,ftop,fint,fcann,meansp,fsp,diam,meanswTL,stdswTL,skewswTL,maxswTL,mxsim\n');
        %             1, 2,  3,  4,  5,  6,    7,    8,    9,  10,  11,  12, 13,  14,  15,   16,   17,  18,  19,  20,    21,    22, 23,    24,   25,    26,   27, 28, 29,   30,    31,   32,  33,  34,       35,     36,      37      38      ,
        prop_file = sprintf('%s/Model%u/properties.csv',Data_dir,jj);
        props = csvread(prop_file,1,0);
        %props = props(:,look_at);
        props(:,2) = (props(:,2)*Nf*Nf + props(:,3)*Nf*Np)/((Nf-(emp_propsC(25)*(Nf+Np)))*(Nf+Np));
        props(:,3) = (props(:,4)*Np*Nf + props(:,5)*Np^2)/(Np*(Nf+Np));
        
        %props(:,[9:12,14:17,19:22]) = [];
        mean_props = mean(props);
        std_props = std(props);
        
        normalized_error(ii,jj+1) = mean((mean_props-emp_propsC)./std_props);
        hits = (mean_props + 2*std_props)>emp_propsC;
        hits = hits&((mean_props -2*std_props)<emp_propsC);
        bin_hits(:,jj+1,ii) = hits;
        
        num_hit(ii,jj+1) = sum(hits);
        for kk = 1:length(look_at)
            figure(kk)
            subplot(2,3,ii)
            hold on
            plot(-1:9,ones(1,11)*emp_propsC(look_at(kk)));
            plot(jj,mean_props(look_at(kk)),'k.')
            errorbar(jj,mean_props(look_at(kk)),2*std_props(kk),'k-')
            title(web_names{ii})
            xlabel('Model')
            ylabel(property_names{look_at(kk)})
            xlim([-0.5,8.5]);
            
        end
        
        %}
        
    end
    for kk = 1:length(look_at)
        h = figure(kk);
        print(h,'-djpeg',sprintf('../../Presentation-4-11-16/propPlot%u.jpg',kk));
        
    end
    
    
    
end