close all
clear all
clc
web_names = cell(6,1);
web_names{5} = 'BahiaFalsa';
web_names{6} = 'Carpinteria';
web_names{4} = 'Flensburg';
web_names{2} = 'Otago';
web_names{3} = 'Punta';
web_names{1} = 'Sylt';
%I don't trust ythan since it doesn't include p->f links or p->p links.
%web_name = 'Ythan';
np_max = 0.3;
%Director

Nf = 1;
Np = 1;
for jj =1:6
    %First, load the empirical web
    web_name = web_names{jj};
    web_dir = sprintf('/Users/Nick/Documents/Research/FoodWebData/Matlab_web_data2/%s.mat',web_name);
    
    web = load(web_dir);
    web=web.web;
    %
    Np(jj) = web.para.taxa.metrics.Np;
    N = web.para.taxa.metrics.N;
    Nf(jj) = N-Np(jj);
    
    Cff(jj) = web.para.taxa.metrics.Cff;
    Cpf(jj) = web.para.taxa.metrics.Cpf;
    Cfp(jj) = web.para.taxa.metrics.Cfp;
    Cpp(jj) = web.para.taxa.metrics.Cpp;
    C(jj) = web.para.taxa.metrics.C;
    
end
figure
bar(Cff)
figure
bar(Cpf)
figure
bar(Cfp)
figure
bar(Cpp)
figure
bar(C)