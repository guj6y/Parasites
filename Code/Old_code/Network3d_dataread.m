close all
clear all
clc

file_name = cell(6,1);
file_name{1} = 'Carp_model1'; %niche iwth increased S
file_name{2} = 'Carp_model2'; %niche with parasites reverse
file_name{3} = 'Carp_model3'; %inverse niche parasites reverse 0-.3
file_name{4} = 'Carp_model3b';%inverse niche parasites reverse .1-.3
file_name{5} = 'Carp_model4'; %inverse niche parasites reverse all C 0-.3
file_name{6} = 'Carp_model4b';%inverse niche parasites reverse all C .1-.3

results=cell(6,1);

averages = cell(6,1);
devs = cell(6,1);

actual = zeros(1,19);

for ii = 1:6
    results{ii} = zeros(1000,19);
     fid = fopen(sprintf('Data_good/%s.txt',file_name{ii}),'r');
    A = textscan(fid,'%s %s\n',1);
    B = textscan(fid,'%s',20);
    C = textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);

    fid = fopen('carp_emp_anal.txt','r');
    textscan(fid,'%s %s\n',1);
    textscan(fid,'%s',20);
    emp = textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);

    for jj  = 2:20
        C{jj}(1:2:2000) = [];
        results{ii}(:,jj-1) = C{jj}(:);
        emp{jj}(1) = [];
        actual(jj-1) = emp{jj};
    end
    
    averages{ii} = mean(results{ii});
    devs{ii} = std(results{ii});
    
    
end

err = cell(6,1);

for ii = 1:6
    err{ii} = (averages{ii} - actual)./devs{ii};
    
end

figure

colors = ['b.';'g.';'y.';'c.';'k.';'m.'];

for jj =1:6
    hold on
    plot(1:6,err{jj}(2:7),colors(jj,:))
end
legend('plain','inverse','inverse restrict1','inverse restrict0','inverse all1','inverse all0')
%axis([1 6 -5 5])
title(sprintf('%s,%s,%s,%s,%s,%s',B{1}{3},B{1}{4},B{1}{5},B{1}{6},B{1}{7},B{1}{8}))
%
figure
for jj = 1:6
  hold on
  plot(1:6,err{jj}(8:13),colors(jj,:))  
    
end
%axis([1 6 -5 5])
legend('plain','inverse','inverse restrict1','inverse restrict0','inverse all1','inverse all0')
title(sprintf('%s,%s,%s,%s,%s,%s',B{1}{9},B{1}{10},B{1}{11},B{1}{12},B{1}{13},B{1}{14}))
figure
for jj = 1:6
    
    hold on
    plot(1:6,err{jj}(14:19),colors(jj,:))
end
%axis([1 6 -5 5])
legend('plain','inverse','inverse restrict1','inverse restrict0','inverse all1','inverse all0')
title(sprintf('%s,%s,%s,%s,%s,%s',B{1}{15},B{1}{16},B{1}{17},B{1}{18},B{1}{19},B{1}{20}))
%}


