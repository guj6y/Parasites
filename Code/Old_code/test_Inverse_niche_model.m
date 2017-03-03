close all
clear all

web = load('Matlab_web_data2/Carpinteria.mat','web');

web = web.web;

%Determine number of each type of link:
%free liver feeds on free liver
%free liver feds on parasite
%Parasite feeds on free liver
%parasite feeds on parasite


Res_para = web.para.taxa.links.Resources;
Con_para = web.para.taxa.links.Consumers;

tic
a = walk1(1,[],Res_para,Con_para);
b = walk1(2,a,Res_para,Con_para);
c = walk1(3,b,Res_para,Con_para);
toc

LLPara = [Res_para,Con_para];

Types_para = web.para.taxa.names.Types;

LL_Types_para = Types_para(LLPara);

LL_ff = LLPara(LL_Types_para(:,1)==1,:);
LL_Types_para_ = LL_Types_para(LL_Types_para(:,1)==1,:);
LL_ff = LL_ff(LL_Types_para_(:,2)==1,:);

LL_fp = LLPara(LL_Types_para(:,1)==1,:);
LL_Types_para_ = LL_Types_para(LL_Types_para(:,1)==1,:);
LL_fp = LL_fp(LL_Types_para_(:,2)==2,:);

LL_pf = LLPara(LL_Types_para(:,1)==2,:);
LL_Types_para_ = LL_Types_para(LL_Types_para(:,1)==2,:);
LL_pf = LL_pf(LL_Types_para_(:,2)==1,:);

LL_pp = LLPara(LL_Types_para(:,1)==2,:);
LL_Types_para_ = LL_Types_para(LL_Types_para(:,1)==2,:);
LL_pp = LL_pp(LL_Types_para_(:,2)==2,:);

L_ff = length(LL_ff);
L_fp = length(LL_fp);
L_pf =length(LL_pf);
L_pp =length(LL_pp);

L = length(LLPara);

n_f = sum(Types_para==1);
n_p = sum(Types_para==2);

C_ff = L_ff/n_f^2;
C_fp = L_fp/n_f/n_p;
C_pf = L_pf/n_f/n_p;
C_pp = L_pp/n_p^2;

%expected number of fp links:
EL_fp = C_ff*n_f*n_p;

%expected number of pp links
EL_pp = C_pf*n_p^2;

%percentage differences
fraction_fp = L_fp/EL_fp;

fraction_pp = L_pp/EL_pp;




%calculate difference in expected
%{
%Looking at vulnerability of niche space-ever been done?
Odegs = zeros(100*n_f,1);
ns = zeros(100*n_f,1);
n_f =50;
tic
for ii = 1:1000
[Res,Cons,n,c,r] = NicheModel_nk(n_f, C_ff); 

odeg = zeros(n_f,1);
if mod(ii,100)==0
    fprintf('%u\n',ii)
end

for jj = 1:n_f
    odeg(jj) = sum(Res==jj);
end

Odegs(((ii-1)*n_f+1):(ii*n_f)) = odeg;
ns(((ii-1)*n_f+1):(ii*n_f)) = n;

end
toc
vec = [ns,Odegs];
vec = sortrows(vec);

%

parts = 100;
avgs = zeros(parts,1);
dp = 1/parts;
tic
for ii = 1:parts
    gthanO = Odegs(ns>=(ii-1)/parts);
    ngthan = ns(ns>=(ii-1)/parts);
    glthanO = gthanO(ngthan<=(ii/parts));
    avgs(ii) = mean(glthanO);
end
toc
plot(0:dp:.99,avgs)
mean(avgs)
%}

%{
n[Res,Cons,n,c,r,Index_par] = NicheModel_nk...
    (n_f, C_ff, n_p, C_pf); 

nicheweb = sparse(Res,Cons,1);
basalsp=find(sum(nicheweb,2)==0);
topsp=find(sum(nicheweb)==0);

n_basal = n(basalsp);
n_top = n(topsp);
n_par = n(Index_par);
n([basalsp;topsp']) = [];


plot(n,zeros(size(n)),'ob',n_basal,zeros(size(n_basal)),'xg',...
    n_top,zeros(size(n_top)),'.r',n_par,zeros(size(n_par)),'ks')

L/(n_f+n_p)^2
C_web = sum(sum(nicheweb))/(n_f+n_p)^2

NL_ff=sum(sum(nicheweb(1:n_f,1:n_f)))
L_ff
NL_fp=sum(sum(nicheweb(1:n_f,(n_f+1):end)))
L_fp
NL_pf=sum(sum(nicheweb((n_f+1):end,1:n_f)))
L_pf
NL_pp=sum(sum(nicheweb((n_f+1):end,(n_f+1):end)))
L_pp
%}
