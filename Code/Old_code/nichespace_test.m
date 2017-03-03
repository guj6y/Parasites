%
close all
clear all



Trials = 100;
S = 100;
C = 0.05:.01:.15;
basalsc = cell(Trials,1);
topsc = cell(Trials,1);
inters = cell(Trials,1);
swTLs = zeros(S,Trials);

basals = zeros(S,Trials);
herbs = zeros(S,Trials);
omns = zeros(S,Trials);
carns = zeros(S,Trials);
tops = zeros(S,Trials);
canns = zeros(S,Trials);
gens = zeros(S,Trials);
vuls = zeros(S,Trials);
sTLs = zeros(S,Trials);

ns = zeros(S,Trials);
t = zeros(Trials,1);

ress = cell(Trials,1);
conss = cell(Trials,1);
hold on

maxswTL = zeros(Trials,1);
avg_max_swTLs = zeros(length(C),1);
for jj = 1:length(C)
for ii = 1:Trials
    
    tic
    [res,cons,n,c,r] = NicheModel_nk(S, C(jj));
    ress{ii} = res;
    conss{ii} = cons;
    mx = sparse(res,cons,1,S,S);
%     [mx,n,c,r]= NicheModel(S,C);
%     mx = mx';
    t(ii) = toc;
    
    
    basal = find(sum(mx)==0);
    top = find(sum(mx,2)==0)';
    slist = 1:S;
    slist([basal,top]) = [];
    inter = slist;
    
    %{
    basalsc{ii} = basal;
    topsc{ii} = top;
    inters{ii} = inter;
    %}
    
    ns(:,ii) = n;
    
    swTL = zeros(S,1);

    
%first find basal species (no prey); They are defined to have a TL of 1.
basal = sum(mx)==0;  %logical indices of basal species
swTL(basal) = 1;

paTL_mx = mx*(diag(1./sum(mx)));

paTL = (eye(S)-paTL_mx')\ones(S,1);
%paTL = 0;

sTL_mx=full(mx);

sTL=zeros(S,1);

sTL(basal) = 1;
k = 1;
while sum(sTL==0)
    sTL_mxk = sTL_mx^k;
    k=k+1;
    contobasal = sum(sTL_mxk(basal,:))~=0;
    sTL(contobasal&(sTL==0)') = k;
        
end
swTL = (paTL+sTL)/2;

swTLs(:,ii) = swTL;
maxswTL(ii) = max(swTL);
%{
carn = full((sum(mx(~basal,:))==sum(mx))~=basal);
herb = swTL==2;
omn = (~carn')&(~herb)&(~basal');
cann = diag(mx)~=0;
top = sum(mx,2)==0;

subplot(2,1,1)
hold on
plot(n(basal),swTL(basal),'g.')
plot(n(herb),swTL(herb),'b.')
plot(n(carn),swTL(carn),'r.')
plot(n(omn),swTL(omn),'m.')
plot(n(top),swTL(top),'k^')
plot(n(cann),swTL(cann),'bo')


subplot(2,1,2)
hold on
plot(n(basal),sTL(basal),'g.')
plot(n(herb),sTL(herb),'b.')
plot(n(carn),sTL(carn),'r.')
plot(n(omn),sTL(omn),'m.')
plot(n(top),sTL(top),'k^')
plot(n(cann),sTL(cann),'bo')


basals(:,ii) = basal';
herbs(:,ii) = herb;
omns(:,ii) = omn;
carns(:,ii) = carn;
tops(:,ii) = top;
canns(:,ii) = cann;
gens(:,ii) = sum(mx)';
vuls(:,ii) = sum(mx,2)';
sTLs(:,ii) = sTL;
%}
% subplot(2,2,3)
% plot(n,sum(mx),'k.')
% subplot(2,2,4)
% hist(sum(mx))
end
%legend('basasl','herb','carn','top','omn')
avg_max_swTLs(jj) = mean(maxswTL);
end
plot(C,avg_max_swTLs,'.')
%%

%
%
%meansb = zeros(Trials,1);
basal_all = [];

count_with_top = 0;
hold on
for ii = 1:Trials
    nii = ns(:,ii);
    %meansb(ii) = mean(nii(basalsc{ii}));
    %basal_all = [basal_all nii(basalsc{ii})'];
    %plot(ii,nii(basals{ii})')
    if ~isempty(topsc{ii})
        count_with_top=count_with_top+1;
    end
end
%title('scatterplot of (n_i,swTL_i)')
% figure
% hist(meansb)
% title('histogram of average niche value of basal species; agrees with CLT for Uniform')
% figure
% hist(basal_all)
% title('histogram of niche value of all basal species')
%}


%%
%{
figure
hist(ns(logical(basals)))
title('basal niche values')
figure
hist(ns(logical(tops)))
title('tops niche values')
figure
hist(ns(logical(carns)))
title('carns niche values')
figure
hist(ns(logical(herbs)))
title('herbs niche values')
figure
hist(ns(logical(omns)))
title('omns niche values')
figure
hist(ns(logical(canns)))
title('cannsv niche values')

figure
hist(ns((swTLs>=2)&(swTLs<3)))
title('2<=swTL<3')

figure
hist(ns((swTLs>=3)&(swTLs<4)))
title('3<=swTL<4')

figure
hist(ns((swTLs>=4)))
title('4<=swTL')

figure
hist(ns((sTLs>=4)))
title('4<=sTL')

figure
hist(ns((sTLs==3)))
title('3=sTL')

nsv = reshape(ns,Trials*S,1);
%
vulsv = reshape(vuls,Trials*S,1);
gensv = reshape(gens,Trials*S,1);

figure
plot(nsv,vulsv,'b.')
title('vulnerability of niche space')

figure
plot(nsv,gensv,'b.')
title('generality of niche space')

figure
hist(gensv)
title('Generality')

figure
hist(vulsv)
title('vulnerability')
%}



