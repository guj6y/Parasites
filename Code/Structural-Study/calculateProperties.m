function [props] = calculateProperties(res,con)

%This functions calculates a lot of standard properties for food webs.
%Also, a few non-standard properties.  Mostly, the non-standard properties
%are higher moments of degree distributions.
%
%
% Depends: MatLabBGL. (boost graph library for matlab)
%     I use "graphconncomp" and "floyd_warshall_all_sp" to calculate
%     connected components and shortest paths in the food web.
% They aren't too bad to code up a reasonably efficient version.  If you
% don't/can't have BGL.
%

N = max(max(res,con));
S = 1:N;
A = sparse(res,con,1.0,N,N);

%Maximum Similarity: need to compare all species to each other -gross! >.<

%Thought - this is obviously a way to measure competition.  Probably a
%reference about that?  Are species with high similiarity (or just
%prey/pred similarity) more or less likely to go extinct?

%NK's attempt to improve max sim with ridiculous 3d arrays:
co_prey_mxs = zeros(N,N,N);
co_pred_mxs = zeros(N,N,N);
any_prey_mxs = zeros(N,N,N);
any_pred_mxs = zeros(N,N,N);

%might be possible to reduce complexity with clever indexing.  Still, less
%than a second for up to about 250x250 matrices.  This might also be slower
%depending on the machine it is on - less RAM might mean this doesn't even
%work.  Probalby fine on any modern laptop.
%Method of Coralie & Perrine (significantly slower when I tested it):
%{
S_ij = zeros(N,N);
for ii = 1:(N-1)
    for jj = (ii+1):N
        %number Prey in common:
        co_prey = sum((A(:,ii)>0) & (A(:,jj))>0);
        %Predators in common:
        co_pred = sum((A(ii,:)>0) & (A(jj,:)>0));
        %Total links for both species:
        tot_prey = sum((A(:,ii)>0) | (A(:,jj))>0) +...
            sum((A(ii,:)>0) | (A(jj,:)>0));
        %Saving the similarity.
        S_ij(ii,jj) = (co_prey + co_pred)/tot_prey;
        S_ij(jj,ii) = (co_prey + co_pred)/tot_prey;
    end
end
%}


for ii = 1:N
    
    %co_prey_mxs(ii,jj,kk) = 1 means species kk (3rd dim) and species jj
    %(2nd dim) both eat species ii (1st dim).
    co_prey_mxs(:,:,ii) = ((A(:,ii)*ones(1,N))>0)&(A>0);
    
    %co_pred_mxs(ii,jj,kk) = 1 means species kk (3rd dim) and species ii
    %(1st dim) both are eaten by species jj (2nd dim).
    co_pred_mxs(:,:,ii) = ((A(ii,:)'*ones(1,N))>0)'&(A>0);
    
    
    %any_prey_mxs(ii,jj,kk) = 1 means species kk (3rd dim) or species jj
    %(2nd dim) bot eat species ii (1st dim).
    any_prey_mxs(:,:,ii) = ((A(:,ii)*ones(1,N)>0)|(A>0));
    
    
    %any_pred_mxs(ii,jj,kk) = 1 means species kk (3rd dim) or species ii
    %(1st dim) is eaten by species jj (2nd dim).
    any_pred_mxs(:,:,ii) = ((A(ii,:)'*ones(1,N)>0)'|(A>0));
    
    
    %Tried reducing the complexity by filling in at 'right angles' ( fill
    %the 3rd dimension, then a strip in one of the others at each time
    %step).  Ended up increasing computation time ever so slightly.  oh
    %well.  ¯\_(?)_/¯
    %{
    %co_prey_mxs(ii,jj,kk) = 1 means species kk (3rd dim) and species jj
    %(2nd dim) both eat species ii (1st dim).
    co_prey_mxs(:,ii:N,ii) = ((A(:,ii)*ones(1,(N-ii+1)))>0)&(A(:,ii:N)>0);
    co_prey_mxs(:,ii,(ii+1):N) = co_prey_mxs(:,(ii+1):N,ii);
    %co_pred_mxs(ii,jj,kk) = 1 means species kk (3rd dim) and species ii
    %(1st dim) both are eaten by species jj (2nd dim).
    co_pred_mxs(ii:N,:,ii) = ((A(ii,:)'*ones(1,(N-ii+1)))>0)'&(A(ii:N,:)>0);
    co_pred_mxs(ii,:,(ii+1):N) = co_pred_mxs((ii+1):N,:,ii)';
    
    %any_prey_mxs(ii,jj,kk) = 1 means species kk (3rd dim) or species jj
    %(2nd dim) bot eat species ii (1st dim).
    any_prey_mxs(:,ii:N,ii) = ((A(:,ii)*ones(1,(N-ii+1)))>0)|(A(:,ii:N)>0);
    any_prey_mxs(:,ii,(ii+1):N) = any_prey_mxs(:,(ii+1):N,ii);
    
    %any_pred_mxs(ii,jj,kk) = 1 means species kk (3rd dim) or species ii
    %(1st dim) is eaten by species jj (2nd dim).
    any_pred_mxs(ii:N,:,ii) = ((A(ii,:)'*ones(1,(N-ii+1)))>0)'|(A(ii:N,:)>0);
    any_pred_mxs(ii,:,(ii+1):N) = any_pred_mxs((ii+1):N,:,ii)';
    %}
end
%Had to fiddle a *bit* to get the dimensions right.. but this agrees with
%our french colleagues' method.
S_ij = (permute(sum(co_prey_mxs),[2,3,1]) +...
    permute(sum(co_pred_mxs,2),[1,3,2]))./(...
    (permute(sum(any_prey_mxs),[2,3,1])) ...
    + permute(sum(any_pred_mxs,2),[1,3,2]))...
    -eye(N);


mean_max_sim = mean(max(S_ij));


%Eigenvector centrality
[~,svd1] = svds(A,1);

%Clustering Coefficient
clust_MN = trace(A^3)/(sum(sum(A^2))-trace(A^2));

%In degree (generality)
gen = full(sum(A));
%And moments:
gen_1 = mean(gen);
gen_2 = moment(gen,2);
gen_3 = moment(gen,3);
gen_4 = moment(gen,4);
gen_5 = moment(gen,5);

%out degree (vulnerability)
vul = full(sum(A,2));
%and moments:
vul_1 = mean(vul);
vul_2 = moment(vul,2);
vul_3 = moment(vul,3);
vul_4 = moment(vul,4);
vul_5 = moment(vul,5);

%Degree (links)
link = gen+vul';
%and moments:
link_1 = mean(link);
link_2 = moment(link,2);
link_3 = moment(link,3);
link_4 = moment(link,4);
link_5 = moment(link,5);

%Number of species in loops is the number of species in strongly connected
%components larger than one species.

[n_comp,comps] = graphconncomp(A,'Directed',1,'Weak',0);
num_comp = zeros(1,n_comp);
sp_in_loop = [];

for ii = 1:n_comp
    num_comp(ii) = sum(comps==ii);
    if num_comp(ii) > 1
        sp_in_loop = [sp_in_loop; S(comps==ii)'];
    end
end
num_sccs = sum(num_comp>1);
f_in_sccs = length(sp_in_loop)/N;

%Fraction basal species
basal =  sum(A) == 0;
f_basal = sum(basal)/N;

%Fraction herbivores
bin_herb = sum(A(basal,:)) == sum(A);
f_herb = sum(bin_herb)/N;

%fraction top species
top = sum(A,2) == 0;
f_top = sum(top)/N;

%Fraction intermediate species
f_int = 1-f_top-f_basal;

%Fraction cannibals
cann = res(res==con);
f_cann = length(cann)/N;

%Mean shortest path between two species.
[sp_mx, ~] = floyd_warshall_all_sp(A);

mean_sp = mean(sp_mx(sp_mx<Inf));

%Fraction of paths between two species.
f_sp = sum(sum(sp_mx < Inf))/(N*(N-1));

%Diameter of the web
diam = max(sp_mx(sp_mx<Inf));

%calculate the Trophic Levels

%Formula for the prey-averaged trophic level.
paTL_mx = A*(diag(1./sum(A)));

paTL = (speye(N)-paTL_mx')\ones(N,1);

%Short trophic level is shortest path to a basal; since I've already
%calculated shortest paths this is easy >;* <<<3 xoxo -nk
sTL = (min(sp_mx(basal,:))+1)';

swTL = (paTL+sTL)/2;

%Some metrics of swTL distribution.
mean_swTL = mean(swTL);
std_swTL = std(swTL);
skew_swTL = skewness(swTL);
max_swTL = max(swTL);

%Omnivory isn't this easy?
omn = (mod(swTL,1)~=0)&(gen>=2)';
f_omn = sum(omn)/N;



props = full([svd1,clust_MN,...
    gen_1,gen_2,gen_3,gen_4,gen_5,...
    vul_1,vul_2,vul_3,vul_4,vul_5,...
    link_1,link_2,link_3,link_4,link_5,...
    num_sccs,f_in_sccs,...
    f_basal,f_herb,f_omn,f_top,f_int,f_cann,...
    mean_sp,f_sp,diam,...
    mean_swTL,std_swTL,skew_swTL,max_swTL...
    mean_max_sim,...
    ]);

end