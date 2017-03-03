function [local_properties] = calculateLocalProperties(res,con,varargin)

if nargin >2
    N = varargin{1};
else
    N = max(max([res,con]));
end


A = full(sparse(res,con,1.0,N,N));
basalsp = full(sum(A)==0);
nonbasalsp = ~basalsp;
nicheweb = full(A);
S = 1:N;
L = sum(sum(nicheweb));


%shortest path to basal
[sp_mx] = johnson_all_sp(sparse(A));

min_sp_tobasal = zeros(N,1);
sp_tobasal_mx = sp_mx(basalsp,nonbasalsp);
sp_tobasal_mx(isinf(sp_tobasal_mx)) = NaN;
min_sp_tobasal(nonbasalsp) = min(sp_tobasal_mx,[],'omitnan');

%number of basal it is connected to
num_basal_con = zeros(N,1);
sp_tobasal_mx(isnan(sp_tobasal_mx)) = 0;

num_basal_con(nonbasalsp) = sum(sp_tobasal_mx>0);
%Clustering Coefficient(s)
cc = clustering_coefficients(sparse(A));

%In degree (generalit)
gen = full(sum(A))/(L/N);

%out degree (vulnerability)
vul = full(sum(A,2))/(L/N);

%short-weighted trophic level
nonCannA = A.*(1-eye(N));

%Formula for the prey-averaged trophic level.
paTL_mx = sparse(nonCannA)*(diag(1./sum(sparse(nonCannA))));

paTL = (speye(N)-paTL_mx')\ones(N,1);

%Short trophic level is shortest path to a basal; since I've already
%calculated shortest paths this is easy >;* <<<3 xoxo -nk
sTL = (min(sp_mx(basalsp,:))+1)';

swTL = (paTL+sTL)/2;

%Mean vulnerability of prey
%Mean topological fraction of preys' consumers
mean_vul_prey = zeros(N,1);
mean_imp_prey = zeros(N,1);
max_vul_prey = zeros(N,1);
min_vul_prey = zeros(N,1);

median_patl_prey_diff = zeros(N,1);
median_swtl_prey_diff = zeros(N,1);
for ii = 1:N
    
    prey_ii = nicheweb(:,ii)>0;    
    vul_prey_ii = sum(nicheweb(prey_ii,:),2);
    mean_vul_prey(ii) = mean(vul_prey_ii);
    patl_prey_ii = paTL(prey_ii);
    swtl_prey_ii = swTL(prey_ii);
    
    if sum(prey_ii) == 0
        median_patl_prey_diff(ii) = nan;
        median_swtl_prey_diff(ii) = nan;
    else
        median_patl_prey_diff(ii) = paTL(ii) - median(patl_prey_ii);
        median_swtl_prey_diff(ii) = swTL(ii) - median(swtl_prey_ii);
    end
    
    if numel(vul_prey_ii) >0
    max_vul_prey(ii) = max(vul_prey_ii);
    min_vul_prey(ii) = min(vul_prey_ii);
    else
        max_vul_prey(ii) = nan;
        min_vul_prey(ii) = nan;
    end
    
    mean_imp_prey(ii) = mean(1./vul_prey_ii);
    
end

%Mean generality of predators
%Mean topological fraction of predators' diet
mean_gen_pred = zeros(N,1);
mean_imp_pred = zeros(N,1);
max_gen_pred = zeros(N,1);
min_gen_pred = zeros(N,1);

median_patl_pred_diff = zeros(N,1);
median_swtl_pred_diff = zeros(N,1);

for ii = 1:N
    
    pred_ii = nicheweb(ii,:)>0;
    gen_pred_ii = sum(nicheweb(:,pred_ii));
    mean_gen_pred(ii) = mean(gen_pred_ii);
    patl_pred_ii = paTL(pred_ii);
    swtl_pred_ii = swTL(pred_ii);
    
    if sum(pred_ii) == 0
        median_patl_pred_diff(ii) = nan;
        median_swtl_pred_diff(ii) = nan;
    else
        median_patl_pred_diff(ii) = paTL(ii) - median(patl_pred_ii);
        median_swtl_pred_diff(ii) = swTL(ii) - median(swtl_pred_ii);
    end
    
    if numel(gen_pred_ii)>0
        max_gen_pred(ii) = max(gen_pred_ii);
        min_gen_pred(ii) = min(gen_pred_ii);
    else
        max_gen_pred(ii) = nan;
        min_gen_pred(ii) = nan;
    end
    mean_imp_pred(ii) = mean(1./gen_pred_ii);
    
end




%Number of species in loops is the number of species in strongly connected
%components larger than one species.
[n_comp,comps] = graphconncomp(sparse(A),'Directed',1,'Weak',0);
num_comp = zeros(1,n_comp);
sp_in_loop = [];

for ii = 1:n_comp
    num_comp(ii) = sum(comps==ii);
    if num_comp(ii) > 1
        sp_in_loop = [sp_in_loop; S(comps==ii)'];
    end
end
inLoop = zeros(N,1);
inLoop(sp_in_loop) = 1;

top_sp = sum(A,2)==0;
%basal?
mean_vul_prey(basalsp) = nan;
max_vul_prey(basalsp) = nan;
min_vul_prey(basalsp)= nan;

mean_imp_prey(basalsp) = nan;

mean_gen_pred(top_sp) = nan;
mean_imp_pred(top_sp) = nan;
max_gen_pred(top_sp) = nan;
min_gen_pred(top_sp) = nan;


                    %1  %2    %3  %4              %5
local_properties = [cc,gen',vul,mean_vul_prey,mean_imp_prey...
    ...%6             %7            %8             %9        %10   %11
    ,mean_gen_pred,mean_imp_pred,min_sp_tobasal,num_basal_con,swTL,inLoop...
    ...12               13          14              15
    ,max_gen_pred, min_gen_pred, max_vul_prey, min_vul_prey...
    ...   %16                 %17
    ,median_patl_prey_diff,median_patl_pred_diff,...
    ...  %18                      %19
    median_swtl_prey_diff,median_swtl_pred_diff];

end