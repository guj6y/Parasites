function [local_properties,global_properties] = calculateProperties(res,con)

N = max(max(res,con));
A = sparse(res,con,1.0,N,N);
basalsp = full(sum(A)==0);
nonbasalsp = ~basalsp;
nicheweb = full(A);
nichewithoutcan = full(A);
nichewithoutcan(eye(N)>0) = 0;

%Betweenness centrality
options.unweighted = true;
bc = betweenness_centrality(A,options);

%Eigenvector centrality
[eigen_vec,eigen_val] = eigs(A,1,'LR');

%Clustering Coefficient(s)
cc = clustering_coefficients(A);

clust_WS = mean(cc);
clust_MN = trace(A^3)/(sum(sum(A^2))-trace(A^2));

%In degree (generalit)
gen = full(sum(A));
%And moments:
mom_in_1 = mean(gen);
mom_in_2 = moment(gen,2);
mom_in_3 = moment(gen,3);
mom_in_4 = moment(gen,4);
mom_in_5 = moment(gen,5);

%out degree (vulnerability)
vul = full(sum(A,2));
%and moments:
mom_out_1 = mean(vul);
mom_out_2 = moment(vul,2);
mom_out_3 = moment(vul,3);
mom_out_4 = moment(vul,4);
mom_out_5 = moment(vul,5);

%{
%Mean Chain Length (Min?)
%Loops (max? min?)
Loop=[]; %list of species involved in loops
loop_here = [];
    tmp=nichewithoutcan;
    m=zeros(length(nonbasalsp),length(nonbasalsp));
    chLen=0;
    chNum=0;
    nichewithoutloop=nichewithoutcan;
    for i=1:length(nonbasalsp) %path of length i
        loop_here=[];
        for j=1:length(nonbasalsp) %only the non basal species can eat other species
            if tmp(nonbasalsp(j),nonbasalsp(j))>0 %there's a loop of length i involving j
                Loop=[Loop nonbasalsp(j)];
                loop_here=[loop_here nonbasalsp(j)];
            end
        end
        nichewithoutloop(loop_here,loop_here)=0;
        tmp=tmp*nichewithoutcan;
    end
    tmp=nichewithoutloop;
    for i=1:length(nonbasalsp) %path of length i
        for j=1:length(nonbasalsp) %only the non basal species can eat other species
            s=sum(tmp(nonbasalsp(j),basalsp)); %number of chains of length i involving j
            m(j,i)=m(j,i)+s;
            chLen=chLen+s*i;
            chNum=chNum+s;
        end
        tmp=tmp*nichewithoutloop;
    end

    Loop=length(unique(Loop))/N;
    chLen=chLen/chNum;
    chNum=log(chNum);
%-------------------------------------


% ChSD --> Standard Deviation of ChLen
    chSD=0;
    for k=1:length(nonbasalsp) %for each non basal species
        p=0; %chain length
        num=0; %number of chains
        for l=1:length(nonbasalsp)
            p=p+l*m(k,l);
            num=num+m(k,l);
        end
        if num==0
            chSD=chSD+chLen^2;
        else
            chSD=chSD+(p/num-chLen)^2;
        end
    end
    chSD=sqrt(chSD/(N-1));
%}
%Mean vulnerability of prey
%Mean importance on prey

    mean_vul_prey = zeros(N,1);
    mean_imp_prey = zeros(N,1);
    
for ii = 1:N
    
    prey_ii = nicheweb(:,ii)>0;
    vul_prey_ii = sum(nicheweb(prey_ii,:),2);
    mean_vul_prey(ii) = mean(vul_prey_ii);
    mean_imp_prey(ii) = mean(1./vul_prey_ii);
    
end

%Mean generality of predators
%Mean importance to predators
mean_gen_pred = zeros(N,1);
mean_imp_pred = zeros(N,1);
for ii = 1:N

    pred_ii = nicheweb(ii,:)>0;
    gen_pred_ii = sum(nicheweb(:,pred_ii));
    mean_gen_pred(ii) = mean(gen_pred_ii);
    mean_imp_pred(ii) = mean(1./gen_pred_ii);
end
        
%shortest path to basal
[sp_mx, pred_mx] = floyd_warshall_all_sp(A);

min_sp_tobasal = zeros(N,1);
sp_tobasal_mx = sp_mx(basalsp,nonbasalsp);
sp_tobasal_mx(isinf(sp_tobasal_mx)) = NaN;
web_D = max(max(sp_tobasal_mx));
min_sp_tobasal(nonbasalsp) = min(sp_tobasal_mx,[],'omitnan');

%number of basal it is connected to
num_basal_con = zeros(N,1);
sp_tobasal_mx(isnan(sp_tobasal_mx)) = 0;

num_basal_con(nonbasalsp) = sum(sp_tobasal_mx>0);

local_properties = [cc,abs(eigen_vec),bc,gen',vul,mean_vul_prey,mean_imp_prey...
    mean_gen_pred,mean_imp_pred,min_sp_tobasal,num_basal_con];

global_properties = [clust_MN, clust_WS, abs(eigen_val), web_D, ...
    mom_in_1,mom_in_2,mom_in_3,mom_in_4,mom_in_5,...
    mom_out_1,mom_out_2,mom_out_3,mom_out_4,mom_out_5];
end