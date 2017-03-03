function [global_properties] = calculateGlobalProperties(res,con)

N = max(max(res,con));
A = sparse(res,con,1.0,N,N);
basalsp = full(sum(A)==0);
nonbasalsp = ~basalsp;



%Eigenvector centrality
[~,eigen_val] = eigs(A,1,'LR');

%Clustering Coefficient(s)
cc = clustering_coefficients(A);
clust_WS = mean(cc);
clust_MN = trace((A-eye(N))^3)/(sum(sum(A^2))-trace(A^2));

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
    
        
%shortest path to basal
[sp_mx, ~] = floyd_warshall_all_sp(A);

sp_tobasal_mx = sp_mx(basalsp,nonbasalsp);
sp_tobasal_mx(isinf(sp_tobasal_mx)) = NaN;
web_D = max(max(sp_tobasal_mx));

%connected Components
[Q,R] = graphconncomp(A);
inLoop = zeros(N,1);


for ii = 1:Q
    if sum(R==ii) >1
        inLoop(R==ii) = 1;
    end
end
spLoops = sp_mx + diag(diag(A)) + diag(inLoop);
fsp = sum(sum(spLoops >0))/N^2;
fLoop = sum(inLoop)/N;

global_properties = [clust_MN,...       1
                     clust_WS,...       2
                     abs(eigen_val),... 3
                     web_D, ...         4
                     mom_in_1,...       5
                     mom_in_2,...       6
                     mom_in_3,...       7
                     mom_in_4,...       8
                     mom_in_5,...       9
                     mom_out_1,...      10
                     mom_out_2,...      11
                     mom_out_3,...      12
                     mom_out_4,...      13
                     mom_out_5,...      14
                     fsp,...            15
                     fLoop];...         16
end