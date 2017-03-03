N = 100;
rOutOuts = zeros(N,1);
rOutIns = zeros(N,1);
rInOuts = zeros(N,1);
rInIns = zeros(N,1);

scriptSs = zeros(N,1);
scriptHs = zeros(N,1);

means = zeros(N,1);
for ii = 1:N
    
    [res, con] = NicheModel_nk(147,0.1);
    
    %Assortativity. -- in degree = generaltiy; out degree = vulnerability
    A = sparse(res,con,1,max([res;con]),max([res;con]));
    gen = sum(A)';
    vul = sum(A,2);
    sourceInDegree = gen(res);
    sourceOutDegree = vul(res);
    targetInDegree = gen(con);
    targetOutDegree = vul(con);
    
    rOutOut = corr([sourceOutDegree,targetOutDegree]);
    rOutIn = corr([sourceOutDegree,targetInDegree]);
    rInOut = corr([sourceInDegree,targetOutDegree]);
    rInIn = corr([sourceInDegree,targetInDegree]);
    
    rOutOuts(ii) = rOutOut(2);
    rOutIns(ii) = rOutIn(2);
    rInOuts(ii) = rInOut(2);
    rInIns(ii) = rInIn(2);
    
    scriptS = (mean(gen.*vul)-mean(gen)*mean(vul))/(std(gen)*std(vul));
    scriptH = (std(gen)*std(vul))/mean(vul);
    scriptSs(ii) = scriptS;
    scriptHs(ii) = scriptH;
    
    means(ii) = mean(vul);
end

mean(rOutOuts)
mean(rInIns)
mean(rOutIns)
mean(rInOuts)
mean(scriptSs)
mean(scriptHs)
mean(means)
