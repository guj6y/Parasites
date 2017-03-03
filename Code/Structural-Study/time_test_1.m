t = zeros(100,1);
for ii = 1:1000
    ii
    tic
    NicheModel_nk(100,.12,15);
    t(ii) = toc;
    
    
end

mean(t)