N= 100;
t = zeros(1,N);
num_bad_ONM=0;
for ii = 1:N
    tic
[res0, cons0, n0, c0, r0,] = NicheModel_nk(Nf,C);
        [res2, cons2, n2, c2, r2] = ...
            AddParasites1(n0,c0,r0,Np,[C,Cfp],1);
    t(ii) = toc;
end

mean(t)
num_bad_ONM