t = zeros(100,1);

%
for ii = 1:100
    
    tic;
    [a b c d e] = NicheModel_nk(100,.12);
    AddParasites1(a,b,c,d,e,50,[.12 .15]);
    t(ii) = toc;
    ii

end
%}
%average time of successful webs (time < 10 s) was 1.0834.
%{
for ii = 1:100
    
    tic;
    [a b c d e] = NicheModelParasites([100 50], [.12 .15]);
    t(ii) = toc;
    ii

end
%}