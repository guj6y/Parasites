close all
M = 1000;
C = .12;
Cs = zeros(1,M);
S = 100;
for ii = 1:M
    ii
[Res,~,~,~,~] = NicheModel_nk_(S,C);
Cs(ii) = numel(Res)/(S^2);

end

[N,X] = hist(Cs);
w = X(2)-X(1);
b = (1-2*C)/(2*C); 
hold on
variance = C-(1/(1+b))*(2/(2+b))/3
[H,P] = kstest((Cs-C)/sqrt(variance/M))
plot(-4:.1:4,pdf('normal',-4:.1:4,0,1))
bar((X-C)/sqrt(variance/M),N/M)
