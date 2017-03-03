nS = 8;
sFactor = 1.5;
t1 = zeros(nS,1);
t2 = zeros(nS,1);
t3 = zeros(nS,1);
S0 = 20;
Ss = round(20*sFactor.^(0:nS-1));
counter = 0;

sparseFactor = .1;
for S = Ss
counter = counter+1;
A = rand(S);

% sparsify to 10%:

A = A<sparseFactor;

wij = sparse(A);

B = rand(S,1);
h = .2;
h2 = 6/5;
B0h = .5;
one = zeros(S,1);

[res,con] = find(wij);

%
tic
for ii = 1:1000
     F1=(B.^h)./((one*B0h + A'*(B.^h)))'.*wij;
end

t3(counter) = toc;
%}

tic 
for ii =1:1000
    Bh = B.^(1+h);
    den =((one*B0h + wij'*(Bh)));
    F0 = Bh(res)./den(con);
    F3 = sparse(res,con,F0,S,S);

end
t1(counter) = toc;
%}
tic 
for ii =1:1000
    Bh = B.^(1+h);
    den =((one*B0h + wij'*(Bh)));
    
    F3 = sparse(res,con,Bh(res)./den(con),S,S);

end

t2(counter) = toc;

end
x = Ss.^2.*sparseFactor;
loglog(x,t1,x,t2)
loglog(Ss,t1,Ss,t2,Ss,t3)
%{
tic 
for ii =1 :100
    A.*wij;
    
end
toc

tic 
for ii =1 :100
    A.*A;
    
end
toc

tic 
for ii =1 :100
    wij.*wij;
    
end
toc
%}