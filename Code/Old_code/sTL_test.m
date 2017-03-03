N=20;
M=full(M);


basal = sum(M)==0;

sTL=zeros(20,1);

sTL(basal) = 1;
k = 1;
while sum(sTL==0)
    Mk = M^k;
    k=k+1;
    contobasal = sum(Mk(basal,:))~=0;
    sTL(contobasal&(sTL==0)') = k;
end