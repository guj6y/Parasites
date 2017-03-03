function [SOL,r,M,x,swTL,spTypes,y0] = ATN_model_nk(res,cons)
%ATN_model runs the allometric trophic network model given food web.  Takes
%a vector of resources and vector of consumers (for each link) as input.
%
% 
%
% [SOL,r,M,x,swTL,spTypes] = ATN_model_nk(res,cons) gives as output:
% SOL:  solution structure from ode45
% r:    growth rates of basal species
% M:    body sizes (ratio to basal species)
% x:    metabolix rates
% swTL: short weighted trophic level; ref = williams and martinez 2004 - "limits to
% trophic levels..."
% spTypes:  Cell array with the type of each species (basal,invertebrate,
% or fish)
% y0: the vector of initial biomass densities
%
% res,cons are resources and consumers, respectively.  This is a
% linklist,i.e. it represents feeding links by res(i) -> cons(i).
%
% ATN requires a few parameters.  These are currently hard coded from Martinez et al
% 2012 in AAAI.
% 

%Allometric Trophic Network Model (ATN) is a bio-energetic model that
%allometrically scales predatory interactions and metabolic rates by a
%species' body mass (and ratios with its prey).  metabolic rates scale
%according to the type of species (fish vs. invertebrates, for example;
%these classifications are determined in random networks by a species'
%trophic level.  Trophic level is calculated here using  a "short-weighted
%measure" of trophic Level.(williams and martinez 2004 - "limits to
%trophic levels...")

% s = RandStream('mt19937ar','Seed',3);
% RandStream.setGlobalStream(s);

%Species List:
N = max(max([res cons]));
S = 1:N;

%TL is also calculated using Network3d when reported as a network property.
%This is because I can't streamline the entire process on this machine; Not
%enough time to code all network properties in MatLab before comps.-NK

swTL = zeros(N,1);

nicheweb = sparse(res,cons,1,N,N);

%first find basal species (no prey); They are defined to have a TL of 1.
basal = sum(nicheweb)==0;  %logical indices of basal species
swTL(basal) = 1;

%Short weighted TL requires calculation of prey averaged TL:
paTL_mx = nicheweb*(diag(1./sum(nicheweb)));

paTL = (eye(N)-paTL_mx')\ones(N,1);

clear paTL_mx

%Short weighted TL also requires calculation of shortest path to a basal
%Species.  This is done using powers of the adjacency matrix (probably a
%better way exists)
sTL_mx=full(nicheweb);

sTL=zeros(N,1);

sTL(basal) = 1;
k = 1;
while sum(sTL==0)
    sTL_mxk = sTL_mx^k;
    k=k+1;
    contobasal = sum(sTL_mxk(basal,:))~=0;
    sTL(contobasal&(sTL==0)') = k;
end

clear sTL_mx

%average the two
swTL = (paTL+sTL)/2;

%Body masses:
M = zeros(N,1);

M(basal) = 1;

%Here we need to decide what type each species is.

spTypes = cell(N,1);

%fraction species with TL>3 that are fish:
fish_frac=0.6;
fish_rand = rand(N,1);
fish = (swTL>=3)&(fish_rand<fish_frac);
inverts = (swTL>1)&(swTL<3)+((swTL>=3)&(fish_rand>=fish_frac));

[spTypes{swTL==1}] = deal('Basal');
[spTypes{inverts}]=deal('Invertebrate');
[spTypes{fish}] = deal('Fish');

%%%Draw a lognormal RV for each species.

mu_fish = 100;
sd_fish = 50;

mu_inv = 10;
sd_inv = 100;
%lognrnd takes the mean, stdev of the corresponding normal as input
muz_fish = log(mu_fish^2/sqrt(sd_fish^2+mu_fish^2));
sdz_fish = sqrt(log(sd_fish^2/mu_fish^2+1));

muz_inv = log(mu_inv^2/sqrt(sd_inv^2+mu_inv^2));
sdz_inv = sqrt(log(sd_inv^2/mu_inv^2+1));

Z_fish = lognrnd(muz_fish,sdz_fish,1);
Z_inv = lognrnd(muz_inv,sdz_inv,1);


%Calculate body sizes
M = (Z_fish.*fish).^(swTL-1) + (Z_inv.*inverts).^(swTL-1);
M(basal) = 1;

%Calculate metabolic rates
x = zeros(N,1);

x(inverts) = .314*M(inverts).^(-.24);
x(fish) = .88*M(fish).^(-.24);
diet_sizes = sum(nicheweb);


%calculating preference of j towards i (wij):
%row vector of 1/(i's diet size) (:= 0 if basal).
wi = 1./diet_sizes;
wi(isinf(wi)) = 0;

%%%CAUTION THIS IS THE TRANSPOSE OF HOW IT IS DONE IN OTHER CODES%%%%%%%%%

wij = nicheweb*diag(wi);

%A few more parameters:

%hill exponenet determines the shape of the functional response
h = 1.2;

%ones vector and matrix
one = ones(N,1);
onem = ones(N,N);

%ecosystem carrying capacity
K = 540;

%half-saturation density
B0 = 80;

%pre-calculated B0 to hill exponent.
B0h = B0^h;

%self interference term.  =0or1
c=0;

%growth rates of basal species
r = randn(N,1)*.2+.6;
r(~basal) = 0;

%Maximum consumtion rate of j eating i.
yij = 10*nicheweb;

%maximum assimilation of j eating i;
eij = nicheweb;
%for herbivory
eij(basal,:) = eij(basal,:)/.66;
%for carnivory
eij(~basal,:) = eij(~basal,:)/.85;
%reciprocal because it only occurs as reciprocal in equations.
    
    function [dB] = ATN_rhs(t,B)
        
        B(B<1e-6) = 0;
            
        F = diag(B.^h)*wij*diag(1./(one*B0h + c*B*B0h + wij'*B.^h));
        
        dB = (r.*(one-basal*B/K) - x +(yij.*F)'*one.*x).*B-(F.*yij.*eij)*(x.*B);
        
        
        if~isreal(dB)
            fprintf('non-real derivative\n')
        end
        
    end

% s = RandStream('mt19937ar','Seed','shuffle');
% RandStream.setGlobalStream(s);

y0 = 5+rand(N,1)*45;
tic
SOL = ode45(@ATN_rhs,[0 2000],y0);
toc
figure
hold on
plot(SOL.x,SOL.y(basal,:),'g')
plot(SOL.x,SOL.y(fish,:),'b')
plot(SOL.x,SOL.y(inverts,:),'r')
sum(SOL.y(:,end)>1e-6)

end


