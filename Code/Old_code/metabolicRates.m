function [ATN_param] = metabolicRates(res,cons,N)
%ATNParameterGenerator gives an array with all parameters for ATN.
%   res = resource list
%   cons = consuemr list
%   n = niche values for all species
%   c = diet centers for all species
%   r = diet ranges for all species
%   N = number of species

%Allometric Trophic Network Model (ATN) is a bio-energetic model that
%allometrically scales predatory interactions and metabolic rates by a
%species' body mass (and ratios with its prey).  metabolic rates scale
%according to the type of species (fish vs. invertebrates, for example;
%these classifications are determined in random networks by a species'
%trophic level.  Trophic level is calculated here using  a "short-weighted
%measure" of trophic Level.(williams and martinez 2004a - "limits to
%trophic levels...")

%%%using all parameters from AAAI paper for now.-NK

%ATN has two distinct DiffEq for basal vs. non-basal species
   %(autotrophs vs. heterotrophs)
   %for basal species:
   %dB_i/dt = (Lotka-Volterra com'petition among consumers) - predation
   %parameters:
    %r_i: intrinsic growth rate ~N(0.6,0.2)
    %K: carrying capacity = 540
   %for consumers:
   %dB_i/dt = -metabolic loss + consumption - predation
   %parameters:
    %x_i: i's metabolic rate = f(M_i)
    %y_ij = 10
    %
   %Consumers can also be exploited by humans (future work).

%Species List:
S = 1:N;

%Need to first calculate the "short-weighted Trohpic level of each species"

%TL is also calculated using Network3d when reported as a network property.
%This is because I can't streamline the entire process on this machine; Not
%enough time to code all network properties in MatLab before comps.-NK

swTL = zeros(N,1);

nicheweb = sparse(res,cons,1,N,N);

%first find basal species (no prey); They are defined to have a TL of 1.
basal = sum(nicheweb)==0;  %logical indices of basal species
swTL(basal) = 1;

paTL_mx = nicheweb*(diag(1./sum(nicheweb)));

paTL = (eye(N)-paTL_mx')\ones(N,1);

clear paTL_mx

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

%At this point it is unclear whether we have a lognormal drawn for each
%species or just one for each web(&each type of species).  As of now I draw
%a new lognormal RV for each species.

mu_fish = 100;
sd_fish = 50;

mu_inv = 10;
sd_inv = 100;

muz_fish = log(mu_fish^2/sqrt(sd_fish^2+mu_fish^2));
sdz_fish = sqrt(log(sd_fish^2/mu_fish^2+1));

muz_inv = log(mu_inv^2/sqrt(sd_inv^2+mu_inv^2));
sdz_inv = sqrt(log(sd_inv^2/mu_inv^2+1));

Z_fish = lognrnd(muz_fish,sdz_fish,N,1);
Z_inv = lognrnd(muz_inv,sdz_inv,N,1);

M = (Z_fish.*fish).^(swTL-1) + (Z_inv.*inverts).^(swTL-1);
M(basal) = 1;

x = zeros(N,1);

x(inverts) = .314*M(inverts).^(-.24);
x(fish) = .88*M(fish).^(-.24);


ATN_param = [basal x];

end

