function [web_mx,n,c,r,Cpf_new,Cpp_new] = makeRawWeb_para(Nf,Np,C4,p_range)

%Connectance for f->f:
Cff = C4(1);

%Connectacnes for p->f:
Cpf = C4(2);

%Connectance for f->p:
Cfp = C4(3);

%Connectances for p->p:
Cpp = C4(4);


N=Nf+Np;
np_min = p_range(1);
np_max = p_range(2);
Enp = (np_min+np_max)/2;

%This chunk of code generates our initial web.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Free-liver values
    %assign niche values from a uniform distribution
    nf = rand(Nf,1);  
    
    %designate range for each species

    %parameters for beta distribution (Free livers on free livers):
    alpha = 1;
    beta = (0.5-Cff)/(Cff);
    rf = betarnd(alpha,beta,Nf,1); 
    
    %vector of ranges:
    rf = rf.*nf;  
    % set center of range, uniformly distributed in [r_i/2,n_i]; 
    cf = min(1-rf./2,rand(Nf,1).*(nf-rf./2)+rf./2);

    %sort everything:
    [nf_new, Indx] = sort(nf);                                                
    %n_new: niche values in ascending order
    %Indx: indices of species in descending order 

    %the smallest r to highest index species 
    rf_new = rf(Indx);                                                       %NK: Maybe not how I'd do this
    cf_new = cf(Indx); 

    rf_new(1) = 0; %change the r of highest index species to 0
    %so we have a basal species in every web
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Parasite values
    %assign parasitic niche values from a uniform distribution
    np = rand(Np,1)*(np_max-np_min) + np_min;  
    
    
    %designate range for each parasite

    %parameters for beta distribution:
    alpha = 1;
    beta = (Enp-Cfp)/(Cfp); %%%%%double check
    rp = betarnd(alpha,beta,Np,1); 
    
    %vector of ranges:
    rp = rp.*np;  

    % set center of range, uniformly distributed in [n_p,1-r_p/2]; 
    cp = rand(Np,1).*(1-rp/2-np)+np;

    %sort everything:
    [np_new, Indx] = sort(np);                                                
    %np_new: niche values in ascending order
    %Indx: indices of species in descending order 

    rp_new = rp(Indx);                                                       %NK: Maybe not how I'd do this
    cp_new = cp(Indx); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Put parasites in with free-livers
    n = [nf_new;np_new];
    c = [cf_new;cp_new];
    r = [rf_new;rp_new];
    
    %Generate a niche web; at this point we should be aiming at (and pretty
    %much hitting) Cfp,Cff
    %
    preymins = c-r/2;
    preymaxs = c+r/2;
    
    preymins_mx = ones(N,1)*preymins';
    preymaxs_mx = ones(N,1)*preymaxs';

    %fills the empty matrix with niche ranges:
    n_mx = n*ones(1,N ); 
    
    web_mx=((n_mx>=preymins_mx)+(n_mx<=preymaxs_mx)==2*ones(N));
    
    %%%Looking at parasite->Free liver links
    pf_diets = sum(web_mx((Nf+1):(Nf+Np),1:Nf));
    
    %number of p->f links
    Lpf = sum(pf_diets);
    
    %number of free livers eating parasites
    Npf = sum(pf_diets>0);
    
    %Cpf connectance
    Cpf_web = Lpf/Nf/Np;
    
    %average width of free-liver overlap in parasitic range
    widthpf = Lpf/Npf/Np;
    
    widthpf_target = Cpf*Nf/Npf;
    
    %%%Looking at hyper-parasitism
    pp_diets = sum(web_mx((Nf+1):(Nf+Np),(Nf+1):(Nf+Np)));
    Lpp = sum(pp_diets);
    Npp = sum(pp_diets>0);
    widthpp = Lpp/Npp/Np;
    
    widthpp_target = Cpp*Np/Npp;
    Cpf_new = widthpf_target/widthpf*Cpf;
    Cpp_new = widthpp_target/widthpp*Cpp;
    
    %pf links
    preyminspf = preymins(1:Nf);
    preymaxspf = preymaxs(1:Nf);
    
    preyminspf((pf_diets>0)') = cf_new((pf_diets>0)')-rf_new((pf_diets>0)')/2*widthpf_target/widthpf;
    preymaxspf((pf_diets>0)') = cf_new((pf_diets>0)')+rf_new((pf_diets>0)')/2*widthpf_target/widthpf;
    
    preyminspf_mx = ones(Np,1)*preyminspf';
    preymaxspf_mx = ones(Np,1)*preymaxspf';
    
    %pp links
    preyminspp = preymins((Nf+1):N);
    preymaxspp = preymaxs((Nf+1):N);
    preyminspp((pp_diets>0)') = cp_new((pp_diets>0)')-rp_new((pp_diets>0)')/2*widthpp_target/widthpp;
    preymaxspp((pp_diets>0)') = cp_new((pp_diets>0)')+rp_new((pp_diets>0)')/2*widthpp_target/widthpp;
    
    preyminspp_mx = ones(Np,1)*preyminspp';
    preymaxspp_mx = ones(Np,1)*preymaxspp';
    
    preyminsp_mx = [preyminspf_mx preyminspp_mx];
    preymaxsp_mx = [preymaxspf_mx preymaxspp_mx];
    
    preymins_mx((Nf+1):N,:) = preyminsp_mx;
    preymaxs_mx((Nf+1):N,:) = preymaxsp_mx;
    
    web_mx=((n_mx>=preymins_mx)+(n_mx<=preymaxs_mx)==2*ones(N));

end