function [Res,Cons,n,c,r]= ...
     NicheModelParasites(num_species2, C4)
 
 
%globalStream = RandStream.getGlobalStream;
%reset(globalStream);

%third iteration; redraws all species if there is a problem.

np_min = 0;
np_max = 0.3;
Enp = (np_min+np_max)/2;

%Calculate expected width of diet 
Nf = num_species2(1);
Np = num_species2(2);

Cff = C4(1);
Cpf = C4(2);
Cfp = C4(3);
Cpp = C4(4);

C = (Cff*Nf^2 + (Nf*Np)*(Cfp+Cpf) + Np^2*Cpp)/(Np+Nf)^2;



tries=10000;                                                               
%error on the connectance for free-livers
error_n=0.05;
%error on the connectance for parasites
nicheweb=0;

c=[]; %center values
r=[]; %range values

%validation of the plain niche web:
ok_n=0; 

num_species = Nf+Np;
N = num_species;

while ((num_species>0) && (tries>0) && (ok_n==0))
    tries=tries-1;
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
    preymins_a = c-r/2;
    preymaxs_a = c+r/2;
    
    preymins_a_mx = ones(N,1)*preymins_a';
    preymaxs_a_mx = ones(N,1)*preymaxs_a';
    %}
    
    
    
    
    %Creat each block of the matrix:
    %ff links
    preyminsff = cf_new-rf_new/2;
    preymaxsff = cf_new+rf_new/2;
    
    preyminsff_mx = ones(Nf,1)*preyminsff';
    preymaxsff_mx = ones(Nf,1)*preymaxsff';
    
    %pf links
    preyminspf = cf_new-rf_new/2;
    preymaxspf = cf_new+rf_new/2;
    
    preyminspf_mx = ones(Np,1)*preyminspf';
    preymaxspf_mx = ones(Np,1)*preymaxspf';
    
    %fp links
    preyminsfp = cp_new-rp_new/2;
    preymaxsfp = cp_new+rp_new/2;
    
    preyminsfp_mx = ones(Nf,1)*preyminsfp';
    preymaxsfp_mx = ones(Nf,1)*preymaxsfp';
    
    %pp links
    preyminspp = cp_new-rp_new*Cpp/Cfp/2;
    preymaxspp = cp_new+rp_new*Cpp/Cfp/2;
    
    preyminspp_mx = ones(Np,1)*preyminspp';
    preymaxspp_mx = ones(Np,1)*preymaxspp';

    
    %fills the empty matrix with niche ranges:
    n_mx = n*ones(1,num_species); 
    
    preymins_mx = [preyminsff_mx preyminsfp_mx;preyminspf_mx preyminspp_mx];
    preymaxs_mx = [preymaxsff_mx preymaxsfp_mx;preymaxspf_mx preymaxspp_mx];

    %Construct the web adjacency matrix;
    %if species in the row is in the diet range of the species in the 
    %column, it gets eaten (if species i is eaten by species j, set (i,j) =
    %1).
    
    web_mx=((n_mx>=preymins_mx)+(n_mx<=preymaxs_mx)==2*ones(num_species));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  if there is an isolated species                                    %
    %  or something not connected to a basal                              %
    %  or a disconnected group                                            %
    %  ***Need to add uniqueness condition                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Make sure that all species are weakly connected to the specified basal
    %species
    unconnected_species = (1:num_species)';
    
    connected_species = [];
    
    test_mx = sparse((web_mx+web_mx')>0);
    
    [res_t,cons_t] = find(test_mx);
    
    connected_species = walk1(1,connected_species,res_t,cons_t);
    
    unconnected_species(connected_species) = [];
    
    if isempty(unconnected_species)
        web_connected = 1;
    else
        web_connected = 0;
    end
    
    %make sure that all species have a directed connection from a basal
    %species; only if connected.
    if web_connected
        unconnected_species = (1:num_species)';
        basal_species = unconnected_species(sum(web_mx)==0)';
        
        %if a parasite is a basal species, it has no host; no good
        if sum(basal_species>Nf)
            web_connected = 0;
        else
            %All parasites have a host; check that all energy flow starts
            %at a basal
            connected_species = [];
            [res,cons] = find(web_mx);
            for kk = basal_species
                connected_species = walk1(kk,connected_species,res,cons);
            end

            unconnected_species(connected_species) = [];

            if isempty(unconnected_species)
                web_connected = 1;
            else
                web_connected = 0;
            end
        end
    end     
    
    %need to make sure that the free-liver web has the right connectance.
    web_free = web_mx(1:Nf,1:Nf);
    
    if web_connected
        links = sum(sum(web_free));
        C_web_ff = links/Nf^2;
        if (abs(C_web_ff-Cff)*1.0/Cff) > error_n
            ok_n=0;
            fprintf('miss Cff\n')
        else
            ok_n=1;
        end
    end
    
    if web_connected
        links = sum(sum(web_mx(1:Nf,(Nf+1):(Np+Nf))));
        C_web_fp = links/Nf/Np;
        if (abs(C_web_fp-Cfp)*1.0/Cfp) > error_n
            ok_n=0;
            fprintf('miss Cfp\n')
        else
            ok_n=1;
        end
    end
    
    if web_connected
        links = sum(sum(web_mx((Nf+1):(Np+Nf),1:Nf)));
        C_web_pf = links/Nf/Np;
        if (abs(C_web_pf-Cpf)*1.0/Cpf) > error_n
            ok_n=0;
            fprintf('miss Cpf\n')
        else
            ok_n=1;
        end
    end
    
    if web_connected
        links = sum(sum(web_mx((Nf+1):(Np+Nf),(Nf+1):(Np+Nf))));
        C_web_pp = links/Np^2;
        if (abs(C_web_pp-Cpp)*1.0/Cpp) > error_n
            ok_n=0;
            fprintf('miss Cpp\n')
        else
            ok_n=1;
        end
    end
    
    if web_connected
        % indices of links
        links = sum(sum(web_mx));  
        % Actual connectance
        C_web = links/(num_species^2);  
        if (abs(C_web-C)*1.0/C) > error_n
            ok_n=0;
        else
            ok_n=1;
%             Cff
%             C_web_ff
%             Cfp
%             C_web_fp
%             Cpf 
%             C_web_pf
%             Cpp
%             C_web_pp
        end
    end  
    
end
    [Res,Cons] = find(web_mx);
    sum((((.1-r/2)<c)+(c<(.3+r/2)))==2)
end