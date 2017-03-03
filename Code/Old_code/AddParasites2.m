function [Res,Cons,n,c,r,C_all] = ...
    AddParasites2(n,c,r,Np,C2,np_min,np_max,para_above)
%TODO: add limit to the recursion

%Second iteration.  Match the overall connectance of web with parasites.  
%Restrict parasitic niche values.  Allow all types of 
%interactions.




%np_max = 0.3;
Enp = (np_min+np_max)/2;

%Calculate expected width of diet 
Nf = length(n);
Cff = C2(1);
C = C2(2);



C_new = (C*(Np+Nf)^2 - Cff*(Np*Nf+Nf^2))/(Np^2 + Np*Nf);

tries=275;                                                               
%error on the connectance for free-livers
error_n=0.05;
%error on the connectance for parasites
num_species = Np+Nf;

%validation of the plain niche web:
ok_n=0; 
while ((num_species>0) && (tries>0) && (ok_n==0))
    


%while ((num_species>0) && (tries>0) && (ok_n==0))
    tries=tries-1;
    %assign parasitic niche values from a uniform distribution
    np = rand(Np,1)*(np_max-np_min)+np_min;  
    
    
    %designate range for each parasite

    %parameters for beta distribution:
    alpha = 1;
    beta = (1-Enp-C_new)/(C_new); 
    rp = betarnd(alpha,beta,Np,1); 
    
    %vector of ranges:
    rp = rp.*(1-np);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if para_above ==1
        % set center of range, uniformly distributed in [n_p,1-r_p/2]; 
        cp = rand(Np,1).*(1-rp/2-np)+np;
    else
        % No restriction on parasitic feeding centers
        cp = rand(Np,1).*(1-rp)+rp/2;
    end

    %sort everything:
    [np_new, Indx] = sort(np);                                                
    %np_new: niche values in ascending order
    %Indx: indices of species in descending order 

    rp_new = rp(Indx);                                                       %NK: Maybe not how I'd do this
    cp_new = cp(Indx); 

    %Put parasites in with free-livers
    n_new = [n;np_new];
    c_new = [c;cp_new];
    r_new = [r;rp_new];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %lower border of niche range for every prey:
    preymins = c_new - r_new/2;
    
    %upper border of niche range for every predator:
    preymaxs = c_new + r_new/2;

    %fills the empty matrix with niche ranges:
    n_mx = n_new*ones(1,num_species); 
    %matrix with the lowest points of ranges in every column:
    preymins_mx = ones(num_species,1)*preymins'; 
    %same, with highest:
    preymaxs_mx = ones(num_species,1)*preymaxs';

    %Construct the web adjacency matrix;
    %if species in the row is in the diet range of the species in the 
    %column, it gets eaten (if species i is eaten by species j, set (i,j) =
    %1).
    
    web_mx=((n_mx>=preymins_mx)+(n_mx<=preymaxs_mx)==2*ones(num_species));
    
    %Need to check that:
        %1. all parasites have a host (added parasite isn't a new basal
        %species)
        %2. web is connected
    
    %Make sure that all species are weakly connected to the specified basal
    %species
    unconnected_species = (1:num_species)';
    
    connected_species = [];
    test_mx = (web_mx+web_mx')>0;
    
    connected_species = walk(1,connected_species,test_mx);
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
        elseif sum(...
                sum(web_mx(Nf+1:end,Nf+1:end))>0 &... %hyperparasites have
                (sum(web_mx(:,Nf+1:end)) ... equal diet size on parasites
            == sum(web_mx(Nf+1:end,Nf+1:end))))...  as total diet size
             %then...
        %Parasites have no free-living prey; they ain't parasites.
        web_connected = 0;
        else
            %All parasites have a host; check that all energy flow starts
            %at a basal
            connected_species = [];

            for kk = basal_species
                connected_species = walk(kk,connected_species,web_mx);
            end

            unconnected_species(connected_species) = [];

            if isempty(unconnected_species)
                web_connected = 1;
            else
                web_connected = 0;
            end
        end
    end    

    %TODO: make sure both connectances are hit?    
    if web_connected
        % indices of links (which element in web_mx?)
        links = sum(sum(web_mx));  
        % Actual connectance
        C_web = links/(num_species^2);
        if (abs(C_web-C)/C) > error_n
            ok_n=0;
        else
            ok_n=1;
        end
    end
end
    


if tries == 0
    [~,~,n,c,r,~] = NicheModelPara(Nf,Cff,Np);
    [Res,Cons,n,c,r,C_all] = ...
        AddParasites2(n,c,r,Np,C2,np_min,np_max,para_above);
else
%Calculate Cff:
linksff = sum(sum(web_mx(1:Nf,1:Nf)));
Cff = linksff/(Nf^2);

%Calculate Cfp:
linksfp = sum(sum(web_mx(1:Nf,(Nf+1):(Np+Nf))));
Cfp = linksfp/(Np*Nf);

%Calculate Cpf:
linkspf = sum(sum(web_mx((Nf+1):(Np+Nf),1:Nf)));
Cpf = linkspf/(Np*Nf);

%Calculate Cpp:
linkspp = sum(sum(web_mx((Nf+1):(Np+Nf),(Nf+1):(Np+Nf))));
Cpp = linkspp/(Np*Np);

C_all = [C_web,Cff,Cpf,Cfp,Cpp];    
    [Res,Cons] = find(web_mx);
    n = n_new;
    c = c_new;
    r = r_new;
end

end