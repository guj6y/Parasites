function [Res,Cons,n,c,r,C_all] = ...
    AddParasites_test(n,c,r,Np,C5,np_min,np_max,para_above,match_fp)
%Add parasites to a plain niche web.  allows for many unique scenarios:
% np_min,np_max: define the interval in which parasites live. para_above:
% true if parasites eat above niche value. match_fp: true if match Cfp over
% C.

Nf = length(n);

    %Ensure that the smallest parasitic niche value is above or below the
    % exterme free liver niche value.  This is necessary so that all the
    % parasites can have a host.  This will change the distribution of
    % parasitic niche values - we can expect the new range of parasite
    % values to be in U(n_(1),1) or U(0,n_(Nf)), where we can easily
    % calculate the expected value of the order statistic.  Does this
    % change the mean of n_p and do we need to take that into account when
    % calculating the range for each parasite?

%Ensuring every parasite will have a host.
np_min = max(np_min,min(n));
%E[min(n)] = 1/(Nf+1);
np_max = min(np_max,max(n));
%E[max(n)] = Nf/(Nf+1);

Enp = (np_min+np_max)/2;

%Calculate expected width of diet


C = C5(1);
Cff = C5(2);
Cpf = C5(3); %Not used.
Cfp = C5(4);
Cpp = C5(5); %Not used.


%Defining the expected width of parasitic ranges.  Adjusting for the fact
%that each parasite is guaranteed to have a host.

if match_fp %If we intend to match f->p connectance...
    C_new = (Cfp*Nf*Np - Np)/(Nf*Np);  %Take Cfp as given by the empirical web
else %if we intend to match overall connectance...
    %Define what the needs to be to get the right number of links.
    C_new = (C*(Np+Nf)^2 - Cff*(Np*Nf+Nf^2) -Np)/(Np^2 + Np*Nf);
end

tries=275;  %Low number of tries to get a higher turnover on the base webs.

%error on the connectance we are trying to hit; free-livers are set
error_n=0.05;

num_species = Np+Nf;

%validation of the niche web:
ok_n=0;



freeLivers = true(num_species,1);
freeLivers((Nf+1):num_species) = false;
parasites = ~freeLivers;

Sp = 1:num_species;

cp = zeros(Np,1);
np = zeros(Np,1);
rp = zeros(Np,1);
while ((num_species>0) && (tries>0) && (ok_n==0))
    
    tries=tries-1;
    %assign parasitic niche values from a uniform distribution
    np = rand(Np,1)*(np_max-np_min)+np_min;
   
    %-----designate feeding range for each parasite
    
    %parameters for beta distribution:
    alpha = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if para_above ==1
        % set center of range, uniformly distributed in [n_p,1-r_p/2];
        beta = alpha*(1-Enp)/C_new - alpha;
        rp = betarnd(alpha,beta,Np,1);
        rp = rp.*(1-np);
        
        %For each parasite, pick a free-living niche value lower than
        %itself as the center of its range.  Makes sense since parasites 
        %should usually be pretty highly specialized to their hosts?  If
        %not, could also wiggle a bit to more closely mimic the niche model
        
        for jj = 1:Np
            cp(jj) = datasample(n(n>np(jj)),1);
        end
        
    else
        % Parasites eat below
        beta = alpha*(Enp)/C_new - alpha;
        rp = betarnd(alpha,beta,Np,1);
        rp = rp.*np;

        mean(rp)
        mean(r)
        for jj = 1:Np
            cp(jj) = datasample(n(n<np(jj)),1);
        end
    end
    
    %sort everything:
    [np_new, Indx] = sort(np);
    %np_new: niche values in ascending order Indx: indices of species in
    %descending order
    
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
    
    %Construct the web adjacency matrix; if species in the row is in the
    %diet range of the species in the column, it gets eaten (if species i
    %is eaten by species j, set (i,j) = 1).
    
    web_mx=((n_mx>=preymins_mx)&(n_mx<=preymaxs_mx));
    
    %Need to check that: 1. all parasites have a host (added parasite isn't
    %a new basal species) 2. web is connected
    
    %Make sure that all species are weakly connected to the specified basal
    %species Make sure the graph is weakly connected (no isolated species)
    
    weak_comp = graphconncomp(sparse(web_mx),'Directed',1,'Weak',1);
    web_connected = (weak_comp == 1);
    
    %make sure that all species have a directed connection from a basal
    %species; only if connected.
    if web_connected
        unconnected_species = (1:num_species)';
        basal_species = sum(web_mx)==0;
        
        %if a parasite is a basal species, it has no host; no good
        if sum(basal_species>Nf)>0
            web_connected = 0;
        elseif sum(sum(web_mx(parasites,parasites))==sum(web_mx(:,parasites)))>0
            %then... Parasites have no free-living prey; they ain't
            %parasites.
            web_connected = 0;
        else
            %All parasites have a (non-parasitic) host; check that all
            %energy flow starts at a basal
            connected_species = [];
            
            for kk = Sp(basal_species)
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
    
    if web_connected %Check the connectance if the web is connected.
        if match_fp
            %Links we were trying to match
            links = sum(sum(web_mx(1:Nf,(Nf+1):num_species)));
            % Actual connectance
            C_web = links/(Nf*Np);
            %Within Tolerance?
            if (abs(C_web-C_new)/C_new) > error_n
                ok_n=0;
            else
                ok_n=1;
            end
        else
            %Links that we were trying to match
            links = sum(sum(web_mx));
            % Actual connectance
            C_web = links/(num_species^2);
            %Within Tolerance?
            if (abs(C_web-C_new)/C_new) > error_n
                ok_n=0;
            else
                ok_n=1;
            end
        end
        
    end
end



if tries == 0
    [~,~,n,c,r] = NicheModel_nk(Nf,Cff);
    [Res,Cons,n,c,r,C_all] = ...
        AddParasites(n,c,r,Np,C5,np_min,np_max,para_above,match_fp);
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