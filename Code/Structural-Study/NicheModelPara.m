function [res,con,n_new,c_new,r_new,C_all,parasites]= ...
     NicheModelPara(N, C, Np, para_above)                
%------------------------------------------------------------------

%globalStream = RandStream.getGlobalStream;
%reset(globalStream);

tries=10000;                                                               
%error on the connectances
error_n=0.05;

%validation of the plain niche web:
ok_n=0; 

%Randomly assign parasites:

para = datasample(1:N,Np,'Replace',false);
parasites = zeros(N,1);
parasites(para) = true;
parasites = parasites>0;
freeLivers = ~parasites>0; 
Nf = N - Np;

while ((N>0) && (tries>0) && (ok_n==0))
    tries=tries-1;
    %assign niche values from a uniform distribution
    n = rand(N,1);  
    
    
    %designate range for each species

    %parameters for beta distribution:
    alpha = 1;
    beta = (1-2*C)/(2*C); 
    r = betarnd(alpha,beta,N,1); 
    
    %vector of ranges: 
    r=r.*n;
    if para_above
    r(parasites) = r(parasites).*(1-n(parasites));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set center of range, uniformly distributed in [r_i/2,n_i]; 
    c = min(1-r./2,rand(N,1).*(n-r./2)+r./2);
    
    if para_above
        c(parasites) = max(rand(Np,1).*(1-n(parasites)-r(parasites)/2) ...
            + n(parasites),r(parasites)/2);
    end
    %sort everything:
    [n_new, Indx] = sort(n);                                                
    %n_new: niche values in ascending order
    %Indx: indices of species in descending order 
    %(-> 1 is the index of the smallest niche range, 10 is the index of the
    %largest)

    %the smallest r to highest index species 
    r_new = r(Indx);                                                       %NK: Maybe not how I'd do this
    c_new = c(Indx); 

    r_new(1) = 0; %change the r of highest index species to 0
    %so we have a basal species in every web
    
    parasites = parasites(Indx);
    freeLivers = freeLivers(Indx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    
    %lower border of niche range for every prey:
    preymins = c_new - r_new/2;
    
    %upper border of niche range for every predator:
    preymaxs = c_new + r_new/2;

    %fills the empty matrix with niche ranges:
    n_mx = n_new*ones(1,N); 
    %matrix with the lowest points of ranges in every column:
    preymins_mx = ones(N,1)*preymins'; 
    %same, with highest:
    preymaxs_mx = ones(N,1)*preymaxs';

    %Construct the web adjacency matrix;
    %if species in the row is in the diet range of the species in the 
    %column, it gets eaten (if species i is eaten by species j, set (i,j) =
    %1).
    
    web_mx=(n_mx>=preymins_mx)&(n_mx<=preymaxs_mx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  if there is an isolated species                                    %
    %  or something not connected to a basal                               %
    %  or a disconnected group                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %Make sure the graph is weakly connected (no isolated species)
    weak_comp = graphconncomp(sparse(web_mx),'Directed',1,'Weak',1);
    web_connected = (weak_comp == 1);
    
    %make sure that all species have a directed connection from a basal
    %species; only if connected.
   
     %make sure that all species have a directed connection from a basal
    %species; only if connected.
    if web_connected
        unconnected_species = (1:N)';
        basal = sum(web_mx)==0;
        basal_species = unconnected_species(basal);
        %if a parasite is a basal species, it has no host; no good
        if sum(basal(parasites))>0
            web_connected = 0;
        elseif sum(sum(web_mx(parasites,parasites))==sum(web_mx(:,parasites)))>0
            %Parasites have no free-living prey; they ain't parasites.
            web_connected = 0;
        else
            %All parasites have a (non-parasitic) host; check that all
            %energy flow starts at a basal
            connected_species = [];
            
            for kk = basal_species'
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
    
    
    
    if web_connected
        % indices of links (which element in web_mx?)
        links = sum(sum(web_mx));  
        % Actual connectance
        C_web = links/(N^2);  
        if (abs(C_web-C)*1.0/C) > error_n
            ok_n=0;
        else
            ok_n=1;
        end
    end  

end


%Calculate Cff:
linksff = sum(sum(web_mx(freeLivers,freeLivers)));
Cff = linksff/(Nf^2);

%Calculate Cfp:
linksfp = sum(sum(web_mx(freeLivers,parasites)));
Cfp = linksfp/(Np*Nf);

%Calculate Cpf:
linkspf = sum(sum(web_mx(parasites,freeLivers)));
Cpf = linkspf/(Np*Nf);

%Calculate Cpp:
linkspp = sum(sum(web_mx(parasites,parasites)));
Cpp = linkspp/(Np*Np);

links = sum(sum(web_mx));
C_web = links/N^2;

C_all = [C_web,Cff,Cpf,Cfp,Cpp];

    [res,con] = find(web_mx);
tries
    
end