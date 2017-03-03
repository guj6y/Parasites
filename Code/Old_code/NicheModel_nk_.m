function [Res,Cons,n_new,c_new,r_new]= ...
     NicheModel_nk(num_species, C)                
%------------------------------------------------------------------

%globalStream = RandStream.getGlobalStream;
%reset(globalStream);

tries=10000;                                                               
%error on the connectance for free-livers
error_n=100;
%error on the connectance for parasites
nicheweb=0;

c=[]; %center values
r=[]; %range values

%validation of the plain niche web:
ok_n=0; 

while ((num_species>0) && (tries>0) && (ok_n==0))
    tries=tries-1;
    %assign niche values from a uniform distribution
    n = rand(num_species,1);  
    
    
    %designate range for each species

    %parameters for beta distribution:
    alpha = 1;
    beta = (1-2*C)/(2*C); 
    r = betarnd(alpha,beta,num_species,1); 
    
    %vector of ranges:
    r = r.*n;  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set center of range, uniformly distributed in [r_i/2,n_i]; 
    c = min(1-r./2,rand(num_species,1).*(n-r./2)+r./2);

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  if there is an isolated species                                    %
    %  or something not connected to a basal                              %
    %  or a disconnected group                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
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

    web_connected = 1;
    
    if web_connected
        % indices of links (which element in web_mx?)
        links = sum(sum(web_mx));  
        % Actual connectance
        C_web = links/(num_species^2);  
        if (abs(C_web-C)*1.0/C) > error_n
            ok_n=0;
        else
            ok_n=1;
        end
    end  

end
    [Res,Cons] = find(web_mx);

    
end