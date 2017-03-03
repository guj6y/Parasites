function [web_connected,unconnected_species] = checkConnected(web_mx,N)

%Make sure that all species are weakly connected to the specified basal
    %species (no isolated groups)
    unconnected_species = (1:N)';
    connected_species = [];
    
    %need test matrix since this is weakly connected component; symmetrized
    % (undirected) graph
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
        unconnected_species = (1:N)';
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

end