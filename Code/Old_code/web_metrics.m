function [property_matrix] = web_metrics(res,con)
    
    %Generate adjacency matrix
    %Note that this is in mathematician's convention: 
    %       
    % Aij = { 1 if j eats i (directed link from i to j)
    %       { 0 if no trophic interaction
    A = sparse(res,con,1);
    
    %Get number of species.
    S = max(max([a b]));
    
    %Get number of links.
    L = length(res);
    
   
    
    


end