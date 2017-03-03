function [counts,fractions] = calculateParaStats(order,para,free)

    S = numel(order);
    
    paraOrders = order(para);
    freeOrders = order(free);
    
    nFree = numel(freeOrders);
    nPara = numel(paraOrders);
    
    counts = zeros(S,1);
    fractions = zeros(S,1);
    
    for ii = 1:S
        if para(ii)>0
            
            counts(ii) = sum(freeOrders<order(ii));
            
        elseif free(ii)>0
            
            counts(ii) = sum(paraOrders<order(ii));
            
        else
            
            counts(ii) = nan;
            
        end
        
        
    end
    
    fractions(para) = counts(para)/nFree;
    fractions(free) = counts(free)/nPara;
    
    fractions(~(para|free)) = nan;

end