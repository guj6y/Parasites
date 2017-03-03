function [l, u] = calculateConfidenceIntervals(x,param)

    n=sum(isfinite(x));
    
    xbar = mean(x,'omitnan');
    sd = std(x,'omitnan');
    
    %calculate 1-alpha confidence intervals.
    tcrit = tinv(1-param.alpha/2,n-1);
    
    l = sd./sqrt(n).*tcrit;
    u = sd./sqrt(n).*tcrit;
    


end