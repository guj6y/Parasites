function plotParasiteExperimentSolOutput(sol,opts,para,basal)
    
    %This function plots a solution structure that is returned by the
    %function integrateParasitesExperiment.m  opts is plotting options.
    %
    % More sophisticated option stuff to come.
    figure;
    hold on
    n = sol.n;
    
    for ii = 1:n
        ii
        soln = sol.(sprintf('sol%u',ii));
        plot(soln.x',log10(soln.y(basal,:))','g')
        plot(soln.x',log10(soln.y(~(basal|para)',:)),'b')
        if sum(para)>0
        plot(soln.x',log10(soln.y(para,:)),'r')
        end
        
    end
    
    hold off
end
