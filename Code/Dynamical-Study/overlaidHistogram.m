function [figHandle] = overlaidHistogram(prop,Biomass,varargin)
    
    [~,Edges] = histcounts(prop,varargin{:});
    figHandle = figure;
    hold on
    H1 = histogram(prop(Biomass>0),Edges);
    H1.FaceColor = 'b';
    
    H2 = histogram(prop(Biomass==0),Edges);
    H2.FaceColor = 'r';
    

end
