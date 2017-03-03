function [handle] = plot3Properties(prop1,prop2,prop3,biomass)
handle = figure;
hold on
plot3(prop1(biomass==0),prop2(biomass==0),prop3(biomass==0),'r.')

plot3(prop1(biomass>0),prop2(biomass>0),prop3(biomass>0),'b.')

end