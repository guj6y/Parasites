h = figure;
plotStyle = {'ro-','bx-','gs-','k*-','m+-'};
for ii = 1:5
    hold on
    plot(log10(ZAll),mean(persistenceAll{ii}),plotStyle{ii})
    
    
    
end
hold off

title('Persistence for Different Final Times')
xlabel('log_{10}(Z)')
ylabel('Fraction persistence')
legend('Tf = 2,000',...
    'Tf = 5,000',...
    'Tf = 10,000',...
    'Tf = 20,000',...
    'Tf = 40,000',...
    'Location','SouthEast')

saveas(h,'TfStudy.pdf')