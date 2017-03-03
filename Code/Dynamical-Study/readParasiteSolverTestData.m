
S = 40;
C = .15;

close all

nFPar = 8;
fParAll = linspace(.05,.75,nFPar);



models = fullfact([2,2,2,3]);

%This changes the order of experiments; species types is now the last thing
%tested so that I can get more intersting factors out of the way first..
%ZFree,ZPara,FracFree,Concomittant
models = [models(:,4),models(:,1:3)];

nModels = length(models);

SList = 1:S;

%The abundance level web directory
dataDir = sprintf('~/Desktop/JulySimulations/Web1SolverTests');



codeDir = pwd;
startWeb = 1;
endWeb = 1;
nWeb = endWeb-startWeb+1;
persistenceAll = cell(nModels,1);

solverDirs = {'ode45',...
    'ode15s',...
    'ode23',...
    'ode113'};

nSolver = length(solverDirs);

persistenceAll = cell(nSolver,1);
modelCount = 0;
for solver = 1:nSolver
    
    dataDirSolver = sprintf('%s/%s',dataDir,solverDirs{solver});
    persistenceAll{solver} = zeros(nfPar,nModels);
    totBiomAssAllAll{solver} = zeros(nfPar,nModels);
    meanCoefVarAll{solver} = zeros(nfPar,nModels);
    modelCount = 0;
    for model = models';
        modelCount = modelCount+1;
        
        resultsDir = sprintf('%s/Results%u%u%u%u',dataDirSolver,model');
        
        dataDirS = sprintf('%s/S%03u',S);
        
        %The abundance,connectance level web directory
        dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));
        
        for ii = 1:nWeb
            
            finals1 = csvread(sprintf('%s/finals%04u.csv',resultsDir,ii-1));
            means1 = csvread(sprintf('%s/means%04u.csv',resultsDir,ii-1));
            stds1 = csvread(sprintf('%s/stds%04u.csv',resultsDir,ii-1));
            
            persistenceAll{solver}(:,modelCount) = mean(finals1(:,1:end-1)>1e-30);
            totBiomassAll{solver}(:,modelCount) = sum(finals1(:,1:end-1)>1e-30);
            meanCoefVarAll{solver}(:,modelCount) = mean(stds1(:,1:end-1)./means1(:,1:end-1),'omitnan')';
        end
        
        
    end
    
    
end

ZFrees = false(nModels,3);
ZParas = false(nModels,2);
freeLivings = false(nModels,2);
concs = false(nModels,2);

ZFrees(:,1) = models(:,1) == 1;
ZFrees(:,2) = models(:,1) == 2;
ZFrees(:,3) = models(:,1) == 3;
ZParas(:,1) = models(:,2) == 1;
ZParas(:,2) = models(:,2) == 2;
freeLivings(:,1) = models(:,3) == 1;
freeLivings(:,2) = models(:,3) == 2;
concs(:,1) = models(:,4) == 1;
concs(:,2) = models(:,4) == 2;

solverDirs = {'ode45',...
    'ode15s',...
    'ode23',...
    'ode113'};

figTitles = {'10^{-3}';'10^{-4}'};
AXW = fullfact([2,2]);
ZColors = {'r','k','b'};
modelMarks = {'^-','v-','+-','x-'};
legendEntries = cell(12,1);
modelLegendEntries = solverDirs;
zLegendEntries = {'Z= 10|','Z = 100|','Z = 1,000|'};
for ZParaChoice  = 1:2
fig = figure;
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    count = 0;
    for ZFreeChoice = 1:3;
        for modelChoice = 1:4
            count = count+1;
            plotStyles = strcat(ZColors(ZFreeChoice),modelMarks(modelChoice));
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            y1 = persistenceAll{modelChoice}(:,pickYs');
            plot(fParAll,y1,plotStyles{:})
            %axis([0 1 0 1.05]);
            grid on
            legendEntries(count) = strcat(zLegendEntries(ZFreeChoice),...
                                        modelLegendEntries(modelChoice));
        end
    end
end
    hold off
    subplot(2,2,1);
    title(sprintf('%s\nNo Fraction Free Living Stage',figTitles{ZParaChoice(:)}));
    ylabel('Fraction Persistence')
    leg = legend(legendEntries,'Location','SouthEast');
    leg.FontSize = 10;
    subplot(2,2,2);
    title(sprintf('Fraction Free Living Stage'));
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'No Concomittant','FontWeight','bold')
    subplot(2,2,3);
    ylabel('Fraction Persistence')
    xlabel('f_{par}')
    subplot(2,2,4);
    xlabel('f_{par}')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    
end
saveas(fig,sprintf('ExperimentFigure'),'pdf')
%
% fig = figure;
% subplot(2,2,1);
% hold on
% y = mean(meansAll(ZFrees(:,1),:));
% plot(log10(ZAll),y,'ro-')
% y = mean(meansAll(ZFrees(:,2),:));
% plot(log10(ZAll),y,'bs--')
% y = mean(meansAll(ZFrees(:,3),:));
% plot(log10(ZAll),y,'kv-.')
% legend('h=1','h=1.2','h=2','Location','SouthEast')
% title('Functional Responses')
% ylabel('Fraction Persistence')
% axis([-2 5 0 1.05]);
% grid on
% hold off
% subplot(2,2,2);
% hold on
% y = mean(meansAll(ZParas(:,1),:));
% plot(log10(ZAll),y,'ro-')
% y = mean(meansAll(ZParas(:,2),:));
% plot(log10(ZAll),y,'bs--')
% y = mean(meansAll(ZParas(:,3),:));
% plot(log10(ZAll),y,'kv-.')
% legend('S=20','S=30','S=40','Location','SouthEast')
% title('Diversity')
% axis([-2 5 0 1.05]);
% grid on
% hold off
% subplot(2,2,3);
% hold on
% y = mean(meansAll(freeLivings(:,1),:));
% plot(log10(ZAll),y,'ro-')
% y = mean(meansAll(freeLivings(:,2),:));
% plot(log10(ZAll),y,'bs--')
% legend('Invertebrates','Vertebrates','Location','SouthEast')
% title('Metabolic Types')
% ylabel('Fraction Persistence')
% xlabel('log(Z)')
% axis([-2 5 0 1.05]);
% grid on
% hold off
% subplot(2,2,4);
% hold on
% y = mean(meansAll(concs(:,1),:));
% plot(log10(ZAll),y,'ro-')
% y = mean(meansAll(concs(:,2),:));
% plot(log10(ZAll),y,'bs--')
% legend('Weak','Strong','Location','SouthEast')
% title('Preference Model')
% xlabel('log(Z)')
% axis([-2 5 0 1.05]);
% grid on
% hold off
% saveas(fig,'BroseAll2','pdf')
