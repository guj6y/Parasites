richnessAll = [20 30 40];
%connectanceAll = 0.15;
C = .15;
nBasal = 5;
%funcitonalResponsesAll(1,:) = h

nZ = 8;
hAll = [1, 1.2, 2];
ZAll = logspace(-2,5,nZ);
axAll = [.314,.88];
%wij option 1: weak generalist
%wij option 2: strong generalist
                  %S%h%ax%wij
generalistsAll = {'weak','strong'};
%S%h%ax%wij
models = fullfact([3,3,2,2]);
%models = fullfact([1,1,1,1]);
nModels = length(models);

SList = 1:S;

%The abundance level web directory
dataDir =...
    sprintf('/Volumes/Buster/Research/Parasitism/Dynamics/Brose2006Webs');



codeDir = pwd;

nWeb = 100;

persistenceAll1 = cell(nModels,1);
persistenceAll2= cell(nModels,1);
meansAll = cell(2,1);
meansALl{1} = zeros(nModels,8);
meansALl{2} = zeros(nModels,8);


count = 0;
for model = models';
    S = richnessAll(model(1));
    h = hAll(model(2));
    ax = axAll(model(3));
    generalistsAre = generalistsAll{model(4)};
    
count = count+1;
resultsDir1 = sprintf('%s/Results%u%u%u%u',dataDir,model');
resultsDir2 = sprintf('%s/Results2%u%u%u%u',dataDir,model');
persistenceAll1{count} = zeros(nWeb,nZ);
persistenceAll2{count} = zeros(nWeb,nZ);

dataDirS = sprintf('%s/S%03u',S);

%The abundance,connectance level web directory
dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));

for ii = 1:nWeb
    
    finals1 = csvread(sprintf('%s/finals%04u.csv',resultsDir1,ii-1));
    means1 = csvread(sprintf('%s/means%04u.csv',resultsDir1,ii-1));
    stds1 = csvread(sprintf('%s/stds%04u.csv',resultsDir1,ii-1));
    
    finals2 = csvread(sprintf('%s/finals%04u.csv',resultsDir2,ii-1));
    means2 = csvread(sprintf('%s/means%04u.csv',resultsDir2,ii-1));
    stds2= csvread(sprintf('%s/stds%04u.csv',resultsDir2,ii-1));
    
    persistenceAll1{count}(ii,:) = mean(finals1>1e-30);
    persistenceAll2{count}(ii,:) = mean(finals2>1e-30);
    
    
end

meansAll{1}(count,:) = mean(persistenceAll1{count});
meansAll{2}(count,:) = mean(persistenceAll2{count});

end
FRs = false(nModels,3);
Ss = false(nModels,3);
AXs = false(nModels,2);
Ws = false(nModels,2);

FRs(:,1) = models(:,2) == 1;
FRs(:,2) = models(:,2) == 2;
FRs(:,3) = models(:,2) == 3;
Ss(:,1) = models(:,1) == 1;
Ss(:,2) = models(:,1) == 2;
Ss(:,3) = models(:,1) == 3;
AXs(:,1) = models(:,3) == 1;
AXs(:,2) = models(:,3) == 2;
Ws(:,1) = models(:,4) == 1;
Ws(:,2) = models(:,4) == 2;



AXW = fullfact([2,2]);
FRColors = {'r','g','b'};
FRMarks = {'^-','v-'};
legendEntriesH = {'h=1','h=1.2','h=2'};
legendEntriesE = {'; ph',': otf'};
legendEntries = [strcat(legendEntriesH,legendEntriesE{1}) strcat(legendEntriesH,legendEntriesE{2})];
for sChoice = 1:3
    fig = figure;
    
    for subplotChoice = 1:4;
        subplot(2,2,subplotChoice);
        hold on
        for ModelChoice = 1:2;
            plotStyles = strcat(FRColors,FRMarks(ModelChoice));
            for FRChoice = 1:3;
                pickYs = Ss(:,sChoice)&AXs(:,AXW(subplotChoice,1))&...
                    Ws(:,AXW(subplotChoice,2))&FRs(:,FRChoice);
                y1 = meansAll{ModelChoice}(pickYs,:);
                plot(log10(ZAll),y1,plotStyles{FRChoice})
                axis([-2 5 0 1.05]);
                grid on
            end
        end
        
        hold off
    end
    subplot(2,2,1);
    title(sprintf('Invertebrates; S = %u',richnessAll(sChoice)));
    ylabel('Fraction Persistence')
    legend(legendEntries,'Location','SouthEast')
    subplot(2,2,2);
    title(sprintf('Vertebrates; S = %u',richnessAll(sChoice)));
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Weak Generalist','FontWeight','bold')
    subplot(2,2,3);
    ylabel('Fraction Persistence')
    xlabel('Log_{10}(Z)')
    subplot(2,2,4);
    xlabel('Log_{10}(Z)')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Strong Generalist','FontWeight','bold')
    saveas(fig,sprintf('BroseSplit%u',sChoice),'pdf')
end

fig = figure;
subplot(2,2,1);
hold on
y = mean(meansAll{2}(FRs(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(meansAll{2}(FRs(:,2),:));
plot(log10(ZAll),y,'bs--')
y = mean(meansAll{2}(FRs(:,3),:));
plot(log10(ZAll),y,'kv-.')
legend('h=1','h=1.2','h=2','Location','SouthEast')
title('Functional Responses')
ylabel('Fraction Persistence')
axis([-2 5 0 1.05]);
grid on
hold off
subplot(2,2,2);
hold on
y = mean(meansAll{2}(Ss(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(meansAll{2}(Ss(:,2),:));
plot(log10(ZAll),y,'bs--')
y = mean(meansAll{2}(Ss(:,3),:));
plot(log10(ZAll),y,'kv-.')
legend('S=20','S=30','S=40','Location','SouthEast')
title('Diversity')
axis([-2 5 0 1.05]);
grid on
hold off
subplot(2,2,3);
hold on
y = mean(meansAll{2}(AXs(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(meansAll{2}(AXs(:,2),:));
plot(log10(ZAll),y,'bs--')
legend('Invertebrates','Vertebrates','Location','SouthEast')
title('Metabolic Types')
ylabel('Fraction Persistence')
xlabel('log(Z)')
axis([-2 5 0 1.05]);
grid on
hold off
subplot(2,2,4);
hold on
y = mean(meansAll{2}(Ws(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(meansAll{2}(Ws(:,2),:));
plot(log10(ZAll),y,'bs--')
legend('Weak','Strong','Location','SouthEast')
title('Preference Model')
xlabel('log(Z)')
axis([-2 5 0 1.05]);
grid on
hold off
saveas(fig,'BroseAll','pdf')
