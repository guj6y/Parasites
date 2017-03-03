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



%The abundance level web directory
dataDir =...
    sprintf('/Volumes/Buster/Research/Parasitism/Dynamics/Brose2006Webs');



codeDir = pwd;

nWeb = 100;

persistenceAll = cell(nModels,1);
stabilitiesAll = cell(nModels,1);

meansAll = zeros(nModels,8);
stabsAll = zeros(nModels,8);

count = 0;
for model = models';
    S = richnessAll(model(1));
    h = hAll(model(2));
    ax = axAll(model(3));
    generalistsAre = generalistsAll{model(4)};
    
count = count+1;
resultsDir = sprintf('%s/Results%u%u%u%u',dataDir,model');
persistenceAll{count} = zeros(nWeb,nZ);
stabilitiesAll{count} = zeros(nWeb,nZ);

dataDirS = sprintf('%s/S%03u',S);

%The abundance,connectance level web directory
dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));

for ii = 1:nWeb
    
    finals1 = csvread(sprintf('%s/finals%04u.csv',resultsDir,ii-1));
    means1 = csvread(sprintf('%s/means%04u.csv',resultsDir,ii-1));
    stds1 = csvread(sprintf('%s/stds%04u.csv',resultsDir,ii-1));
    
    CV = stds1./means1;
    CV(finals1<1e-30) = nan;
    stabilitiesAll{count}(ii,:) = mean(CV,'omitnan');
    persistenceAll{count}(ii,:) = mean(finals1>1e-30);
    
    
end

meansAll(count,:) = mean(persistenceAll{count},'omitnan');
stabsAll(count,:) = -mean(stabilitiesAll{count},'omitnan');

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
legendEntries = {'h=1','h=1.2','h=2'};
for sChoice = 1:3
    fig = figure;
    
    for subplotChoice = 1:4;
        subplot(2,2,subplotChoice);
        hold on
        for ModelChoice = 1:2
            plotStyles = strcat(FRColors,FRMarks(ModelChoice));
            for FRChoice = 1:3;
                pickYs = Ss(:,sChoice)&AXs(:,AXW(subplotChoice,1))&...
                    Ws(:,AXW(subplotChoice,2))&FRs(:,FRChoice);
                y1 = meansAll(pickYs,:);
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
    legend(legendEntries,'Location','Best')
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
    saveas(fig,sprintf('BroseSplitPer%u',sChoice),'jpg')
end

for sChoice = 1:3
    fig = figure;
    
    for subplotChoice = 1:4;
        subplot(2,2,subplotChoice);
        hold on
        for ModelChoice = 1:2
            plotStyles = strcat(FRColors,FRMarks(ModelChoice));
            for FRChoice = 1:3;
                pickYs = Ss(:,sChoice)&AXs(:,AXW(subplotChoice,1))&...
                    Ws(:,AXW(subplotChoice,2))&FRs(:,FRChoice);
                y1 = stabsAll(pickYs,:);
                plot(log10(ZAll),y1,plotStyles{FRChoice})
                grid on
            end
        end
        hold off
    end
    subplot(2,2,1);
    title(sprintf('Invertebrates; S = %u',richnessAll(sChoice)));
    ylabel('Population Stabilitys')
    legend(legendEntries,'Location','Best')
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
    ylabel('Populatio Stability')
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
    saveas(fig,sprintf('BroseSplitStab%u',sChoice),'jpg')
end

fig = figure;
subplot(2,2,1);
hold on
y = mean(meansAll(FRs(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(meansAll(FRs(:,2),:));
plot(log10(ZAll),y,'bs--')
y = mean(meansAll(FRs(:,3),:));
plot(log10(ZAll),y,'kv-.')
legend('h=1','h=1.2','h=2','Location','Best')
title('Functional Responses')
ylabel('Fraction Persistence')
axis([-2 5 0 1.05]);
grid on
hold off
subplot(2,2,2);
hold on
y = mean(meansAll(Ss(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(meansAll(Ss(:,2),:));
plot(log10(ZAll),y,'bs--')
y = mean(meansAll(Ss(:,3),:));
plot(log10(ZAll),y,'kv-.')
legend('S=20','S=30','S=40','Location','Best')
title('Diversity')
axis([-2 5 0 1.05]);
grid on
hold off
subplot(2,2,3);
hold on
y = mean(meansAll(AXs(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(meansAll(AXs(:,2),:));
plot(log10(ZAll),y,'bs--')
legend('Invertebrates','Vertebrates','Location','Best')
title('Metabolic Types')
ylabel('Fraction Persistence')
xlabel('log(Z)')
axis([-2 5 0 1.05]);
grid on
hold off
subplot(2,2,4);
hold on
y = mean(meansAll(Ws(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(meansAll(Ws(:,2),:));
plot(log10(ZAll),y,'bs--')
legend('Weak','Strong','Location','Best')
title('Preference Model')
xlabel('log(Z)')
axis([-2 5 0 1.05]);
grid on
hold off
saveas(fig,'BroseAllPer2','jpg')

fig = figure;
subplot(2,2,1);
hold on
y = mean(stabsAll(FRs(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(stabsAll(FRs(:,2),:));
plot(log10(ZAll),y,'bs--')
y = mean(stabsAll(FRs(:,3),:));
plot(log10(ZAll),y,'kv-.')
legend('h=1','h=1.2','h=2','Location','Best')
title('Functional Responses')
ylabel('Population Stability')

grid on
hold off
subplot(2,2,2);
hold on
y = mean(stabsAll(Ss(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(stabsAll(Ss(:,2),:));
plot(log10(ZAll),y,'bs--')
y = mean(stabsAll(Ss(:,3),:));
plot(log10(ZAll),y,'kv-.')
legend('S=20','S=30','S=40','Location','Best')
title('Diversity')

grid on
hold off
subplot(2,2,3);
hold on
y = mean(stabsAll(AXs(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(stabsAll(AXs(:,2),:));
plot(log10(ZAll),y,'bs--')
legend('Invertebrates','Vertebrates','Location','Best')
title('Metabolic Types')
ylabel('Population Stability')
xlabel('log(Z)')

grid on
hold off
subplot(2,2,4);
hold on
y = mean(stabsAll(Ws(:,1),:));
plot(log10(ZAll),y,'ro-')
y = mean(stabsAll(Ws(:,2),:));
plot(log10(ZAll),y,'bs--')
legend('Weak','Strong','Location','Best')
title('Preference Model')
xlabel('log(Z)')

grid on
hold off
saveas(fig,'BroseAllStab2','jpg')

