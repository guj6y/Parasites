S = 40;
C = .15;

close all
%A stupid thing


models = fullfact([2,2,2,2]);

%This changes the order of experiments; species types is now the last thing
% tested so that I can get more intersting factors out of the way first..
%ZFree,ZPara,FracFree,Concomittant
models = [models(:,4),models(:,1:3)];

nModels = length(models);

SList = 1:S;

%The abundance level web directory
dataDir = sprintf('/Volumes/Buster/Research/Parasitism/Dynamics/JulySimulations');

kFree = [2 3];
kPara = [-3,-4];
%
codeDir = pwd;
startWeb = 1;
endWeb = 100;
nWeb = endWeb-startWeb+1;
persistenceAll = cell(nModels,1);
persistenceFree = cell(nModels,1);
persistencePara = cell(nModels,1);

xAll = cell(nModels,nWeb);

avgBiomassAll = cell(nModels,1);
avgBiomassFree = cell(nModels,1);
avgBiomassPara = cell(nModels,1);

extctStatisticParaAll = cell(nModels,1);
extctStatisticFreeAll = cell(nModels,1);

fracConExtctFirstAll = cell(nModels,1);
fracResExtctFirstAll = cell(nModels,1);

fracParaExtctFirstAll = cell(nModels,1);
fracHostExtctFirstAll = cell(nModels,1);

fracConExtctFirstAll_genAveraged = cell(nModels,1);
fracResExtctFirstAll_vulAveraged = cell(nModels,1);

fracParaExtctFirstAll_genAveraged = cell(nModels,1);
fracHostExtctFirstAll_vulAveraged = cell(nModels,1);


nAliveAll = cell(nModels,1);
nAliveFree = cell(nModels,1);
nAlivePara = cell(nModels,1);

count = startWeb-1;

bsrAll = cell(nModels,1);
extctAll = cell(nModels,1);

LLAll = cell(nModels,1);
resExtctAll = cell(nModels,1);
conExtctAll = cell(nModels,1);
groups = cell(nModels,1);
resBodySizeAll = cell(nModels,1);
conBodySizeAll = cell(nModels,1);
vulBasalAll = cell(nModels,1);
meanVulBasalAll = zeros(nWeb,nModels);
sumVulBasalAll = zeros(nWeb,nModels);
nBasal = zeros(nWeb,nModels);


dataDirS = sprintf('%s/S%03u',dataDir,S);

%The abundance,connectance level web directory
dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));

model = models(1,:);
resultsDir = sprintf('%s/ResultsPATL%u%u%u%u',dataDir,model');

finals = csvread(sprintf('%s/finals%04u.csv',resultsDir,0))';
[~,nFPar0] = size(finals);
%Checks for doubled data

fPar0 = unique(sort(finals(1,:)));
if nFPar0 > numel(fPar0)
    
    nFPar0 = numel(fPar0);
end

meansAll = zeros(nModels,nFPar0);
stdsAll = zeros(nModels,nFPar0);

meansAllFree = zeros(nModels,nFPar0);
stdsAllFree = zeros(nModels,nFPar0);

meansAllPara = zeros(nModels,nFPar0);
stdsAllPara = zeros(nModels,nFPar0);

b_All = zeros(nModels,nFPar0);
b_AllFree = zeros(nModels,nFPar0);
b_AllPara = zeros(nModels,nFPar0);



biomassRatiosAll = zeros(0,nFPar0);

paraMxBinAll = zeros(0,nFPar0);
patlAll = [];
basalAll = [];
fParAll = fPar0;
for model = models'
    resultsDir = sprintf('%s/ResultsPATL%u%u%u%u',dataDir,model');
    count = count+1;
    try
    finals = csvread(sprintf('%s/finals%04u.csv',resultsDir,0))';
    [~,nFPar] = size(finals);
    catch
    end
    if nFPar > nFPar0
        
        nFPar = nFPar/2;
        
    end
    
    persistenceAll{count} = zeros(nWeb,nFPar);
    persistenceFree{count} = zeros(nWeb,nFPar);
    persistencePara{count} = zeros(nWeb,nFPar);
    
    fracConExtctFirstAll{count} = zeros(nWeb,nFPar);
    fracResExtctFirstAll{count} = zeros(nWeb,nFPar);
    
    fracParaExtctFirstAll{count} = zeros(nWeb,nFPar);
    fracHostExtctFirstAll{count} = zeros(nWeb,nFPar);
    
    fracConExtctFirstAll_genAveraged{count} = zeros(nWeb,nFPar);
    fracResExtctFirstAll_vulAveraged{count} = zeros(nWeb,nFPar);
    
    dataDirS = sprintf('%s/S%03u',dataDir,S);
    
    %The abundance,connectance level web directory
    dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));
    
    bsrAll{count} = zeros(0,nFPar);
    extctAll{count} = zeros(0,nFPar);
    LLAll{count} = zeros(0,2);
    
    resExtctAll{count} = zeros(0,nFPar);
    conExtctAll{count} = zeros(0,nFPar);
    resBodySizeAll{count} = zeros(0,nFPar);
    conBodySizeAll{count} = zeros(0,nFPar);
    
    extctStatisticParaAll{count} = zeros(nWeb,nFPar);
    extctStatisticFreeAll{count} = zeros(nWeb,nFPar);
    
    
    for ii = 1:nWeb
        
        fid = fopen(sprintf('%s/idxPar%04u.txt',dataDirSC,ii-1),'r');
        idxPar = textscan(fid,'%u\n','CollectOutput',1);
        idxPar = idxPar{1};
        fclose(fid);
        
        LL = csvread(sprintf('%s/web%04u.csv',dataDirSC,ii-1));
        try
            finals = csvread(sprintf('%s/finals%04u.csv',resultsDir,ii-1))';
            means = csvread(sprintf('%s/means%04u.csv',resultsDir,ii-1))';
            stds = csvread(sprintf('%s/stds%04u.csv',resultsDir,ii-1))';
        catch
            fprintf('Couldn''t load something.\nModel: %u%u%u%u\n',model)
            fprintf('Web: %u',ii)
            finals = zeros(nRows,nCols);
            means = zeros(nRows,nCols);
            stds = zeros(nRows,nCols);
        end
        
        [nRows,nCols] = size(finals);
        [fParSorted,fParOrder] = sort(finals(1,:));
        finals = finals(:,fParOrder);
        means= means(:,fParOrder);
        stds = stds(:,fParOrder);
        %This fixes doubles.  Becareful though!  it re-writes data!!
        if nCols > nFPar
            finals = finals(:,1:end/2);
            means = means(:,1:end/2);
            stds = stds(:,1:end/2);
            
            csvwrite(sprintf('%s/finals%04u.csv',resultsDir,ii-1),finals');
            csvwrite(sprintf('%s/means%04u.csv',resultsDir,ii-1),means');
            csvwrite(sprintf('%s/stds%04u.csv',resultsDir,ii-1),stds');
            
        end
        
        %Extract the number and fraction of parasites for each column.
        nPars = finals(2,:);
        fPars = finals(1,:);
        
        paraMxBin = false(S,nFPar);
        freeMxBin = paraMxBin;
        extctOrderAll = zeros(S,nFPar);
        
        
        %Useful information extracted, take only the data (row of zeros at
        %the end! It's a blank entry, ooops!
        finals = finals(3:end-1,:);
        means = means(3:end-1,:);
        stds = stds(3:end-1,:);
       
        
        
        props = calculateLocalProperties(LL(:,1),LL(:,2));
        
        gen = props(:,2)*numel(LL)/2/S;
        vul = props(:,3)*numel(LL)/2/S;
        
        linkDensity = length(LL(:,1))/S;
        
        swtl = props(:,10);
        sptb = props(:,8) + 1;
        patl = 2*swtl - sptb;
        
        basal = props(:,10) == 1;
        
        onlyMeans = means;
        
        for kk = 1:nFPar
            paraMxBin(idxPar(1:nPars(kk)),kk) = true;
            freeMxBin(:,kk) = ~(paraMxBin(:,kk)|basal);
            try
            timeSeriesDir = sprintf('%s/TimeSeries/web%04u',resultsDir,ii-1);
            timeData = csvread(sprintf('%s/fpar%4f.csv',timeSeriesDir,fPars(kk)));
            
            
            times = timeData(1,:);
            biomasses = timeData(2:end,2:end);
            
            extctOrder = sum(biomasses==0,2);
            nExtinct = max(extctOrder);
            extctOrder(extctOrder==0) = nan;
            extctOrder = nExtinct-extctOrder+1;
            
            extctOrderAll(:,kk) = extctOrder;
            
            [counts,fractions] = calculateParaStats(extctOrder,...
                paraMxBin(:,kk),freeMxBin(:,kk));
            
            end
        end
        
        
        extctOrderPara = extctOrderAll;
        extctOrderPara(~paraMxBin) = nan;
        
        extctOrderFree = extctOrderAll;
        extctOrderFree(~freeMxBin) = nan;
        
        extctStatisticPara = mean(extctOrderPara-1,'omitnan')./max(extctOrderAll);
        extctStatisticFree = mean(extctOrderFree-1,'omitnan')./max(extctOrderAll);
        
        extctStatisticParaAll{count}(ii,:) = extctStatisticPara;
        extctStatisticFreeAll{count}(ii,:) = extctStatisticFree;
        
        extctOrderAll(isnan(extctOrderAll)) = inf;
        
        conExtctFirst = (extctOrderAll(LL(:,1),:) > extctOrderAll(LL(:,2),:))*1;
        resExtctFirst = (extctOrderAll(LL(:,1),:) < extctOrderAll(LL(:,2),:))*1;
        
        conExtctFirst(basal(LL(:,1)),:) = nan;
        resExtctFirst(basal(LL(:,1)),:) = nan;
        
        
        conExtctFirst(isinf(extctOrderAll(LL(:,1),:))&isinf(extctOrderAll(LL(:,2),:))) = nan;
        resExtctFirst(isinf(extctOrderAll(LL(:,1),:))&isinf(extctOrderAll(LL(:,2),:))) = nan;
        
        paraExtctFirst = conExtctFirst;
        hostExtctFirst = resExtctFirst;
        
        paraExtctFirst(~paraMxBin(LL(:,2),:)) = nan;
        hostExtctFirst(~paraMxBin(LL(:,2),:)) = nan;
        
        
        fracParaExtctFirstAll{count}(ii,:) = mean(paraExtctFirst,'omitnan');
        fracHostExtctFirstAll{count}(ii,:) = mean(hostExtctFirst,'omitnan');
        
        fracConExtctFirstAll_genAveraged{count}(ii,:) = sum(conExtctFirst./gen(LL(:,2)),'omitnan')/S;
        fracResExtctFirstAll_vulAveraged{count}(ii,:) = sum(conExtctFirst./vul(LL(:,1)),'omitnan')/S;
        
        fracConExtctFirstAll{count}(ii,:) = mean(conExtctFirst,'omitnan');
        fracResExtctFirstAll{count}(ii,:) = mean(resExtctFirst,'omitnan');
        
        patlMx = repmat(patl,1,nFPar);
        bodySizeMx = zeros(S,nFPar);
        
        bodySizeMx(~paraMxBin) = (10^kFree(model(1))).^...
            (patlMx(~paraMxBin) - 1);
        
        
        bodySizeMx(paraMxBin) = 10.^(kPara(model(2)) + ...
            kFree(model(1))*(patlMx(paraMxBin) - 2));
        
        xAllSpecies = .314*bodySizeMx.^-.25;
        
        bodySizeRatios = bodySizeMx(LL(:,2),:)./bodySizeMx(LL(:,1),:);
        
        extctMx = onlyMeans<1e-10;
        
        xAllWebs{count,ii} = xAllSpecies;
        orderAllWebs{count,ii} = extctOrderAll;
        
        bsrAll{count} = [bsrAll{count};bodySizeRatios];
        LLAll{count} = [LLAll{count};LL];
        extctAll{count} = [extctAll{count};extctMx];
        
        resExtctAll{count} = [resExtctAll{count};...
            extctAll{count}(LL(:,1),:)];
        resBodySizeAll{count} = [resBodySizeAll{count};...
            bodySizeMx(LL(:,1),:)];
        
        conExtctAll{count} = [conExtctAll{count};...
            extctAll{count}(LL(:,2),:)];
        conBodySizeAll{count} = [conBodySizeAll{count};...
            bodySizeMx(LL(:,2),:)];
        
        freeMeans = onlyMeans;
        freeMeans(paraMxBin) = nan;
        
        
        freeAlive = (freeMeans>1e-10)*1;
        freeAlive(paraMxBin) = nan;
        
        paraMeans = onlyMeans;
        paraMeans(~paraMxBin) = nan;
        
        paraAlive = (paraMeans>1e-10)*1;
        paraAlive(~paraMxBin) = nan;
        
        persistenceAll{count}(ii,:) = mean(finals>1e-10);
        persistenceFree{count}(ii,:) = mean(freeAlive(~basal,:),'omitnan');
        persistencePara{count}(ii,:) = mean(paraAlive(~basal,:),'omitnan');
        
        means(means==0) = nan;
        freeMeans(freeMeans==0) = nan;
        paraMeans(paraMeans==0) = nan;
        
        avgBiomassAll{count}(ii,:) = sum(log10(means),'omitnan')./persistenceAll{count}(ii,:)./S;
        avgBiomassFree{count}(ii,:) = sum(log10(freeMeans(~basal,:)),'omitnan')./persistenceFree{count}(ii,:)./S;
        avgBiomassPara{count}(ii,:) = sum(log10(paraMeans(~basal,:)),'omitnan')./persistencePara{count}(ii,:)./S;
        
        nAliveAll{count}(ii,:)  = sum(means>1e-10);
        nAliveFree{count}(ii,:) = sum(freeAlive(~basal,:)>1e-10);
        nAlivePara{count}(ii,:) = sum(paraAlive(~basal,:)>1e-10);
        
        vulBasal = props(basal,3)*linkDensity;
        meanVulBasal = mean(vulBasal);
        sumVulBasalAll(ii,count)  = sum(vulBasal);
        meanVulBasalAll(ii,count) = meanVulBasal;
        nBasal(ii,count) = sum(basal);
        if count == 1
            paraMxBinAll = [paraMxBinAll;paraMxBin];
            patlAll = [patlAll;patl];
        end
        
        %couple things to add:
        %1. extinction orders: how many of your prey go extinct before you?
        % how many of your predators?  Do parasites or free-livers tend to
        % go extinct first?  That kind of thing.
        %2. how your properties change as species go extinct.  related to
        %previous point, but requires re-calculation of local (and global!)
        %properties.
        
        if (numel(idxPar)+sum(basal))~=40
            fprintf('break\n')
        end
        
    end
    meanPersistenceAll = mean(persistenceAll{count});
    meanPersistenceFree = mean(persistenceFree{count});
    meanPersistencePara = mean(persistencePara{count});
    
    stdPersistenceAll = std(persistenceAll{count});
    stdPersistenceFree = std(persistenceFree{count});
    stdPersistencePara = std(persistencePara{count});
    
    if numel(meanPersistenceAll) == nFPar0
        meansAll(count,:) = meanPersistenceAll;
        stdsAll(count,:) = stdPersistenceAll;
        
        meansAll(models(:,1) == model(1),1) = meanPersistenceAll(1);
        stdsAll(models(:,1) == model(1),1) = stdPersistenceAll(1);
        
        meansAllFree(models(:,1) == model(1),1) = meanPersistenceFree(1);
        stdsAllFree(models(:,1) == model(1),1) = stdPersistenceFree(1);
        
        meansAllPara(models(:,1) == model(1),1) = meanPersistencePara(1);
        stdsAllPara(models(:,1) == model(1),1) = stdPersistencePara(1);
        
        meansAllPara(count,:) = meanPersistencePara;
        stdsAllPara(count,:) = stdPersistencePara;
        
        meansAllFree(count,:) = meanPersistenceFree;
        stdsAllFree(count,:) = stdPersistenceFree;
        
        
    else
        
        meansAll(count,2:end) = meanPersistenceAll;
        stdsAll(count,2:end) = stdPersistenceAll;
        
        meansAllPara(count,2:end) = meanPersistencePara;
        stdsAllPara(count,2:end) = stdPersistencePara;
        
        meansAllFree(count,2:end) = meanPersistenceFree;
        stdsAllFree(count,2:end) = stdPersistenceFree;
    end
    
    %Too look at a single fraction of parasites, get rid of these reshape
    %commands, then modify the use of conExtctAll below.
    %bsrAll{count} = reshape(bsrAll{count},[],1);
    %conExtctAll{count} = reshape(conExtctAll{count},[],1);
    %resExtctAll{count} = reshape(resExtctAll{count},[],1);
    
    %conBodySizeAll{count} = reshape(conBodySizeAll{count},[],1);
    %resBodySizeAll{count} = reshape(resBodySizeAll{count},[],1);
    
    groups{count} = conExtctAll{count} + 2*resExtctAll{count};
    
end
%}
ZFrees = false(nModels,2);
ZParas = false(nModels,2);
freeLivings = false(nModels,2);
concs = false(nModels,2);

ZFrees(:,1) = models(:,1) == 1; %Small Free
ZFrees(:,2) = models(:,1) == 2; %big Free
ZParas(:,1) = models(:,2) == 1; %Big Para
ZParas(:,2) = models(:,2) == 2; % SMall Para
freeLivings(:,1) = models(:,3) == 1; %Don't include free LIving
freeLivings(:,2) = models(:,3) == 2; %INclude Free Living
concs(:,1) = models(:,4) == 1; %Don't include concomittant
concs(:,2) = models(:,4) == 2; %INclude concomittant
AXW = fullfact([2,2]);

%{
%{
for ii = 1:9
    figure('Position',[0,0,1440,900]);
    %subplot by model type.
    for jj = 1:4
        subplot(2,2,jj)
        %Group by Z pairs
        hold on
        for kk = 1:4
            pickYs = ZFrees(:,AXW(kk,1))&...
                ZParas(:,AXW(kk,2))&...
                freeLivings(:,AXW(jj,1))&...
                concs(:,AXW(jj,2));
            y = persistenceAllFree{pickYs}(:,ii);
            x = nBasal(:,pickYs);
            plot(x,y,'o')
        end
        axis([1 10 0 1]);
    end
end
%}
%{
marks = {'ko','ro';'ks','rs'};
legendEntries = {'No FF and No Conc',...
    'No FF and Conc',...
    'FF and No Conc',...
    'FF and Conc'...
    };

for Zpairs = 1:4;
    ZPairBin = ZFrees(:,AXW(Zpairs,1))&ZParas(:,AXW(Zpairs,2));
    
    

    bsScatterHistFig = figure('Position', [0,0,1440,900]);
    
    %Modify conBodySizeAll{1}(:,1) to conBodySizeAll{1}(:,fPar) after
    %making changes to conBodySizeAll above to plot different fractions of
    %parasties.
    [~,pickFewer] = datasample(conBodySizeAll{1}(:,4),1000);
    pickFewer = ones(size(conBodySizeAll{1}(:,1)))==1;
    for subplotChoice = 1:4;

        subplot(3,4,subplotChoice);
        
        pickModel = freeLivings(:,AXW(subplotChoice,1))&...
            concs(:,AXW(subplotChoice,2))&...
            ZPairBin;
        %Modify conBodySizeAll{pickModel} to
        %conBodySizeAll{pickModel}(:,fPar) after making changes to
        %conBodySizeAll above to plot different fractions of parasties.
        %similarly for resBodySizeAll
        
        y = conBodySizeAll{pickModel}(:,4);
        x = resBodySizeAll{pickModel}(:,4);
        
        if mod(subplotChoice,2)==1
            g = groups{pickModel}(:,4);
        else
            g = groups{pickModel}(:,4);
        end
        g = g==0;
        gNames = cell(size(g));
        gNames(g) = {'alive'};
        gNames(~g) = {'dead'};
        h = gscatter(log10(x(pickFewer)),log10(y(pickFewer)),g(pickFewer));
        
        
        refline(1,0)
        h(1).Color = 'r';
        h(1).DisplayName = 'Dead';
        h(2).Color = 'b';
        h(2).DisplayName = 'Alive';
        bsr = log10(y./x);
        minBsr = min(bsr);
        maxBsr = max(bsr);
        edges = linspace(minBsr,maxBsr,50);
        
        
        subplot(3,4,subplotChoice+4)
        histogram(bsr(g==1),edges,'FaceColor','b');
        ylabel(sprintf('Count (Total = %u)',sum(g==1)))
        xlabel(sprintf('log_{10}(BSR,\\mu=%.2f',mean(bsr(g==1))));
        title('Surviving Links')
        
        subplot(3,4,subplotChoice+8)
        histogram(bsr(g==0),edges,'FaceColor','r');
        ylabel(sprintf('Count (Total = %u)',sum(g==0)))
        xlabel(sprintf('log_{10}(BSR,\\mu=%.2f',mean(bsr(g==0))));
        title('Dead Links')
        
        
    end
    hold off
    subplot(3,4,1);
    title(sprintf('No Free Stage,No Conc.'));
    ylabel('Consumer Body Size')
    xlabel('Resource Body Size')
    subplot(3,4,2);
    title(sprintf('Free Stage, No Conc.'));
    ylabel('Consumer Body Size')
    xlabel('Resource Body Size')
%     axesPosition = get(gca,'Position');
%     hNewAxes = axes('Position',axesPosition,...
%         'Color','none',...
%         'YAxisLocation','right',...
%         'XTick',[],...
%         'YTick',[],...
%         'Box','off');
%     ylabel(hNewAxes,'No Concomittant','FontWeight','bold')
    subplot(3,4,3);
    title(sprintf('No Free Stage, Conc.'));
    ylabel('Consumer Body Size')
    xlabel('Resource Body Size')
    subplot(3,4,4);
    title(sprintf('Free Stage, Conc.'));
    xlabel('Resource Body Size')
    ylabel('Consumer Body Size')
%     axesPosition = get(gca,'Position');
%     hNewAxes = axes('Position',axesPosition,...
%         'Color','none',...
%         'YAxisLocation','right',...
%         'XTick',[],...
%         'YTick',[],...
%         'Box','off');
%     ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    saveas(bsScatterHistFig,sprintf('bsScatterHist%u',Zpairs),'jpg');
end
%hsaveas(fig1,sprintf('ExperimentFigure1b'),'jpg')
%}

%{
noCnoFF = concs(:,1)&freeLivings(:,1);
CnoFF = concs(:,2)&freeLivings(:,1);

p1x = noCnoFF&ZFrees(:,1)&ZParas(:,1);
p2x = noCnoFF&ZFrees(:,2)&ZParas(:,1);
p3x = noCnoFF&ZFrees(:,1)&ZParas(:,2);
p4x = noCnoFF&ZFrees(:,2)&ZParas(:,2);

p1y = CnoFF&ZFrees(:,1)&ZParas(:,1);
p2y = CnoFF&ZFrees(:,2)&ZParas(:,1);
p3y = CnoFF&ZFrees(:,1)&ZParas(:,2);
p4y = CnoFF&ZFrees(:,2)&ZParas(:,2);
%}
%{
for ii = 1:9
    biomassfig = figure('Position', [0,0,1440,900]);
    subplot(12,3,[1 4 7 10]);
    
    x1 = persistenceAll{p1x}(:,ii);
    y1 = persistenceAll{p1y}(:,ii);
    x2 = persistenceAll{p2x}(:,ii);
    y2 = persistenceAll{p2y}(:,ii);
    x3 = persistenceAll{p3x}(:,ii);
    y3 = persistenceAll{p3y}(:,ii);
    x4 = persistenceAll{p4x}(:,ii);
    y4 = persistenceAll{p4y}(:,ii);
    
    x = [x1;x2;x3;x4];
    y = [y1;y2;y3;y4];
    
    [nx1,~] = histcounts(x1,linspace(0,1,41));
    [nx2,~] = histcounts(x2,linspace(0,1,41));
    [nx3,~] = histcounts(x3,linspace(0,1,41));
    [nx4,edgesx] = histcounts(x4,linspace(0,1,41));
    
    [ny1,~] = histcounts(y1,linspace(0,1,41));
    [ny2,~] = histcounts(y2,linspace(0,1,41));
    [ny3,~] = histcounts(y3,linspace(0,1,41));
    [ny4,edgesy] = histcounts(y4,linspace(0,1,41));
    
    groups = [ones(size(x1));
              2*ones(size(x2));
              3*ones(size(x3));
              4*ones(size(x4));];
    gscatter(x,y,groups,[],'sox+')
    line(mean(x)*ones(size(0:.1:1)),0:.1:1,'Color','red','LineStyle','--','LineWidth',0.01)
    line(0:.1:1,mean(y)*ones(size(0:.1:1)),'Color','red','LineStyle','--','LineWidth',0.01)
axis([0 1 0 1])
    refline(1,0)
    
    
    titleStr = sprintf(' fpar = %.2f\nPeristence of all Species at Final State',fParAll(ii));
    title(titleStr)
    xlabel('Persistence without Concomittant Links')
    ylabel('Persistence with concomittant links')
    subplot(12,3,13)
    histogram(x1,edgesx)
    subplot(12,3,16)
    histogram(x2,edgesx)
    subplot(12,3,19)
    histogram(x3,edgesx)
    subplot(12,3,22)
    histogram(x4,edgesx)
    xlabel('(Persistence Without Concomittant Predation')
    subplot(12,3,25)
    histogram(y1,edgesy)
    subplot(12,3,28)
    histogram(y2,edgesy)
    subplot(12,3,31)
    histogram(y3,edgesy)
    subplot(12,3,34)
    histogram(y4,edgesy)
    xlabel('Persistence With Concomittant Predation')
    subplot(12,3,[2 5 8 11]);
    
    x1 = persistenceAllFree{p1x}(:,ii);
    y1 = persistenceAllFree{p1y}(:,ii);
    x2 = persistenceAllFree{p2x}(:,ii);
    y2 = persistenceAllFree{p2y}(:,ii);
    x3 = persistenceAllFree{p3x}(:,ii);
    y3 = persistenceAllFree{p3y}(:,ii);
    x4 = persistenceAllFree{p4x}(:,ii);
    y4 = persistenceAllFree{p4y}(:,ii);
    
    x = [x1;x2;x3;x4];
    y = [y1;y2;y3;y4];
    
    [~,edgesx] = histcounts(x,'BinLimits',[0,1]);
    [~,edgesy] = histcounts(y,'BinLimits',[0,1]);
    
    [nx1] = histcounts(x1,edgesx);
    [nx2] = histcounts(x2,edgesx);
    [nx3] = histcounts(x3,edgesx);
    [nx4] = histcounts(x4,edgesx);
    
    [ny1] = histcounts(y1,edgesy);
    [ny2] = histcounts(y2,edgesy);
    [ny3] = histcounts(y3,edgesy);
    [ny4] = histcounts(y4,edgesy);
    
    groups = [ones(size(x1));
              2*ones(size(x2));
              3*ones(size(x3));
              4*ones(size(x4));];
    gscatter(x,y,groups,[],'sox+')
    line(mean(x)*ones(size(0:.1:1)),0:.1:1,'Color','red','LineStyle','--','LineWidth',0.01)
    line(0:.1:1,mean(y)*ones(size(0:.1:1)),'Color','red','LineStyle','--','LineWidth',0.01)
        axis([0 1 0 1])
        refline(1,0)

    title('Persistence of Free Living Consumers')
    xlabel('Persistence without Concomittant Links')
    ylabel('Persistence with concomittant links')
    subplot(12,3,14)
    histogram(x1,edgesx)
    subplot(12,3,17)
    histogram(x2,edgesx)
    subplot(12,3,20)
    histogram(x3,edgesx)
    subplot(12,3,23)
    histogram(x4,edgesx)
    xlabel('(Persistence Without Concomittant Predation')
    subplot(12,3,26)
    histogram(y1,edgesy)
    subplot(12,3,29)
    histogram(y2,edgesy)
    subplot(12,3,32)
    histogram(y3,edgesy)
    subplot(12,3,35)
    histogram(y4,edgesy)
    xlabel('Persistence With Concomittant Predation')
    
    
    subplot(12,3,[3 6 9]);
    
    x1 = persistenceAllPara{p1x}(:,ii);
    y1 = persistenceAllPara{p1y}(:,ii);
    x2 = persistenceAllPara{p2x}(:,ii);
    y2 = persistenceAllPara{p2y}(:,ii);
    x3 = persistenceAllPara{p3x}(:,ii);
    y3 = persistenceAllPara{p3y}(:,ii);
    x4 = persistenceAllPara{p4x}(:,ii);
    y4 = persistenceAllPara{p4y}(:,ii);
    
    x = [x1;x2;x3;x4];
    y = [y1;y2;y3;y4];
    
    [nx,edgesx] = histcounts(x,'BinLimits',[0,1]);
    [ny,edgesy] = histcounts(y,'BinLimits',[0,1]);
    
    [nx1,edgesx1] = histcounts(x1,edgesx);
    [nx2,edgesx2] = histcounts(x2,edgesx);
    [nx3,edgesx3] = histcounts(x3,edgesx);
    [nx4,edgesx4] = histcounts(x4,edgesx);
    
    [ny1,edgesy1] = histcounts(y1,edgesy);
    [ny2,edgesy2] = histcounts(y2,edgesy);
    [ny3,edgesy3] = histcounts(y3,edgesy);
    [ny4,edgesy4] = histcounts(y4,edgesy);
    
    groups = [ones(size(x1));
              2*ones(size(x2));
              3*ones(size(x3));
              4*ones(size(x4));];
    gscatter(x,y,groups,[],'sox+')
    line(mean(x)*ones(size(0:.1:1)),0:.1:1,'Color','red','LineStyle','--','LineWidth',0.01)
    line(0:.1:1,mean(y)*ones(size(0:.1:1)),'Color','red','LineStyle','--','LineWidth',0.01)
    axis([0 1 0 1])
    refline(1,0)
    
    title('Persistence of Parasites')
    xlabel('Persistence without Concomittant Links')
    ylabel('Persistence with concomittant links')
    subplot(12,3,15)
    histogram(x1,edgesx)
    subplot(12,3,18)
    histogram(x2,edgesx)
    subplot(12,3,21)
    histogram(x3,edgesx)
    subplot(12,3,24)
    histogram(x4,edgesx)
    xlabel('(Persistence Without Concomittant Predation')
    subplot(12,3,27)
    histogram(y1,edgesy)
    subplot(12,3,30)
    histogram(y2,edgesy)
    subplot(12,3,33)
    histogram(y3,edgesy)
    subplot(12,3,36)
    histogram(y4,edgesy)
    xlabel('Persistence With Concomittant Predation')
    
    saveas(biomassfig,sprintf('persistenceFig%u',ii),'jpg')
    close(biomassfig);
end
%}

%{
noCnoFF = concs(:,1)&freeLivings(:,1);
CnoFF = concs(:,2)&freeLivings(:,1);

p1x = noCnoFF&ZFrees(:,1)&ZParas(:,1);
p2x = noCnoFF&ZFrees(:,2)&ZParas(:,1);
p3x = noCnoFF&ZFrees(:,1)&ZParas(:,2);
p4x = noCnoFF&ZFrees(:,2)&ZParas(:,2);

p1y = CnoFF&ZFrees(:,1)&ZParas(:,1);
p2y = CnoFF&ZFrees(:,2)&ZParas(:,1);
p3y = CnoFF&ZFrees(:,1)&ZParas(:,2);
p4y = CnoFF&ZFrees(:,2)&ZParas(:,2);
%}
%{
for ii = 1:9
    biomassfig = figure('Position', [0,0,1200,300]);
    subplot(1,3,1);
    
    x1 = totalBiomassAll{p1x}(:,ii)./(40*persistenceAll{p1y}(:,ii));
    y1 = totalBiomassAll{p1y}(:,ii)./(40*persistenceAll{p1y}(:,ii));
    x2 = totalBiomassAll{p2x}(:,ii)./(40*persistenceAll{p2x}(:,ii));
    y2 = totalBiomassAll{p2y}(:,ii)./(40*persistenceAll{p2y}(:,ii));
    x3 = totalBiomassAll{p3x}(:,ii)./(40*persistenceAll{p3x}(:,ii));
    y3 = totalBiomassAll{p3y}(:,ii)./(40*persistenceAll{p3y}(:,ii));
    x4 = totalBiomassAll{p4x}(:,ii)./(40*persistenceAll{p4x}(:,ii));
    y4 = totalBiomassAll{p4y}(:,ii)./(40*persistenceAll{p4y}(:,ii));
    
    x = log10([x1;x2;x3;x4]);
    y = log10([y1;y2;y3;y4]);
    
    groups = [ones(size(x1));
              2*ones(size(x2));
              3*ones(size(x3));
              4*ones(size(x4));];
    gscatter(x,y,groups)
    refline(1,0)
    titleStr = sprintf('All Species at Final State, fpar = %.2f',fParAll(ii));
    title(titleStr)
    xlabel('Biomass per species witout Concomittant Links')
    ylabel('Biomass per species wit concomittant links')
    subplot(1,3,2);
    
    x1 = totalBiomassFree{p1x}(:,ii)./(40*persistenceAllFree{p1y}(:,ii));
    y1 = totalBiomassFree{p1y}(:,ii)./(40*persistenceAllFree{p1y}(:,ii));
    x2 = totalBiomassFree{p2x}(:,ii)./(40*persistenceAllFree{p2x}(:,ii));
    y2 = totalBiomassFree{p2y}(:,ii)./(40*persistenceAllFree{p2y}(:,ii));
    x3 = totalBiomassFree{p3x}(:,ii)./(40*persistenceAllFree{p3x}(:,ii));
    y3 = totalBiomassFree{p3y}(:,ii)./(40*persistenceAllFree{p3y}(:,ii));
    x4 = totalBiomassFree{p4x}(:,ii)./(40*persistenceAllFree{p4x}(:,ii));
    y4 = totalBiomassFree{p4y}(:,ii)./(40*persistenceAllFree{p4y}(:,ii));
    
    x = log10([x1;x2;x3;x4]);
    y = log10([y1;y2;y3;y4]);
    
    groups = [ones(size(x1));
              2*ones(size(x2));
              3*ones(size(x3));
              4*ones(size(x4));];
    gscatter(x,y,groups)
    refline(1,0)
    title('Free Living Consumers at Final State')
    xlabel('Biomass per species witout Concomittant Links')
    ylabel('Biomass per species wit concomittant links')
    subplot(1,3,3);
    
    x1 = totalBiomassPara{p1x}(:,ii)./(40*persistenceAllPara{p1y}(:,ii));
    y1 = totalBiomassPara{p1y}(:,ii)./(40*persistenceAllPara{p1y}(:,ii));
    x2 = totalBiomassPara{p2x}(:,ii)./(40*persistenceAllPara{p2x}(:,ii));
    y2 = totalBiomassPara{p2y}(:,ii)./(40*persistenceAllPara{p2y}(:,ii));
    x3 = totalBiomassPara{p3x}(:,ii)./(40*persistenceAllPara{p3x}(:,ii));
    y3 = totalBiomassPara{p3y}(:,ii)./(40*persistenceAllPara{p3y}(:,ii));
    x4 = totalBiomassPara{p4x}(:,ii)./(40*persistenceAllPara{p4x}(:,ii));
    y4 = totalBiomassPara{p4y}(:,ii)./(40*persistenceAllPara{p4y}(:,ii));
    
    x = log10([x1;x2;x3;x4]);
    y = log10([y1;y2;y3;y4]);
    
    groups = [ones(size(x1));
              2*ones(size(x2));
              3*ones(size(x3));
              4*ones(size(x4));];
    gscatter(x,y,groups)
    refline(1,0)
    title('Parasites at Final State')
    xlabel('Biomass per species witout Concomittant Links')
    ylabel('Biomass per species wit concomittant links')
    saveas(biomassfig,sprintf('biomassFig%u',ii),'jpg')
    close(biomassfig);
end
%}
%{
AXW = fullfact([2,2]);

marks = {'ko','ro';'ks','rs'};
legendEntries = {'No FF and No Conc',...
    'No FF and Conc',...
    'FF and No Conc',...
    'FF and Conc'...
    };

fig1 = figure('Position', [0,0,1440,900]);

for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for freeLivingsChoice = 1:2;
        for concsChoice = 1:2
            plotStyles = marks{freeLivingsChoice,concsChoice};
            pickYs = concs(:,concsChoice)&...
                ZParas(:,AXW(subplotChoice,1))&...
                ZFrees(:,AXW(subplotChoice,2))&...
                freeLivings(:,freeLivingsChoice);
            y1 = meansAll(pickYs,:);
            plot(fParAll,y1,plotStyles)
            axis([0 0.8 0 1.05]);
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('Big Parasites'));
    ylabel('Fraction All Species Persistence')
    legend(legendEntries,'Location','Best')
    subplot(2,2,2);
    title(sprintf('Small Parasites'));
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Small Free Livers','FontWeight','bold')
    subplot(2,2,3);
    ylabel('Fraction All Species Persistence')
    xlabel('Fraction of Consumers as Parasites')
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Big Free Livers','FontWeight','bold')
    
end
saveas(fig1,sprintf('ExperimentFigure1b'),'jpg')

fig2 = figure('Position', [0,0,1440,900]);
title('Free Liver Persistence')
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for freeLivingsChoice = 1:2;
        for concsChoice = 1:2
            plotStyles = marks{freeLivingsChoice,concsChoice};
            pickYs = concs(:,concsChoice)&...
                ZParas(:,AXW(subplotChoice,1))&...
                ZFrees(:,AXW(subplotChoice,2))&...
                freeLivings(:,freeLivingsChoice);
            y1 = meansAllFree(pickYs,:);
            plot(fParAll,y1,plotStyles)
            axis([0 0.8 0 1.05]);
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('Big Parasites'));
    ylabel('Fraction Free Consumer Persistence')
    legend(legendEntries,'Location','Best')
    subplot(2,2,2);
    title(sprintf('Small Parasites'));
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Small Free Livers','FontWeight','bold')
    subplot(2,2,3);
    ylabel('Fraction Free-Consumer Persistence')
    xlabel('Fraction of Consumers as Parasites')
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Big Free Livers','FontWeight','bold')
    
end
saveas(fig2,sprintf('ExperimentFigure2b'),'jpg')

fig3 = figure('Position', [0,0,1440,900]);
title('Parasite Persistence')
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for freeLivingsChoice = 1:2;
        for concsChoice = 1:2
            plotStyles = marks{freeLivingsChoice,concsChoice};
            pickYs = concs(:,concsChoice)&...
                ZParas(:,AXW(subplotChoice,1))&...
                ZFrees(:,AXW(subplotChoice,2))&...
                freeLivings(:,freeLivingsChoice);
            y1 = meansAllPara(pickYs,:);
            plot(fParAll,y1,plotStyles)
            axis([0 0.8 0 1.05]);
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('Big Parasites'));
    ylabel('Fraction parasite Persistence')
    legend(legendEntries,'Location','Best')
    subplot(2,2,2);
    title(sprintf('Small Parasites'));
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Small Free Livers','FontWeight','bold')
    subplot(2,2,3);
    ylabel('Fraction Parasite Persistence')
    xlabel('Fraction of Consumers as Parasites')
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Big Free Livers','FontWeight','bold')
    
end
saveas(fig3,sprintf('ExperimentFigure3b'),'jpg')
%}
%}

%creates figures organized in a way I didn't like any more.  Better to do
%it a different way?
% fParAll = linspace(0,1,5);

CI = true;
sigLevel = .05;
nullChoice = ZFrees(:,2)+1;

xmin = min(fParAll)-.05;
xmax = max(fParAll)+.05;
ymin = -0.05;
ymax = 1.05;

middlePercentage = .5;
lowerQuantile = (1-middlePercentage)/2;
upperQuantile = 1 - lowerQuantile;

AXW = fullfact([2,2]);
ZParaColors = {'r','b'};
ZFreeMarks = {'o','x'};

marksMatrix = {{'k^','MarkerFaceColor','k'},{'x','Color',[.75,0,.75]};
    {'v','MarkerFaceColor',[0,.5,.5]}, {'+','Color',[0 0.75 0]}};

legendEntries = {'(Z_{Free},Z_{Para})=(10,1e-3)'...
    ,'(Z_{Free},Z_{Para})=(10,1e-4)'...
    ,'(Z_{Free},Z_{Para})=(100,1e-3)'...
    ,'(Z_{Free},Z_{Para})=(100,1e-4)'...
    };

param = struct('xVals',fParAll...
    ,'makePlot',true...
    ,'savePlot',false...
    ,'saveData',true...
    ,'statFcn',@mean...
    ,'statFcnOpt','omitnan'...
    ,'errBars',struct...
    ,'plotParam',struct...
    ,'dataParam',struct...
    ,'modelCodes',models...
    ,'dataAnalyzed','persistenceAll');

param.errBars = struct('errBarsFcn',@calculateConfidenceIntervals...
    ,'param',struct);

param.errBars.param = struct('alpha',.05);

param.plotParam = struct('marksMatrix',{marksMatrix} ...
    ,'legendEntries',{legendEntries} ...
    ,'title1',sprintf('No Free-Living Stage')...
    ,'title2',sprintf('Free-Living Stage')  ...
    ,'y2label1',sprintf('No Concomittant') ...
    ,'y2label2',sprintf('Concomittant') ...
    ,'ylabel',sprintf('Fraction of All Species %s'...
    ,param.dataAnalyzed) ...
    ,'xlabel',sprintf('Fraction of Consumers as Parasites') ...
    ,'plotRange',[xmin,xmax,ymin,ymax] ...
    ,'filename','AllSpeciesPersistenceFig'...
    ,'format','jpg');

param.dataParam =...
    struct('headers','x,y,yp,ym' ...
    ,'filenames',{{'AllSpecies00';...
    'AllSpecies0F';...
    'AllSpeciesC0';...
    'AllSpeciesCF';}}...
    ,'filePath','../../Data/plotData/'...
    );

param.regressionParam = ...
    struct('doRegression',true...
           ,'cmd',@fitlm...
           ,'saveReg',true...
           ,'regOpts',struct);

processExperimentOutput(persistenceAll,param);

param.dataAnalyzed = 'free-persistence';
param.plotParam.ylabel = 'Persistence of Free Livers';
processExperimentOutput(persistenceFree,param);

param.dataAnalyzed = 'para-persistence';
param.plotParam.ylabel = 'Persistence of Parasites';
processExperimentOutput(persistencePara,param);

param.dataAnalyzed = 'biomass';
param.plotParam.ylabel = 'Average Ecosystem Biomass';
processExperimentOutput(avgBiomassAll,param);

param.dataAnalyzed = 'parasitic-biomass';
param.plotParam.ylabel = 'Averge Parasitic Biomass';
processExperimentOutput(avgBiomassPara,param);

param.dataAnalyzed = 'free-living-biomass';
param.plotParam.ylabel = 'Averge Free-Living Biomass';
processExperimentOutput(avgBiomassFree,param);

param.dataAnalyzed = 'para-extct-order-stat';
param.plotParam.ylabel = 'Average Fraction of Species Extinct Before Parasites';
processExperimentOutput(extctStatisticParaAll,param);

param.dataAnalyzed = 'free-living-biomass';
param.plotParam.ylabel = 'Averge Fraction of Species Extinct Before Free';
processExperimentOutput(extctStatisticFreeAll,param);

param.dataAnalyzed = 'frac-con-first';
param.plotParam.ylabel = 'Fraction of links with consumer extinct first';
processExperimentOutput(fracConExtctFirstAll,param);

param.dataAnalyzed = 'frac-res-first';
param.plotParam.ylabel = 'Fraction of links with resource extinct first';
processExperimentOutput(fracResExtctFirstAll,param);

param.dataAnalyzed = 'frac-res-first-avgd';
param.plotParam.ylabel = 'avgd Fraction of links with resource extinct first';
processExperimentOutput(fracResExtctFirstAll_vulAveraged,param);

param.dataAnalyzed = 'frac-con-first-avgd';
param.plotParam.ylabel = 'avgd Fraction of links with resource extinct first';
processExperimentOutput(fracConExtctFirstAll_genAveraged,param);

param.dataAnalyzed = 'frac-para-first';
param.plotParam.ylabel = 'Fraction of para links with para extinct first';
processExperimentOutput(fracParaExtctFirstAll,param);

param.dataAnalyzed = 'frac-Host-first';
param.plotParam.ylabel = 'Fraction of para links with host extinct first';
processExperimentOutput(fracHostExtctFirstAll,param);

%{
fig1 = figure('Position', [0,0,1440,900]);

med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);
y10 = zeros(2,1);
ME10 = zeros(2,1);
ymax = 0;
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            y0 = totalBiomassAll{pickYs}./nAliveAll{pickYs};
            y1 = mean(y0,'omitnan');
            s1 = std(y0,'omitnan');
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(y0,0.5);
            ql = abs(quantile(y0,lowerQuantile)-med);
            qu = abs(quantile(y0,upperQuantile)-med);
            
            
            if numel(med) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                y10(nullChoice(pickYs)) = y1(1);
                ME10(nullChoice(pickYs)) = s1(1);
                
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
                
                y1 = [y10(nullChoice(pickYs)) y1];
                ME = [ME10(nullChoice(pickYs)) ME];
            end
            
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            ymax = max(ymax,max(y1));
            grid on
        end
    end
    ymax = 1.2*ymax;
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
   ylabel('Per-Capita Biomass')
    legend(legendEntries,'Location','Best')
    axis([xmin xmax ymin ymax]);
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
    axis([xmin xmax ymin ymax]);
    subplot(2,2,3);
   ylabel('Per-Capita Biomass')
    xlabel('Fraction of Consumers as Parasites')
    axis([xmin xmax ymin ymax]);
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    axis([xmin xmax ymin ymax]);
    
end
saveas(fig1,sprintf('BiomassCIPerAll'),'jpg')

fig2 = figure('Position', [0,0,1440,900]);

med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);
y10 = zeros(2,1);
ME10 = zeros(2,1);
ymax = 0;

for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            y0 =  totalBiomassFree{pickYs}./nAliveFree{pickYs};
            y1 = mean(y0,'omitnan');
            s1 = std(y0,'omitnan');
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(y0,0.5);
            ql = abs(quantile(y0,lowerQuantile)-med);
            qu = abs(quantile(y0,upperQuantile)-med);
            
            
            if numel(med) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                y10(nullChoice(pickYs)) = y1(1);
                ME10(nullChoice(pickYs)) = s1(1);
                
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
                
                y1 = [y10(nullChoice(pickYs)) y1];
                ME = [ME10(nullChoice(pickYs)) ME];
            end
            
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            ymax = max(ymax,max(y1));
            grid on
        end
    end
    ymax = 1.2*ymax;
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
   ylabel('Per-Capita Free-Living Biomass')
    legend(legendEntries,'Location','Best')
    axis([xmin xmax ymin ymax]);
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
    axis([xmin xmax ymin ymax]);
    subplot(2,2,3);
   ylabel('Per-Capita Free-Living Biomass')
    xlabel('Fraction of Consumers as Parasites')
    axis([xmin xmax ymin ymax]);
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    axis([xmin xmax ymin ymax]);
    
end
saveas(fig2,sprintf('BiomassCIPerFree'),'jpg')

fig3 = figure('Position', [0,0,1440,900]);

med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);
y10 = zeros(2,1);
ME10 = zeros(2,1);
ymax = 0;
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            y0 = totalBiomassPara{pickYs}./nAlivePara{pickYs};
            y0(isnan(y0)) = 0;
            
            
            y1 = mean(y0,'omitnan');
            s1 = std(y0,'omitnan');
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(y0,0.5);
            ql = abs(quantile(y0,lowerQuantile)-med);
            qu = abs(quantile(y0,upperQuantile)-med);
            
            
            if numel(med) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                y10(nullChoice(pickYs)) = y1(1);
                ME10(nullChoice(pickYs)) = s1(1);
                
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
                
                y1 = [y10(nullChoice(pickYs)) y1];
                ME = [ME10(nullChoice(pickYs)) ME];
            end
            
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            ymax = max(ymax,max(y1));
            grid on
        end
    end
end
    ymax = 1.2*ymax;
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
   ylabel('Per-Capita Parasite Biomass')
    legend(legendEntries,'Location','Best')
    axis([xmin xmax ymin ymax]);
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
    axis([xmin xmax ymin ymax]);
    subplot(2,2,3);
   ylabel('Per-Capita Parasite Biomass')
    xlabel('Fraction of Consumers as Parasites')
    axis([xmin xmax ymin ymax]);
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    axis([xmin xmax ymin ymax]);
    
saveas(fig3,sprintf('BiomassCIPerPara'),'jpg')
%}
%{
%Plotting 95% C-I For mean total Biomass
%{
CI = true;
sigLevel = .05;
nullChoice = ZFrees(:,2)+1;

xmin = -.05;
xmax = 1.05;
ymin = -0.05;
ymax = 1.05;

middlePercentage = .5;
lowerQuantile = (1-middlePercentage)/2;
upperQuantile = 1 - lowerQuantile;

AXW = fullfact([2,2]);
ZParaColors = {'r','b'};
ZFreeMarks = {'o','x'};

marksMatrix = {{'k^','MarkerFaceColor','k'},{'x','Color',[.75,0,.75]};
    {'v','MarkerFaceColor',[0,.5,.5]}, {'+','Color',[0 0.75 0]}};

legendEntries = {'(Z_{Free},Z_{Para})=(10,1e-3)'...
    ,'(Z_{Free}),Z_{Para}=(10,1e-4)'...
    ,'(Z_{Free}),Z_{Para}=(100,1e-3)'...
    ,'(Z_{Free}),Z_{Para}=(100,1e-4)'...
    };

fig1 = figure('Position', [0,0,1440,900]);

med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);
y10 = zeros(2,1);
ME10 = zeros(2,1);
ymax = 0;
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            
            y1 = mean(totalBiomassAll{pickYs});
            s1 = std(totalBiomassAll{pickYs});
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(totalBiomassAll{pickYs},0.5);
            ql = abs(quantile(totalBiomassAll{pickYs},lowerQuantile)-med);
            qu = abs(quantile(totalBiomassAll{pickYs},upperQuantile)-med);
            
            
            if numel(med) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                y10(nullChoice(pickYs)) = y1(1);
                ME10(nullChoice(pickYs)) = s1(1);
                
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
                
                y1 = [y10(nullChoice(pickYs)) y1];
                ME = [ME10(nullChoice(pickYs)) ME];
            end
            
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            ymax = max(ymax,max(y1));
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
   ylabel('Totalies,'Location','Best')
    axis([xmin xmax ymin ymax]);
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
    axis([xmin xmax ymin ymax]);
    subplot(2,2,3);
   ylabel('Total Biomas')
    xlabel('Fraction of Consumers as Parasites')
    axis([xmin xmax ymin ymax]);
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    axis([xmin xmax ymin ymax]);
    
end
saveas(fig1,sprintf('BiomasCIAll'),'jpg')

fig2 = figure('Position', [0,0,1440,900]);

med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);
y10 = zeros(2,1);
ME10 = zeros(2,1);
ymax = 0;

for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            
            y1 = mean(totalBiomassFree{pickYs});
            s1 = std(totalBiomassFree{pickYs});
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(totalBiomassFree{pickYs},0.5);
            ql = abs(quantile(totalBiomassFree{pickYs},lowerQuantile)-med);
            qu = abs(quantile(totalBiomassFree{pickYs},upperQuantile)-med);
            
            
            if numel(med) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                y10(nullChoice(pickYs)) = y1(1);
                ME10(nullChoice(pickYs)) = s1(1);
                
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
                
                y1 = [y10(nullChoice(pickYs)) y1];
                ME = [ME10(nullChoice(pickYs)) ME];
            end
            
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            ymax = max(ymax,max(y1));
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
   ylabel('Total Biomas')
    legend(legendEntries,'Location','Best')
    axis([xmin xmax ymin ymax]);
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
    axis([xmin xmax ymin ymax]);
    subplot(2,2,3);
   ylabel('Total Biomas')
    xlabel('Fraction of Consumers as Parasites')
    axis([xmin xmax ymin ymax]);
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    axis([xmin xmax ymin ymax]);
    
end
saveas(fig2,sprintf('BiomassCIPara'),'jpg')

fig3 = figure('Position', [0,0,1440,900]);

med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);
y10 = zeros(2,1);
ME10 = zeros(2,1);
ymax = 0;
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            
            y1 = mean(totalBiomassPara{pickYs});
            s1 = std(totalBiomassPara{pickYs});
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(totalBiomassPara{pickYs},0.5);
            ql = abs(quantile(totalBiomassPara{pickYs},lowerQuantile)-med);
            qu = abs(quantile(totalBiomassPara{pickYs},upperQuantile)-med);
            
            
            if numel(med) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                y10(nullChoice(pickYs)) = y1(1);
                ME10(nullChoice(pickYs)) = s1(1);
                
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
                
                y1 = [y10(nullChoice(pickYs)) y1];
                ME = [ME10(nullChoice(pickYs)) ME];
            end
            
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            ymax = max(ymax,max(y1));
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
   ylabel('Total Biomas')
    legend(legendEntries,'Location','Best')
    axis([xmin xmax ymin ymax]);
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
    axis([xmin xmax ymin ymax]);
    subplot(2,2,3);
   ylabel('Total Biomas')
    xlabel('Fraction of Consumers as Parasites')
    axis([xmin xmax ymin ymax]);
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    axis([xmin xmax ymin ymax]);
    
end
saveas(fig3,sprintf('BiomassCIPara'),'jpg')
%}
%Plotting IQR of PErsistence
%{
CI = false;
sigLevel = .05;
nullChoice = ZFrees(:,2)+1;

xmin = -.05;
xmax = 1.05;
ymin = -0.05;
ymax = 1.05;

middlePercentage = .5;
lowerQuantile = (1-middlePercentage)/2;
upperQuantile = 1 - lowerQuantile;

AXW = fullfact([2,2]);
ZParaColors = {'r','b'};
ZFreeMarks = {'o','x'};

marksMatrix = {{'k^','MarkerFaceColor','k'},{'x','Color',[.75,0,.75]};
    {'v','MarkerFaceColor',[0,.5,.5]}, {'+','Color',[0 0.75 0]}};

legendEntries = {'(Z_{Free},Z_{Para})=(10,1e-3)'...
    ,'(Z_{Free}),Z_{Para}=(10,1e-4)'...
    ,'(Z_{Free}),Z_{Para}=(100,1e-3)'...
    ,'(Z_{Free}),Z_{Para}=(100,1e-4)'...
    };

fig1 = figure('Position', [0,0,1440,900]);

med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);

for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            
            y1 = meansAll(pickYs,:);
            s1 = stdsAll(pickYs,:);
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(persistenceAll{pickYs},0.5);
            ql = abs(quantile(persistenceAll{pickYs},lowerQuantile)-med);
            qu = abs(quantile(persistenceAll{pickYs},upperQuantile)-med);
            
            
            if numel(med) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
            end
            
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            axis([xmin xmax ymin ymax]);
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
    ylabel('Fraction All Species Persistence')
    legend(legendEntries,'Location','Best')
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
    ylabel('Fraction All Species Persistence')
    xlabel('Fraction of Consumers as Parasites')
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    
end
saveas(fig1,sprintf('PersistenceAllQuant'),'jpg')

fig2 = figure('Position', [0,0,1440,900]);
med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);
title('Free Liver Persistence')
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            
            y1 = meansAllFree(pickYs,:);
            s1 = stdsAllFree(pickYs,:);
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(persistenceAllFree{pickYs},0.5);
            ql = abs(quantile(persistenceAllFree{pickYs},lowerQuantile)-med);
            qu = abs(quantile(persistenceAllFree{pickYs},upperQuantile)-med);
            
            
            if numel(ql) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
            end
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            axis([xmin xmax ymin ymax]);
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
    ylabel('Fraction Free Consumer Persistence')
    legend(legendEntries,'Location','Best')
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
    ylabel('Fraction Free-Consumer Persistence')
    xlabel('Fraction of Consumers as Parasites')
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    
end
saveas(fig2,sprintf('PersistenceFreeQuant'),'jpg')

fig3 = figure('Position', [0,0,1440,900]);
med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);
title('Parasite Persistence')
fParAllPara = fParAll(fParAll~=0);
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
             plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
           %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            
            y1 = meansAllPara(pickYs,:);
            s1 = stdsAllPara(pickYs,:);
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(persistenceAllPara{pickYs},0.5);
            ql = abs(quantile(persistenceAllPara{pickYs},lowerQuantile)-med);
            qu = abs(quantile(persistenceAllPara{pickYs},upperQuantile)-med);
            
            
            if numel(ql) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
            end
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            axis([xmin xmax ymin ymax]);
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
    ylabel('Fraction parasite Persistence')
    legend(legendEntries,'Location','Best')
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
    ylabel('Fraction Parasite Persistence')
    xlabel('Fraction of Consumers as Parasites')
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    
end
saveas(fig3,sprintf('PersistenceParaQuant'),'jpg')

CI = true;
sigLevel = .05;
nullChoice = ZFrees(:,2)+1;

middlePercentage = .25;
lowerQuantile = (1-middlePercentage)/2;
upperQuantile = 1 - lowerQuantile;

AXW = fullfact([2,2]);


legendEntries = {'(Z_{Free},Z_{Para})=(10,1e-3)'...
    ,'(Z_{Free}),Z_{Para}=(10,1e-4)'...
    ,'(Z_{Free}),Z_{Para}=(100,1e-3)'...
    ,'(Z_{Free}),Z_{Para}=(100,1e-4)'...
    };
%}
%Plotting Persistence Confidence intervals
%{
fig1 = figure('Position', [0,0,1440,900]);

med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);

for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            
            %plotStyles = marksMatrix(ZParaChoice,ZFreeChoice);
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            
            y1 = meansAll(pickYs,:);
            s1 = stdsAll(pickYs,:);
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
           %Middle x% of Data:
            med = quantile(persistenceAllPara{pickYs},0.5);
            ql = abs(quantile(persistenceAllPara{pickYs},lowerQuantile)-med);
            qu = abs(quantile(persistenceAllPara{pickYs},upperQuantile)-med);
            
            
            if numel(ql) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
            end
            
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            axis([xmin xmax ymin ymax]);
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
    ylabel('Fraction All Species Persistence')
    legend(legendEntries,'Location','Best')
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
    ylabel('Fraction All Species Persistence')
    xlabel('Fraction of Consumers as Parasites')
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    
end
saveas(fig1,sprintf('PersistenceAllCI'),'jpg')

fig2 = figure('Position', [0,0,1440,900]);
med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);
title('Free Liver Persistence')
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
            %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            
            y1 = meansAllFree(pickYs,:);
            s1 = stdsAllFree(pickYs,:);
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(persistenceAllPara{pickYs},0.5);
            ql = abs(quantile(persistenceAllPara{pickYs},lowerQuantile)-med);
            qu = abs(quantile(persistenceAllPara{pickYs},upperQuantile)-med);
            
            
            if numel(ql) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
            end
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            axis([xmin xmax ymin ymax]);
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
    ylabel('Fraction Free Consumer Persistence')
    legend(legendEntries,'Location','Best')
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
    ylabel('Fraction Free-Consumer Persistence')
    xlabel('Fraction of Consumers as Parasites')
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    
end
saveas(fig2,sprintf('PersistenceFreeCI'),'jpg')

fig3 = figure('Position', [0,0,1440,900]);
med0 = zeros(2,1);
ql0 = zeros(2,1);
qu0 = zeros(2,1);
title('Parasite Persistence')
fParAllPara = fParAll(fParAll~=0);
for subplotChoice = 1:4;
    subplot(2,2,subplotChoice);
    hold on
    
    for ZFreeChoice = 1:2;
        for ZParaChoice = 1:2
            %plotStyles = strcat(ZParaColors(ZParaChoice),ZFreeMarks(ZFreeChoice));
            plotStyles = marksMatrix{ZParaChoice,ZFreeChoice};
            pickYs = ZParas(:,ZParaChoice)&...
                freeLivings(:,AXW(subplotChoice,1))&...
                concs(:,AXW(subplotChoice,2))&...
                ZFrees(:,ZFreeChoice);
           %This is a mean - want to give 95% confidence interval for it.
            %Better use a t-distribution with 99 dof.
            
            y1 = meansAllPara(pickYs,:);
            s1 = stdsAllPara(pickYs,:);
            
            %Margin of error:
            tCritical = tinv(1-sigLevel/2,length(y1)-1);
            ME = s1/sqrt(nWeb)*tCritical;
            
            %Middle x% of Data:
            med = quantile(persistenceAllPara{pickYs},0.5);
            ql = abs(quantile(persistenceAllPara{pickYs},lowerQuantile)-med);
            qu = abs(quantile(persistenceAllPara{pickYs},upperQuantile)-med);
            
            
            if numel(ql) == 9
                med0(nullChoice(pickYs)) = med(1);
                ql0(nullChoice(pickYs)) = ql(1);
                qu0(nullChoice(pickYs)) = qu(1);
                
            else
                ql = [ql0(nullChoice(pickYs)) ql];
                qu = [qu0(nullChoice(pickYs)) qu];
                med = [med0(nullChoice(pickYs)) med];
            end
            
            if CI
            %Plotting 95% confidence interval for the mean:
            errorbar(fParAll,y1,ME,plotStyles{:})
            else
            %Plotting center 50% of data for persistence
            errorbar(fParAll,med,ql,qu,plotStyles{:})
            end
            axis([xmin xmax ymin ymax]);
            grid on
        end
    end
    
    hold off
    subplot(2,2,1);
    title(sprintf('No Fraction Free Living Stage'));
    ylabel('Fraction parasite Persistence')
    legend(legendEntries,'Location','Best')
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
    ylabel('Fraction Parasite Persistence')
    xlabel('Fraction of Consumers as Parasites')
    subplot(2,2,4);
    xlabel('Fraction of Consumers as Parasites')
    axesPosition = get(gca,'Position');
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,'Concomittant','FontWeight','bold')
    
end
saveas(fig3,sprintf('PersistenceParaCI'),'jpg')

figure;
hold on
kFree = 1;
kPara = -3;
for ii = 1:nFPar0
    para = paraMxBinAll(:,ii)>0;
    nonBasal = patlAll > 1;
    free = ~para&nonBasal;
    x = zeros(size(para));
    
    mFree = 10.^(kFree.*(patlAll(free)-1));
    mPara = 10.^(kPara-kFree)*10.^(kFree*(patlAll(para)-1));
    
    x(para) = .314*mPara.^(-0.25);
    x(free) = .314*mFree.^(-0.25);
    
    histogram(log10(x),'DisplayStyle','stair');
    
end

legend(cellstr(num2str(fParAll')));

%}
%}

%
% fig = figure('Position', [0,0,1440,900]);
% fntsz = 24;
% hold on
% set(gca,'FontSize',15);
%
%
% title('Persistence of all Species','FontSize',fntsz)
% pickYs = concs(:,1)&...
%     ZFrees(:,1);
% y1 = mean(meansAll(pickYs,:));
% plot([0 fParAll],[mean([nullZ10,nullZ100]),y1],'r-^',...
%     'LineWidth',4,'MarkerSize',15,'MarkerFaceColor','r')
% pickYs = concs(:,2)&...
%     ZFrees(:,1);
% y1 = mean(meansAll(pickYs,:));
% plot([0 fParAll],[mean([nullZ10,nullZ100]),y1],'b--o'f,...
%     'LineWidth',4,'MarkerSize',10,'MarkerFaceColor','b')

%             pickYs = concs(:,1)&...
%                 ZFrees(:,1);
%             y1 = mean(meansAllPara(pickYs,:));
%             plot([0 fParAll],[mean([nullZ10,nullZ100]),y1],'b-x')
%             pickYs = concs(:,2)&...
%                 ZFrees(:,1);
%             y1 = mean(meansAllPara(pickYs,:));
%             plot([0 fParAll],[mean([nullZ10,nullZ100]),y1],'r-x')

% grid on
%
% hold off
% axis([0 0.8 0 1.05]);
% ylabel('Species persistence','FontSize',fntsz)
% xlabel('Initial Fraction of Parasites','FontSize',fntsz)
% leg = legend('No Concomittant Links','Concomittant Links');
% leg.FontSize = 18;
% saveas(fig,sprintf('ConcomittantFigure.jpg'),'jpg')

%}

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
