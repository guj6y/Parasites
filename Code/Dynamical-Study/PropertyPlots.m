% Comparing nuggets and persistent subsets in terms neo understands better

nuggetDir = '~/Desktop/JulySimulations/WebNuggets/';

OGNuggetDir = '~/Desktop/JulySimulations/WebNuggets/OGWebs';

ExperimentsDir = 'Volumes/Buster/Research/Parasitism/Dynamics/JulySimulations';

%Load Z10 and Z100 nuggets



    
nuggetZ10Dir = sprintf('%s/Z%iNuggets',nuggetDir,10);
nuggetZ100Dir = sprintf('%s/Z%iNuggets',nuggetDir,100);

webIndex10 = csvread(sprintf('%s/webIndex.csv',nuggetZ10Dir),1);
webIndex100 = csvread(sprintf('%s/webIndex.csv',nuggetZ100Dir),1);

minS = 20;
maxS = 80;
SList = minS:maxS;
nS = numel(SList);

minC = 0.01;
maxC = 0.25;

nProperties = 15;

S10 = webIndex10(:,5);
C10 = webIndex10(:,6);

used10 = (S10>=minS)&(C10>=minC)&(C10<=maxC);

S100 = webIndex100(:,5);
C100 = webIndex100(:,6);

used100 = (S100>=minS)&(C100>=minC)&(C100<=maxC);

SBoth = [S10(used10);S100(used100)];
CBoth = [C10(used10);C100(used100)];

close all
grpMeans = grpstats(CBoth,SBoth);
errorTemp = (repmat(grpMeans,1,2)-...
    grpstats(CBoth,SBoth,'meanci'))./...
    repmat(grpMeans,1,2);
rangeError = grpstats(CBoth,SBoth,'range')/2./grpMeans;
errors = max(min(errorTemp(:,1),rangeError),.05);

CErrors = zeros(80,1);
CErrors(unique(SBoth)) = errors;
CErrors = max(CErrors,.05);

CErrors = CErrors(20:end);
CAll = zeros(80,1);
CAll(unique(SBoth)) = grpMeans;
CAll = CAll(20:end);
CAll(CAll==0) = mean(CAll(CAll<.16&CAll>0));
CAll(end) = 0.15;

nWebs = 1000;
%Making nmProperties:
try
    load('nmProperties')
catch
nmProperties = zeros(nS,nWebs,nProperties);
tic
for ii = 1:nS
    CError = CErrors(ii);
    S = minS-1+ii;
    C = CAll(ii);

    for jj = 1:nWebs
        [res,con,n,c,r] = NicheModel_nk_Cerr(S,C,CError);
        properties = calculateLocalProperties(res,con,S);
        globalProperties = calculateGlobalProperties(res,con);
        fsp = globalProperties(15);
        fLoop = globalProperties(16);
        
        simMx = calculateSimilarity(res,con);
        
        gen = properties(:,2);
        vul = properties(:,3);
        swtl = properties(:,10);
        sptb = properties(:,8) +1;
        patl = 2*swtl - sptb;
        
        fTop = sum(vul==0)/S;
        fBasal = sum(gen==0)/S;
        fInt = 1-fTop-fBasal;
        
        fCann = sum(res==con)/S;
        fHerb = sum(patl==2)/S;
        fOmn = sum(mod(patl,1)>0)/S;
        
        TL = mean(patl);
        
        MaxSim = mean(max(simMx));
        density = mean(vul);
        vulSD = std(vul);
        genSD = std(gen);
        A = sparse(res,con,1,S,S);
        ccs = clustering_coefficients(A);
        clust1 = mean(ccs);
        
        clust2 = (trace((A)^3)- trace(A) - trace(A^2))/(sum(sum(A^2))-trace(A^2));
        nmProperties(ii,jj,:) = [fTop,fInt,fBasal,fCann,fHerb,fOmn,...
            TL,MaxSim,vulSD,genSD,density,clust1,clust2,fsp,fLoop];
        
    end
    
    
end
toc
save('nmProperties','nmProperties')
end
%NB: Consider saving swtl and patl.  also, other properties?  Fraction of
%shortest paths?
nmProperties(:,:,6) = nmProperties(:,:,6)./repmat((20:80)',1,1000);

propertyNames = {'fTop';...     1
                 'fInt';...     2
                 'fBasal';...   3
                 'fCann';...    4
                 'fHerb';...    5
                 'fOmn';...     6
                 'TL';...       7
                 'MaxSim';...   8
                 'vulSD';...    9
                 'genSD';...    10
                 'density';...  11
                 'clust1';...   12
                 'clust2';...   13
                 'fsp';...      14
                 'fLoop'};...    15
 
alpha = .025;
tcrit = tinv(1-alpha,nWebs-1);
meanProperties = mean(nmProperties,2);
stdProperties = std(nmProperties,0,2);

meanProperties = reshape(meanProperties,61,nProperties);
stdProperties = reshape(stdProperties,61,nProperties);


per025 = quantile(nmProperties,.025,2);
per975 = quantile(nmProperties,.975,2);

per025 = reshape(per025,61,[]);
per975 = reshape(per975,61,[]);

per025Scaled = (per025 -meanProperties)./stdProperties;
per975Scaled = (per975 -meanProperties)./stdProperties;

meProperties = tcrit*stdProperties/sqrt(nWebs);
u = meanProperties + meProperties;
l =  meanProperties - meProperties;

%SOmething clearly! wrong here.  fBasal should be high, not low.
try 
    load('nuggetProperties')
catch
nNuggets = length(webIndex10(:,1));

nuggetProperties = zeros(nNuggets,nProperties,2);
for k = [1,2]
    Z = 10^k;
    ZWebDir = sprintf('%s/Z%uNuggets',dataDir,Z);
    ZWebIndex = csvread(sprintf('%s/webIndex.csv',ZWebDir),1,0);
    webNumbersZ = ZWebIndex(:,1);    
    
    count = 0;
    
    for web = webNumbersZ'


        count = count + 1;
        S0 = ZWebIndex(count,2);
        S = ZWebIndex(count,5);
        C = ZWebIndex(count,6);
        if (S < 20)||(C<.01)||(C>.25)
            continue
        end

        sListOld = 1:S0;
        sListNew = 1:S;
        sRename = sListOld;
        nExtinct =  S0-S;
        
        OGNodeData = csvread(sprintf('%s/nodeData%06u.csv',OGWebDir,web));
        nodeData = csvread(sprintf('%s/nodeData%06u.csv',ZWebDir,web));
        linkData = csvread(sprintf('%s/linkData%06u.csv',ZWebDir,web));
        
        
        
        gr0 = OGNodeData(:,4);
        n0 = OGNodeData(:,1);
        c0 = OGNodeData(:,2);
        r0 = OGNodeData(:,3);
        
        alive = isinf(nodeData(:,5));
        
        n = n0(alive);
        c = c0(alive);
        r = r0(alive);
        
        res0 = linkData(:,1);
        con0 = linkData(:,2);
        L = sum(~linkData(:,3));
        S = numel(n);

        try
        [res,con,~,~,~,~] = AddTrophicConsumers(n,c,r,0,C);
        catch
            continue
        end
        
        properties = calculateLocalProperties(res,con,S);
        globalProperties = calculateGlobalProperties(res,con);
        fsp = globalProperties(15);
        fLoop = globalProperties(16);
        
        simMx = calculateSimilarity(res,con);
        
        gen = properties(:,2);
        vul = properties(:,3);
        swtl = properties(:,10);
        sptb = properties(:,8) +1;
        patl = 2*swtl - sptb;
        
        fTop = sum(vul==0)/S;
        fBasal = sum(gen==0)/S;
        fInt = 1-fTop-fBasal;
        
        fCann = sum(res==con)/S;
        fHerb = sum(patl==2)/S;
        fOmn = sum(mod(patl,1)>0)/S;
        
        TL = mean(patl);
        
        MaxSim = mean(max(simMx));
        density = mean(vul);
        vulSD = std(vul);
        genSD = std(gen);
        A = sparse(res,con,1,S,S);
        ccs = clustering_coefficients(A);
        clust1 = mean(ccs);
        
        clust2 = (trace((A)^3)- trace(A) - trace(A^2))/(sum(sum(A^2))-trace(A^2));
        nuggetProperties(count,:,k) = [fTop,fInt,fBasal,fCann,fHerb,fOmn,...
            TL,MaxSim,vulSD,genSD,density,clust1,clust2,fsp,fLoop];
        
        
    end
end
save('nuggetProperties','nuggetProperties')
end



S10 = round(webIndex10(:,5));
S100 = round(webIndex100(:,5));


meanPropertiesNuggets10 = grpstats(nuggetProperties(:,:,1),S10);
meanPropertiesNuggets100 = grpstats(nuggetProperties(:,:,2),S100);

nuggetMeans10 = nan(80,nProperties);
nuggetMeans100 = nan(80,nProperties);

nuggetMeans10(unique(S10),:) = meanPropertiesNuggets10;
nuggetMeans100(unique(S100),:) = meanPropertiesNuggets100;

nuggetMeans10 = nuggetMeans10(20:80,:);
nuggetMeans100 = nuggetMeans100(20:80,:);


ciPropertiesNuggets10 = grpstats(nuggetProperties(:,:,1),S10,'meanci','Alpha',.4);
ciPropertiesNuggets100 = grpstats(nuggetProperties(:,:,2),S100,'meanci','Alpha',.4);

ciNuggets10 = nan(80,nProperties,2);
ciNuggets100 = nan(80,nProperties,2);

ciNuggets10(unique(S10),:,:) = ciPropertiesNuggets10;
ciNuggets100(unique(S100),:,:) = ciPropertiesNuggets100;

ciNuggets10 = ciNuggets10(20:80,:,:);
ciNuggets100 = ciNuggets100(20:80,:,:);

nuggetMeans10Scaled = (nuggetMeans10 - meanProperties)./stdProperties;
nuggetMeans100Scaled = (nuggetMeans100 - meanProperties)./stdProperties;

ciL10 = abs(ciNuggets10(:,:,1)-nuggetMeans10)./stdProperties;
ciU10 = abs(ciNuggets10(:,:,2)-nuggetMeans10)./stdProperties;
ciL100 = abs(ciNuggets100(:,:,1)-nuggetMeans100)./stdProperties;
ciU100 = abs(ciNuggets100(:,:,2)-nuggetMeans100)./stdProperties;


sList = 1:80;
try
    load('deletedNMProperties')
catch
deletedNMProperties = zeros(nS,nWebs,nProperties);
tic
for ii = 1:nS

    S = minS-1+ii
    nToDelete = 80-S;
    C = CAll(ii);
jj = 0;
nGood = 0;
    while nGood < nWebs
        jj
        [res,con,n,c,r] = NicheModel_nk(80,.15);
        basal = true(80,1);
        basal(con) = false;
        cons = sList(~basal);
        
        
        if nToDelete > length(cons)
            continue
        end
        damned = datasample(cons,nToDelete,'Replace',false);
        blessed = true(80,1);
        blessed(damned) = false;
        try
        [res,con,~,~,~,~] = AddTrophicConsumers(n(blessed),c(blessed),r(blessed),0,.15);
        catch
            continue
        end
        try
        properties = calculateLocalProperties(res,con,S);
        globalProperties = calculateGlobalProperties(res,con);
        catch
            jj = jj-1;
            nGood = nGood -1;
            continue
        end
        fsp = globalProperties(15);
        fLoop = globalProperties(16);
        
        simMx = calculateSimilarity(res,con);
        
        gen = properties(:,2);

        vul = properties(:,3);
        swtl = properties(:,10);
        sptb = properties(:,8) +1;
        patl = 2*swtl - sptb;
        
        fTop = sum(vul==0)/S;
        fBasal = sum(gen==0)/S;
        fInt = 1-fTop-fBasal;
        
        fCann = sum(res==con)/S;
        fHerb = sum(patl==2)/S;
        fOmn = sum(mod(patl,1)>0)/S;
        
        TL = mean(patl);
        
        MaxSim = mean(max(simMx));
        density = mean(vul);
        vulSD = std(vul);
        genSD = std(gen);
        A = sparse(res,con,1,S,S);
        ccs = clustering_coefficients(A);
        clust1 = mean(ccs);
        
        clust2 = (trace((A)^3)- trace(A) - trace(A^2))/(sum(sum(A^2))-trace(A^2));
        jj = jj+1;
            nGood = nGood+1;
        deletedNMProperties(ii,jj,:) = [fTop,fInt,fBasal,fCann,fHerb,fOmn,...
            TL,MaxSim,vulSD,genSD,density,clust1,clust2,fsp,fLoop];
    end
    
    
end
toc
save('deletedNMProperties','deletedNMProperties')
end

deletedMeans = mean(deletedNMProperties,2);
deletedMeans = reshape(deletedMeans,61,[]);

alpha = .025;
tcrit = tinv(1-alpha,nWebs-1);

deletedSTD = std(deletedNMProperties,0,2);
deletedSTD = reshape(deletedSTD,61,[]);




per025Deleted = quantile(nmProperties,.025,2);
per975Deleted = quantile(nmProperties,.975,2);

per025Deleted = reshape(per025Deleted,61,[]);
per975Deleted = reshape(per975Deleted,61,[]);

per025DeletedScaled = (per025Deleted -meanProperties)./stdProperties;
per975DeletedScaled = (per975Deleted -meanProperties)./stdProperties;

mePropertiesDeleted = tcrit*deletedSTD/sqrt(nWebs);
uDeleted = mePropertiesDeleted;
lDeleted = mePropertiesDeleted;


deletedMeansScaled = (deletedMeans - meanProperties)./stdProperties;

uDeletedScaled = (uDeleted)./stdProperties;
lDeletedScaled = (lDeleted)./stdProperties;

fractionalPropeties = [1:6,8,12,13,14,15];

for plotProp = 1:nProperties;
    figure
    %errorbar(20:80,meanProperties(:,:,plotProp),meProperties(:,:,plotProp))
    hold on
    plot(20:80,meanProperties(:,plotProp)-meanProperties(:,plotProp),'k');
    errorbar(20:80,nuggetMeans10Scaled(:,plotProp),ciL10(:,plotProp),ciU10(:,plotProp),'o')
    errorbar(20:80,nuggetMeans100Scaled(:,plotProp),ciL100(:,plotProp),ciU100(:,plotProp),'rx')
    
    errorbar(20:80,deletedMeansScaled(:,plotProp),lDeletedScaled(:,plotProp),uDeletedScaled(:,plotProp),'gp')
    
    
    plot(20:80,per025Scaled(:,plotProp),'k--')
    plot(20:80,per975Scaled(:,plotProp),'k--')
    YL = ylim;
    YL(1) = max(YL(1),-20);
    YL(2) = min(YL(2),20);
    plot(20:80,-meanProperties(:,plotProp)./stdProperties(:,plotProp),'r');
    plot(20:80,(1-meanProperties(:,plotProp))./stdProperties(:,plotProp),'r');
    xlim([20,80]);
    ylim(YL);
    hold off
    title(propertyNames{plotProp})
    
    

end



