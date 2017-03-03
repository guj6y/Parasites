%Peeking at the data from the parasitic/free-liver invasions.

clear all
fclose('all');
clc

%Starting Values
S0 = 80;
C0 = 0.15;

%numer of webs plus/minus of the median persistence.
nInvaders = 100;
ZAll = [10,100];
nZ = length(ZAll);

ratios = [-3,-4,1];
nRatios = numel(ratios);

models = fullfact([2,2]);
nModels = length(models);
dataDir = '~/Desktop/JulySimulations/WebNuggets';

OGWebDir = sprintf('%s/OGWebs',dataDir);
codeDir = pwd;

progressData = dlmread(...
    sprintf('%s/ExperimentProgress.txt',dataDir));
ZStart = 1;
webStart = 1;
kInvaderStart = 1;
modelStart = 1;

webEnd = progressData(2)-1;
ZEnd = progressData(1);
invaderEnd = nInvaders;
modelEnd = nModels;
kInvaderEnd = nRatios;

persistenceAll = zeros(webEnd,nModels,nRatios,nInvaders);

phi = 0.15;

for Zidx = ZStart:ZEnd
    Z = ZAll(Zidx);
    kFree = log10(Z);
    webDir = sprintf('%s/Z%uNuggets',dataDir,Z);
    invasionDir = sprintf('%s/Z%uInvasions',dataDir,Z);
    
    try
        cd(invasionDir)
        cd(codeDir)
    catch
        mkdir(invasionDir)
    end
    
    webIndex = csvread(sprintf('%s/webIndex.csv',webDir),1,0);
    
    nWeb = numel(webIndex(:,1));
    
    webNumberUsed = csvread(sprintf('%s/webNumbersUsed.csv',invasionDir));
    
    count = 0;
    
    for kWeb = webStart:webEnd
        web = webNumberUsed(kWeb);
        kWeb
        count = count + 1;
        
        S0 = webIndex(web+1,2);
        S = webIndex(web+1,5);
        
        sListNew = 1:S;
        sRename = zeros(S0,1);
        
        
        %Need to load og web to calculate the appropriate x values
        OGNodeData = csvread(sprintf('%s/nodeData%06u.csv',OGWebDir,web));
        linkData = csvread(sprintf('%s/linkData%06u.csv',OGWebDir,web));
        
        
        
        res0 = linkData(:,1);
        con0 = linkData(:,2);
        C = length(res0)/S0^2;
        
        
        localProperties0 = calculateLocalProperties(res0,con0);
        swtl0 = localProperties0(:,10);
        sptb0 = localProperties0(:,8)+1;
        patl0 = 2*swtl0 - sptb0;
        
        basal0 = swtl0 == 1;
        nBasal0 = sum(basal0);
        
        %Now load the filtered web data to remove dead things.
        nodeData = csvread(sprintf('%s/NodeData%06u.csv',webDir,web));
        
        %Find the extinct Species
        extctSp = nodeData(:,end)<inf;
        extntSp = ~extctSp;
        
        %Note that we don't need the new link list data since we need to
        %add species anyway.  These aren't the final n,c,r,gr vectors since
        %I still need to add an invader.
        n1 = OGNodeData(extntSp,1);
        c1 = OGNodeData(extntSp,2);
        r1 = OGNodeData(extntSp,3);
        
        %Growth rates
        gr = [OGNodeData(extntSp,4);0];
        
        basal = [basal0(extntSp);0]>0;
        
        %Prey-averaged trophic levels of living species; keep these from OG
        %web so that the metabolic rates don't change.
        patl1 = patl0(extntSp);
        
        invaderNiches = csvread(...
            sprintf('%s/web%uInvaderNicheParameters.txt',...
            invasionDir,web));
        for jj = 1:nInvaders
            jjInvaderNiche = invaderNiches(:,1);
            nInv = [n1;jjInvaderNiche(1)];
            cInv = [c1;jjInvaderNiche(2)];
            rInv = [r1;jjInvaderNiche(3)];
            
            [resInv,conInv] = AddTrophicConsumers(nInv,cInv,rInv,0,C0);
            localProperties = calculateLocalProperties(resInv,conInv);
            
        end
        
        
        
        for modelIdx = modelStart:modelEnd;
            model = models(modelIdx,:);
            
            param.modelCode = model;
            modelDir = sprintf('%s/model %u%u',invasionDir,model);
            
            try
                cd(modelDir)
                cd(codeDir)
            catch
                mkdir(modelDir)
            end
            count2 = 0;
            for kInvaderIdx = kInvaderStart:kInvaderEnd
                kInvader = ratios(kInvaderIdx);
                
                if (kInvader == 1) && modelIdx>=2
                    continue
                end
                
                
                invaderZDir = sprintf('%s/kInvader%i',modelDir,kInvader);
                
                
                try
                    cd(invaderZDir)
                    cd(codeDir)
                catch
                    mkdir(invaderZDir)
                end
                
                
                means = csvread(sprintf('%s/meansWeb%u.csv',invaderZDir,web));
                finals= csvread(sprintf('%s/finalsWeb%u.csv',invaderZDir,web));
                stds = csvread(sprintf('%s/stdsWeb%u.csv',invaderZDir,web));
                    
                finals = finals(1:nInvaders,:);
                persistenceAll(kWeb,modelIdx,kInvaderIdx,:) = mean(finals>10^-12,2);
                
                
            end
        end
    end
end

%4x3 matrix with persistence for the different invasion scenarios
meanPersistenceOverall = mean(mean(persistenceAll),4);
%50x4x3 array with mean persistence by web for each scenario; think of it
%as 50 (4x3) matricesstacked up.
meanPersistenceByWeb = squeeze(mean(persistenceAll,4));

%88% quantiles for persistence of each web (for 50 webs, it's
%basically the middle 44 (3 outliers discarded).  We could probably
%calculate prediction intervals, but I believe that is a little risky/
%biased towards data.  can we assume these are normally distributed?
uQuantile = quantile(persistenceAll,.06,4);
lQuantile = quantile(persistenceAll,.94,4);

persistenceSortedByInvader = sort(persistenceAll,4);
meanPersistenceSortedByWeb = sort(meanPersistenceByWeb);

