%This code adds a single species to a subset of my web nuggets; I'm testing
%the invasibility of parasites vs. free-livers in stable web nuggets.  I am
%replicating the experiment 100 times for each web, for each model.  I am
%using the middle 20% of persistence webs (median +/- 10%).  Can we predict
%which webs are more likely to be invaded?  The invaders start at 1/100 of
%the smallest biomass in the system.  The invaders are defined to be
%extinct if they pass below min(extinction threshold,IC/1000).  Simulation
%end time is
% min(5000,t(invader extinction)) if invader first extinction
% max(5000,t(last extinction)+2000) if invader ~first extinction
%(Assume that the stable states are fairly stable to IC- not a bad
%approximation, but worth systematic verification).

%rng(1);

clear all
fclose('all');
clc
%Just a reminder..... ;)
fprintf('DO NOT CTRL-C TO END SIMULATION.  USE CANCEL ON THE WAITBAR\n')
fprintf('DO NOT CTRL-C TO END SIMULATION.  USE CANCEL ON THE WAITBAR\n')
fprintf('DO NOT CTRL-C TO END SIMULATION.  USE CANCEL ON THE WAITBAR\n')
fprintf('DO NOT CTRL-C TO END SIMULATION.  USE CANCEL ON THE WAITBAR\n')
fprintf('DO NOT CTRL-C TO END SIMULATION.  USE CANCEL ON THE WAITBAR\n')
fprintf('DO NOT CTRL-C TO END SIMULATION.  USE CANCEL ON THE WAITBAR\n')
fprintf('DO NOT CTRL-C TO END SIMULATION.  USE CANCEL ON THE WAITBAR\n')
fprintf('DO NOT CTRL-C TO END SIMULATION.  USE CANCEL ON THE WAITBAR\n')
fprintf('DO NOT CTRL-C TO END SIMULATION.  USE CANCEL ON THE WAITBAR\n')
fprintf('DO NOT CTRL-C TO END SIMULATION.  USE CANCEL ON THE WAITBAR\n')

%Starting Values
S0 = 80;
C0 = 0.15;

%numer of webs plus/minus of the median persistence.
pmWebs = 50;
nWebs = 2*pmWebs + 1;
nInvaders = 100;
ZAll = [10,100];
nZ = length(ZAll);

ratios = [-3,-4,1];
nRatios = numel(ratios);

models = fullfact([2,2]);
nModels = length(models);
dataDir = '~/Desktop/JulySimulations/WebNuggets';
totalIters = nZ*nWebs*nInvaders*nModels*(nRatios-1) + ...
    nZ*nWebs*nInvaders;



waitbarHandle = waitbar(0,sprintf('Setting Up\n\n\n\n\n'),'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(waitbarHandle,'canceling',0);



OGWebDir = sprintf('%s/OGWebs',dataDir);
codeDir = pwd;

try
    progressData = dlmread(...
        sprintf('%s/ExperimentProgress.txt',dataDir));
    ZStart = progressData(1);
    webStart = progressData(2);
    invaderStart = progressData(3);
    kInvaderStart = progressData(4);
    modelStart = progressData(5);
    countAll = progressData(6);
    
catch ME
    
    ZStart = 1;
    webStart = 1;
    invaderStart = 1;
    modelStart = 1;
    kInvaderStart = 1;
    countAll = 0;
end
ZEnd = nZ;
webEnd = nWebs;
invaderEnd = nInvaders;
modelEnd = nModels;
kInvaderEnd = nRatios;


try
    cd(dataDir)
    cd(codeDir)
    
catch
    mkdir(dataDir)
end

try
    cd(OGWebDir)
    cd(codeDir)
    
catch
    mkdir(OGWebDir)
end

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
    webIndexSorted = sortrows(webIndex,5);
    
    medWebIdx = ceil(nWeb/2);
    minWebIdxUsed = medWebIdx-pmWebs;
    maxWebIdxUsed = medWebIdx+pmWebs;
    
    webNumberUsed = webIndexSorted(minWebIdxUsed:maxWebIdxUsed,1);
    
    
    csvwrite(sprintf('%s/webNumbersUsed.csv',invasionDir),webNumberUsed);
    
    
    count = 0;
    
    
    for kWeb = webStart:webEnd
        web = webNumberUsed(kWeb);
        
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
        
        %Initial (stable)conditions; these are the final values not the
        %mean of the final state; good since averaging a periodic solution
        %could destabilize it!
        IC1 = nodeData(extntSp,4);
        B0 = [IC1;1e-10];
        
        %using invertebrates, for now.
        ax = 0.314;
        
        %Carrying capacity of 5
        K = 5;
        
        %Using Type II.2 response
        h = 1.2;
        
        %half-saturation
        halfSat = 0.5;
        
        %Minimum Final Time
        Tf = 5000;
        
        
        for invader = invaderStart:invaderEnd
            try
                nInvader = progressData(7);
                cInvader = progressData(8);
                rInvader = progressData(9);
                [res,con,n,c,r,~] = AddTrophicConsumers(...
                    [n1;nInvader],[c1;cInvader],[r1;rInvader],0,C0);
                clear progressData
            catch
                [res,con,n,c,r,~] = AddTrophicConsumers(n1,c1,r1,1,C0);
                nicheParamInvader = [n(end) c(end) r(end)];
                invaderNicheFid = fopen(sprintf('%s/web%uInvaderNicheParameters.txt',invasionDir,web),'a');
                fprintf(invaderNicheFid,'%.9e,%.9e,%.9e\n',nicheParamInvader);
                fclose(invaderNicheFid);
            end
            localProperties = calculateLocalProperties(res,con);
            swtl = localProperties(:,10);
            sptb = localProperties(:,8)+1;
            patl = 2*swtl - sptb;
            
            %Assimilation efficiency is given in Brose et al as .45 and .85 for
            %herbivores and carnivores, respectively.  This would correspond to
            %the links we talkin bout.  Easy enough to figure:
            eij = zeros(size(res));
            
            %Future: add an eij for parasitic links
            eij(basal(res)) = .45;
            eij(~basal(res)) = .85;
            
            
            
            %Strong generalist model for preference:
            wij = ones(size(res));
            
            
            %Makes all basal species have a body size = 1; more or less.
            M1 = Z.^(patl1 - 1);
            
            
            
            
            y = 8*ones(S+1,1);
            
            yij = y(con);
            
            param = struct( 'S',S+1 ...
                ,'C',C0...
                ,'res',res...
                ,'con',con...
                ,'K',K...
                ,'eij',eij...
                ,'wij',wij...
                ,'basal',basal...
                ,'r',gr...
                ,'B0',B0...
                ,'h',h...
                ,'x',zeros(S+1,1)...
                ,'halfSat',halfSat...
                ,'Tf',Tf ...
                ,'extctThresh',1e-15...
                ,'AbsTol',1e-10...
                ,'RelTol',1e-7...
                ,'yij',yij...
                ,'nInvaders',1 ...
                ...,'modelCode',model ...
                );
            results = zeros(S+1,5);
            results(:,1) = 1:(S+1);
            
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
                    
                    countAll = countAll + 1;
                    waitbarMsg = sprintf('Integrating...\n');
                    waitbarMsg = sprintf('%sNuggets with Z = %u\n',waitbarMsg,Z);
                    waitbarMsg = sprintf('%sWeb %u out of %u\n',...
                        waitbarMsg,kWeb,length(webNumberUsed));
                    waitbarMsg = sprintf('%sInvader Number %u\n',...
                        waitbarMsg,invader);
                    waitbarMsg = sprintf('%sModel code %u;%u\n',waitbarMsg,model);
                    waitbarMsg = sprintf('%sInvader Ratio %i',waitbarMsg,kInvader);
                    waitbar(countAll/totalIters,waitbarHandle,waitbarMsg);
                    
                    
                    invaderZDir = sprintf('%s/kInvader%i',modelDir,kInvader);
                    
                    
                    try
                        cd(invaderZDir)
                        cd(codeDir)
                    catch
                        mkdir(invaderZDir)
                    end
                    
                    
                    meansFid = fopen(sprintf('%s/meansWeb%u.csv',invaderZDir,web),'a');
                    finalsFid = fopen(sprintf('%s/finalsWeb%u.csv',invaderZDir,web),'a');
                    stdsFid = fopen(sprintf('%s/stdsWeb%u.csv',invaderZDir,web),'a');
                    
                    if kInvader <0
                        para = [zeros(S,1);1];
                        param.para = para;
                        param.phi = para*phi;
                    else
                        param.para = zeros(S+1,1);
                        param.phi = zeros(S+1,1);
                    end
                    
                    M = [M1;10.^(kInvader + kFree*(patl(end)-2))];
                    param.M = M;
                    x = ax*M.^(-0.25);
                    x(basal) = 0;
                    
                    param.x = x;
                    
                    %             fprintf('Working on Web %u (%u out of %u)\n',web,count,...
                    %                         length(webNumberUsed));
                    %             fprintf('Integrating invader %u\n',invader)
                    %             fprintf('Integrating Model%u%u\n',model);
                    %             fprintf('Integrating Invader %u with Z=%i\n',invader,kInvader)
                    
                    sol = integrateParasiteExperiments(param);
                    
                    fprintf(meansFid,'%9e,',sol.mean);
                    fprintf(finalsFid,'%.9e,',sol.y(:,end));
                    fprintf(stdsFid,'%.9e,',sol.std);
                    
                    fprintf(meansFid,'\n');
                    fprintf(finalsFid,'\n');
                    fprintf(stdsFid,'\n');
                    
                    fclose(meansFid);
                    fclose(finalsFid);
                    fclose(stdsFid);
                    if getappdata(waitbarHandle,'canceling')==1
                        break
                    end
                end
                
                if kInvaderIdx==nRatios
                    kInvaderidx = 1;
                    kInvaderStart = 1;
                    kInvaderReset = true;
                else
                    kInvaderIdx = kInvaderIdx+1;
                    kInvaderStart = 1;
                end
                
                if getappdata(waitbarHandle,'canceling')==1
                    break
                end
            end
            
            if (modelIdx == nModels)&&(kInvaderReset)
                modelReset = true;
                modelIdx = 1;
                modelStart = 1;
            end
            
            if getappdata(waitbarHandle,'canceling')==1
                break
            end
        end
        if (invader == nInvaders)&&modelReset
            invaderReset = true;
            invader = 1;
            invaderStart = 1;
            nicheParamInvader = [];
        end
        
        if getappdata(waitbarHandle,'canceling')==1
            break
        end
    end
    if (kWeb == nWebs)&&(invaderReset)
        kWeb = 1;
        webStart = 1;
        nicheParamInvader = [];
    end
    
    if getappdata(waitbarHandle,'canceling')==1
        break
    end
end

if getappdata(waitbarHandle,'canceling')==1
    

    
    expProgFid = ...
        fopen(sprintf('%s/ExperimentProgress.txt',dataDir),'w');
    fprintf(expProgFid,'%u\n%u\n%u\n%u\n%u\n%u\n%.9e\n%.9e\n%.9e',...
        Zidx,kWeb,invader,kInvaderIdx,modelIdx,countAll,nicheParamInvader);
    fclose(expProgFid);
end
delete(waitbarHandle);