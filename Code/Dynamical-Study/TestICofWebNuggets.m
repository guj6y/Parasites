

clear all
fclose('all');

%Starting Values
S0 = 80;
C0 = 0.15;

%We can test on nuggets or not.
useNuggets = false;

nICCode = 2;

%Three simulations.
%1. 100 webs, 10 IC each
%2. 10 webs, 100 IC each
%3. 1 web, 1000 IC

%1: Pick a good sample. of webs (i.e. just choose randomly)
%2: Again, pick a good sample but fixed, i.e.: 90,80,70,60,...
%percentiles.
%3. Pick a web at the median.


dataDir = '~/Desktop/JulySimulations/WebNuggets';
OGWebDir = sprintf('%s/OGWebs',dataDir);

codeDir = pwd;


for Z = [10,100]
    ZWebDir = sprintf('%s/Z%uNuggets',dataDir,Z);
    ZWebIndex = csvread(sprintf('%s/webIndex.csv',ZWebDir),1,0);
    webNumbersZ = ZWebIndex(:,5);
    [~,websSorted] = sort(webNumbersZ);
    
    ZWebDir = sprintf('%s/Z%uNuggets',dataDir,Z);
    ICTestDir = sprintf('%s/ICTestsZ%uICCode%u',dataDir,Z,nICCode);
    
    try
        cd(ICTestDir);
        cd(codeDir);
    catch
        mkdir(ICTestDir)
    end
    
    if nICCode == 1
        websUsed = websSorted(linspace(100,900,81));
        nWebs = 81;
        nIC = 10;
    elseif nICCode == 2
        websUsed = websSorted(linspace(101,901,9));
        nWebs = 9;
        nIC = 100;
    else
        websUsed = websSorted(502);
        nWebs = 1;
        nIC = 1000;
    end
    
    csvwrite(sprintf('%s/webNumbers.csv',ICTestDir),websUsed);
    
    deadWebs = zeros(size(webNumbersZ));
    count = 0;
    
    
    
    for web = websUsed'
        %Loading web.
        
        count = count + 1;
        fprintf('Working on Web %u (%u out of %u)\n',web,count,...
            length(websUsed));
        
        nodeData = csvread(sprintf('%s/nodeData%06u.csv',ZWebDir,web));
        nodeDataOG = csvread(sprintf('%s/nodeData%06u.csv',OGWebDir,web));
        
        n = nodeData(:,1);
        c = nodeData(:,2);
        r = nodeData(:,3);
        if useNuggets
            extnt = isinf(nodeData(:,5));
            [res0,con0,~,~,~] = AddTrophicConsumers(n,c,r,0,.15);
            
            gr = nodeDataOG(extnt,4);
            localPropertiesOG = calculateLocalProperties(res0,con0);
            localPropertiesOG = localPropertiesOG(extnt,:);
            
            swtl = localPropertiesOG(:,10);
            sptb = localPropertiesOG(:,8) + 1;
            patl = 2*swtl-sptb;
            
            n = n(extnt);
            c = c(extnt);
            r = r(extnt);
            
            [res,con,~,~,~] = AddTrophicConsumers(n,c,r,0,.15);
            S = sum(extnt);
        else
            [res,con,~,~,~] = AddTrophicConsumers(n,c,r,0,.15);
            gr = nodeDataOG(:,4);
            
            localProperties = calculateLocalProperties(res,con);
            
            swtl = localProperties(:,10);
            sptb = localProperties(:,8)+1;
            patl = 2*swtl - sptb;
            
            S = 80;
        
        end
        
        basal = swtl == 1;
        nBasal = sum(basal);
        
        %using invertebrates, for now.
        ax = 0.314;
        
        %Assimilation efficiency is given in Brose et al as .45 and .85 for
        %herbivores and carnivores, respectively.  This would correspond to
        %the links we talkin bout.  Easy enough to figure:
        eij = zeros(size(res));
        
        %Future: add an eij for parasitic links
        eij(basal(res)) = .45;
        eij(~basal(res)) = .85;
        
        %Carrying capacity of 5
        K = 5;
        
        %Strong generalist model for preference:
        wij = ones(size(res));
        
        %Using Type II.2 response
        h = 1.2;
        
        %half-saturation
        halfSat = 0.5;
        
        %Minimum Final Time
        Tf = 10000;
        
        %Initial Conditions
        B0mx = .95*rand(S0,nIC)+.05;
        
        csvwrite(sprintf('%s/%uICWeb%uB0mx.csv',ICTestDir,nIC,web),B0mx);
        
        %Growth rate of basal species.  Not sure if this variability is
        %reasonable or not.
        %Growth rate is loaded up above.
        
        %Makes all basal species have a body size = 1; more or less.
        M = Z.^(patl - 1);
        
        %Prey-averaged trophic level does pretty well at reproducing
        %empirical body size ratio distributions.
        
        x = ax*M.^(-0.25);
        x(basal) = 0;
        y = 8*ones(S,1);
        
        yij = y(con);
        param = struct( 'S',S0...
            ,'C',C0...
            ,'res',res...
            ,'con',con...
            ,'K',K...
            ,'swtl',swtl...
            ,'eij',eij...
            ,'wij',wij...
            ,'basal',basal...
            ,'r',gr...
            ,'B0',zeros(S,1)...
            ,'h',h...
            ,'x',x...
            ,'halfSat',halfSat...
            ,'Tf',Tf ...
            ,'extctThresh',1e-10...
            ,'AbsTol',1e-10...
            ,'RelTol',1e-7...
            ,'yij',yij...
            );
        fcMx = zeros(S,nIC);
        for kk = 1:nIC
            fprintf('Integrating IC %02u\n',kk)
            param.B0 = B0mx(:,kk);
            
            sol = integrateToFilterWebs(res,con,param);
             
            fcMx(:,kk) = sol.y(:,end);
            
        end
        csvwrite(sprintf('%s/%uICWeb%uResults.csv',ICTestDir,nIC,web),fcMx);
        clc
    end
    
end


