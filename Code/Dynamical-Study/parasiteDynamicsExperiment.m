%Script for testing parasite inclusion in niche webs.  In the spirit of
%BWM2006.
%
% Experimental Design (goal/full): 1. Z_f \in {10^1,10^2,10^3} 2. Z_p \in
% {10^-1,10^-2,10^-3} 3.   Metabolic classes: (p,f) = (I,I); (V,V); (I,V)
% 4. Model includes Fraction of time free-living 5. Model includes special
% scaling (Hechinger 2012/2013); lower y? 6. Concomittant Predation
%
% With 100 webs, this gives 3x3x3x2x2x2 = 216x100 simulations.  We also
% need to take into account the fraction of species that are parasites;
% that will be the x-axis for any plots that I make.  I would like at least
% as much resolution as in BWM 2006, which is 8.  Thus, we have: 216 x 100
% x 8 = 172800 simulation runs.  I want to go out to Tf = 10000. The Null
% model something like 8-12 hours (? I don't remember and don't have it
% written down >.<) with 3x2x100= 600 simulations.. so I should expect this
% to tak 288*(8 or 12) = 2304 or 3456 hours =  96 or 144 days. If the null
% model only took one hour, I hsould expect it to talke 288/24 = 12 days
%NB
% Might be necessary to get this running on a larger computer.  Or temper
% my experiment slightly: Maybe: 1. Z_f and Z_p tied to metabolic class.
% decrease time by a factor of 9, so 11-16 days to run.  My BWM2006 code
% took 23 hours to run, they had 8*36*100 = 28800 model runs.  Tf was 2000,
% and the web size varied so.. that doesn't help things.  reduced
% experiment would have 19200 model runs.  Increasing time to runby a
% factor of 5 will slow things down, but it shouldn't be as bad as an
% entire factor of 5.  Plus, I'm running what are usually the fastest
% values of Z(at least for the free livers.. mixing might make things take
% longer with more time scales in the problem).  I may need to use a
% different integrator, i.e. ode15s since we have such disparate time
% scales in this system.
%
% Experimental Desgin (to implement): 1a. Metabolic Classes: Z = 100 for V,
% Z = 10 for I; 1b.two values of Z for parasites: 10^-1 & 10^-2.  need to
% expand null model to include these as well; establish a proper baseline.
%   (p,f) = (I,I); (V,V); (I,V)
% 2. Including fraction of free-living time (Y/N). 3. Concomittant
% predation 4. Fraction of parasites: n=8
%     Separate Run:
% 4. Special scaling on.
%
% This gives 3x2x2x2x2x8 = 192*100 =38400 simulation runs per experiment.
% fairly reasonable; should hopefully be comparable to the BWM experiment.
try
%close all;
saveResults = false;

S = 40;
C = .15;
K = 5;

%The order they appear in correspodns to the order theya re varied.  I
%wrote it differently so easier now to do it like this, then re-arrange the
%coloumns
models = fullfact([2,2,2,2]);

%This changes the order of experiments; species types is now the last thing
%tested so that I can get more intersting factors out of the way first..
models = [models(:,4),models(:,1:3)];

%Models(1): Species Types Models(2): Parasite Z Models(3): Fraction Free
%Living in DE Models(4): concomittant LInks

%For another Experiment.. Level one: metabolic types Future: add an ax for
%parasites
% axFrees = [.314, .88, .88]; axParas = [.314, .88, .314]; spTypesFree = [1
% 2 2]; spTypesPara = [1 2 1]; spTypeNames = {'Inv','EVe'};
%Z for free livers are tied to empirical approximations

%Leel one: BSR exponents for free livers
kFrees = [1,2];

%Level two: BSR exponents for free livers
kParas = [-3,-4];

%Level three: including fraction of free living time
fracFree = [0, 1];

%LEvel four: including concomittant links
Concomittant = [0, 1];

halfSat = 0.5;

Tf = 10000;
nWeb = 100;

dataDir = sprintf(...
    '/Volumes/Buster/Research/Parasitism/Dynamics/JulySimulations/'...
    );

%Fractions that I have run:
%0,.25,.5,.75,1
fParMin = 0;
fParMax = 0.5;
nfPar = 11;
fParAll0 = linspace(fParMin,fParMax,nfPar);
%Future: add a special y value for parasites ys = [8, 4];
h = .2;
count = 0;

SList = 1:S;

axFree = .314;
axPara = .314;

%This tells it to only run the simulations where parasites have a refuge.
%models = models(models(:,3) == 2,:);

tic
%need to do 1-6 now.
for  model = [1;1;1;1]% models'
    %The abundance level web directory
    %     axFree = axFrees(model(1)); axPara = axParas(model(1));
    
    %     typeFreeId = spTypesFree(model(1)); typeParaId =
    %     spTypesPara(model(1));
    %
    %     typeNamesFree = spTypeNames(typeFreeId); typeNamesPara =
    %     spTypeNames(typeParaId);
    
    
    kFree = kFrees(model(1));
    kPara = kParas(model(2));
    
    ZFree =10^kFree;
    
    %     yFree = ys(typeFreeId); yPara = ys(typeParaId);
    yFree = 8;
    yPara = 8;
    
    modelCode = model(3:4);
    
    if sum(modelCode)==2
        fParAll = fParAll0;
        nfPar = numel(fParAll);
    else
        fParAll = fParAll0(fParAll0 ~= 0);
        nfPar = numel(fParAll);
    end
    dataDirS = sprintf('%s/S%03u',dataDir,S);
    
    %The abundance,connectance level web directory
    dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));
    resultsDir = ...
        sprintf('%s/ResultsPATL%u%u%u%u',...
        dataDir,model(1),model(2),model(3),model(4));
    
    codeDir = pwd;
    
    try
        cd(dataDirS)
        cd(codeDir)
    catch
        mkdir(dataDirS)
    end
    
    try
        cd(dataDirSC)
        cd(codeDir)
    catch
        mkdir(dataDirSC)
    end
    
    try
        cd(resultsDir)
        cd(codeDir)
    catch
        mkdir(resultsDir)
    end
    
    for ii = 2
        
        %Try to load the (ii)-th Niche Web:
        try
            LL= csvread(sprintf('%s/web%04u.csv',dataDirSC,ii-1));
            res = LL(:,1);
            con = LL(:,2);
            %otherwise, make one.
        catch MELL
            %Generate the nicheweb:
            webBad = true;
            while webBad
                [res, con,~,~,~] = NicheModel_nk(S,C);
                simMx = calculateSimilarity(res,con);
                mx = sparse(res,con,1,S,S);
                webBad = max(max(simMx))==1;
            end
            csvwrite(sprintf('%s/web%04u.csv',dataDirSC,ii-1),[res,con]);
        end
        %Calculate global and local properties of the web -- this is of
        %trivial complexity compared to the dynamic part.  Could **maybe**
        %save some time by saving these results?  I don't think it's
        %necessarily worth it.
        
        localProperties = calculateLocalProperties(res,con);
        swtl = localProperties(:,10);
        basal = swtl == 1;
        sptb = localProperties(:,8);
        patl = 2*swtl - sptb;
        
        nBasal = sum(basal);
        globalProperties = calculateGlobalProperties(res,con);
        nFree = S-nBasal;
        
        %Assimilation efficiency is given in Brose et al as .45 and .85 for
        %herbivores and carnivores, respectively.  This would correspond to
        %the links we talkin bout.  Easy enough to figure:
        eij = zeros(size(res));
        
        %Future: add an eij for parasitic links
        eij(basal(res)) = .45;
        eij(~basal(res)) = .85;
        K = 5;
        
        %Strong generalist model for preference:
        wij = ones(size(res));
        
        
        try
            nodeData = ...
                csvread(sprintf('%s/IC%04u.csv',dataDirSC,ii-1));
            
            B0 = nodeData(:,1);
            gr = nodeData(:,2);
            
        catch MEmx
            
            B0 = .95*rand(S,1)+.05;
            gr = basal.*(randn(S,1)*.1+1);
            
            csvwrite(sprintf('%s/IC%04u.csv',dataDirSC,ii-1)...
                ,[B0,gr]);
        end
        
        %Order of parasites should be the same.. could save this in the
        %node data?
        try
            fid = fopen(sprintf('%s/idxPar%04u.txt',dataDirSC,ii-1),'r');
            idxPar = textscan(fid,'%u\n','CollectOutput',1);
            idxPar = idxPar{1};
            fclose(fid);
            
        catch MEidxPar
            idxPar = datasample(SList(~basal),nFree,'Replace',false);
            fid = fopen(sprintf('%s/idxPar%04u.txt',dataDirSC,ii-1),'w');
            fprintf(fid,'%u\n',idxPar);
            fclose(fid);
        end
        
        jacobianPattern = full(sparse(res,con,1,S,S));
        jacobianPattern = jacobianPattern + full(sparse(con,res,1,S,S));
        
        jacobianPattern(basal,:) = repmat(basal',sum(basal),1);
        jacobianPattern(:,basal') = repmat(basal,1,sum(basal));
        jacobianPattern(eye(S)>0) = 1;
        opts = odeset('JPattern',sparse(jacobianPattern>0));
        x=zeros(S,1);
        param = struct( 'S',S...
            ,'C',C...
            ,'res',res...
            ,'con',con...
            ,'K',K...
            ,'swtl',swtl...
            ,'eij',eij...
            ,'wij',wij...
            ,'basal',basal...
            ,'r',gr...
            ,'B0',B0...
            ,'h',h...
            ,'x',x...
            ,'halfSat',halfSat...
            ,'Tf',Tf ...
            ,'phi',.15...
            ,'extctThresh',1e-10...
            ,'AbsTol',1e-6...
            ,'RelTol',1e-3...
            ,'modelCode',modelCode...
            ,'odeSolver',@ode45...
            ,'options',opts...
            );
        
        count = 0;
        
        yFinals = zeros(S,1);
        yMeans = zeros(S,1);
        yStds = zeros(S,1);
        
        clc
        curTime = toc;
        
        fprintf('You are currently integrating web %u with %u ',ii,nfPar);
        fprintf('different fractions of parasites.\n');
        fprintf('Parasites have BSR = 10^(%i).\n',kPara);
        fprintf('Free Livers have Z = %.2e.\n',ZFree);
        fprintf('Model code is %u%u.\n',model(3),model(4));
        fprintf('Approximately %.2f minutes have elapsed.\n',toc/60)
        fprintf('fPar = |')
        fprintf(' %.3f |',fParAll)
        fprintf('\n      ')
        TSDir = sprintf('%s/TimeSeries',resultsDir);
        
            
        try
            cd(TSDir);
            cd(codeDir);
        catch
            mkdir(TSDir);
        end
        finalCol = zeros(S,1);
        
        
        finalsFid = fopen(...
            sprintf('%s/finals%04u.csv',resultsDir,ii-1)...
            ,'a');
        
        meansFid = fopen(...
            sprintf('%s/means%04u.csv',resultsDir,ii-1)...
            ,'a');
        
        stdsFid = fopen(...
            sprintf('%s/stds%04u.csv',resultsDir,ii-1)...
            ,'a');
        
        timeSeriesWebDir = sprintf('%s/web%04u',TSDir,ii-1);
        
        try
            cd(timeSeriesWebDir)
            cd(codeDir)
        catch
            mkdir(timeSeriesWebDir)
        end
        
        for fPar = [.20]%fParAll
            timeNow = toc;
            timeSoFar = timeNow - curTime;
            fprintf('\b\b\b\b\b\b-------->%4.0fs',timeSoFar)
            para = false(S,1);
            ax = zeros(S,1);
            y = zeros(S,1);
            M = zeros(S,1);
            
            count = count+1;
            %Choose parasites now. Could also bias selection towards higher
            %trophic levels & low/high mean gen of prey (approximating a 2d
            %RV?)
            
            %I think fraction of free livers as parasites is the best way
            %to go.  Seems better than saying a fixed number of parasites
            %at each level, but I'd have to justify that either way.
            
            nPar = round(fPar*nFree);
            %Actual (realized) fraction of parasites can be found fairly
            %easily as we analyze the webs.  The disconnect from target vs.
            %realized makes analysis at a particular level hard, but maybe
            %makes a regression somewhat more interesting?
            
            para(idxPar(1:nPar)) = true;
            finalCol = finalCol + para;
            free = ~para;
            
            y(free) = yFree;
            y(para) = yPara;
            
            
            ax(free) = axFree;
            ax(para) = axPara;
            
            
            M(free) = ZFree.^(patl(free)-1);
            M(para) = 10.^(kPara + kFree*(patl(para)-2));
            
            x = ax.*M.^(-0.25);
            x(basal) = 0;
            
            
            yij = y(con);
            
            param.x = x;
            param.para = para;
            param.yij = yij;
            param.M = M;
            
            sol = integrateParasiteExperiments(param);
            [order ,sorted ] = sort(sol.extctOrder);
            
            csvwrite(sprintf('%s/fPar%4f.csv',timeSeriesWebDir,fPar),...
                [sol.extctTime(2:end-1);sol.B0s(:,2:end-1)]);
            
            yFinals = sol.y(:,end);
            yMeans = sol.mean;
            yStds = sol.std;
            
            
            if saveResults
                fprintf(finalsFid,'%.2f,%u,',fPar,nPar);
                fprintf(finalsFid,'%.9e,',yFinals);
                fprintf(finalsFid,'\b\n');
                
                fprintf(meansFid,'%.2f,%u,',fPar,nPar);
                fprintf(meansFid,'%.9e,',yMeans);
                fprintf(meansFid,'\b\n');
                
                fprintf(stdsFid,'%.2f,%u,',fPar,nPar);
                fprintf(stdsFid,'%.9e,',yStds);
                fprintf(stdsFid,'\b\n');
                
                
            end
            
        end
        fclose all;
    end
end
finalTime = toc;
fprintf('\nAll Done... Only took %.2f Minutes ;*\n',toc/60);

catch HowItFailed
    WhenItFailed = datestr(clock);
    WhenItFailed(isspace(WhenItFailed)) = '_';
    fileName = sprintf('FuckedUpRun%s%s.mat',mfilename,WhenItFailed);
    fprintf('You Blew it.\n')
    fprintf('At least you had to foresight to save your workspace in:\n')
    fprintf('%s\n',fileName);
    fprintf('Why don''t you repair your data and go fix your darn code?\n')
    save(sprintf('FuckedUpRun%s.mat',WhenItFailed))
end






