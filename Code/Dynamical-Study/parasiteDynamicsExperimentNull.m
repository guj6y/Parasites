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
close all;
saveResults = true;

S = 40;
C = .15;
K = 5;

%The order they appear in correspodns to the order theya re varied.  I
%wrote it differently so easier now to do it like this, then re-arrange the
%coloumns
models = [1 1 0 0;
          2 1 0 0];


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

%Want to include no parasites as the null model?  No -- have those data
%already, also, want to re-run with the lower ratios.
fParAll = 0;
nfPar = 1;

%Future: add a special y value for parasites ys = [8, 4];
h = 1.2;
count = 0;

SList = 1:S;

axFree = .314;
axPara = .314;
tic
for  model = models'
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
    
    dataDirS = sprintf('%s/S%03u',dataDir,S);
    
    %The abundance,connectance level web directory
    dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));
    resultsDir = ...
        sprintf('%s/ResultsPATL%u%u%u%u',...
        dataDir,model(1),model(2),model(3),model(4));
    model(2:4) = 1;
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
    
    for ii = 1:(nWeb);
        
        %Try to load the (nWeb)-th Niche Web:
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
            B0 = ...
                csvread(sprintf('%s/IC%04u.csv',dataDirSC,ii-1));
            
            
        catch MEmx
            
            B0 = .95*rand(S,1)+.05;
            
            csvwrite(sprintf('%s/IC%04u.csv',dataDirSC,ii-1)...
                ,B0);
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
        
        
        r = basal.*(randn(S,1)*.1+1);
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
            ,'r',r...
            ,'B0',B0...
            ,'h',h...
            ,'x',x...
            ,'halfSat',halfSat...
            ,'Tf',Tf ...
            ,'phi',.15...
            ,'extctThresh',1e-10...
            ,'AbsTol',1e-10...
            ,'RelTol',1e-7...
            ,'modelCode',modelCode...
            );
        
        count = 0;
        
        yFinals = zeros(S,nfPar+1);
        yMeans = zeros(S,nfPar+1);
        yStds = zeros(S,nfPar+1);
        
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
        fprintf('\n ')
        TSDir = sprintf('%s/TimeSeries',resultsDir);
        
        try
            cd(TSDir);
            cd(codeDir);
        catch
            mkdir(TSDir);
        end
        finalCol = zeros(S,1);
        
        try
            csvread(sprintf('%s/finals%04u.csv',resultsDir,ii-1));
        catch
            
            for fPar = fParAll
                
                fprintf('\b-------->')
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
                
                
                
                yFinals(:,count) = sol.y(:,end);
                yMeans(:,count) = sol.mean;
                yStds(:,count) = sol.std;
                save(sprintf('%s/sol%04ufPar%u.mat',TSDir,ii-1,count),...
                    '-struct','sol');
                %Get this time series by sol = load(sprintf(%s/sol%...);
            end
            if saveResults
                yFinals(:,end) = finalCol;
                yMeans(:,end) = finalCol;
                yStds(:,end) = finalCol;
                csvwrite(sprintf('%s/finals%04u.csv',resultsDir,ii-1),yFinals);
                csvwrite(sprintf('%s/means%04u.csv',resultsDir,ii-1),yMeans);
                csvwrite(sprintf('%s/stds%04u.csv',resultsDir,ii-1),yStds);
                
            end
        end
    end
end
finalTime = toc;
fprintf('\nAll Done... Only took %.2f Minutes ;*\n',toc/60);








