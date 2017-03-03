
saveResults = true;
%average of the 6 empirical webs: S = 142, C = .082
S = 50;
C = 0.1;
SList = 1:S;
model = 2;
%The abundance level web directory
dataDirS = sprintf('/Volumes/Buster/Research/Parasitism/Dynamics/NMWebs/S%03u',S);

%The abundance,connectance level web directory
dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));
resultsDir = sprintf('%s/Results%u',dataDirSC,model);

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

%The fraction of parasites will be the primary explanatory variable;
%Parasites will be randomly chosen among the consumers. (realistically, the
%should be preferntially placed at high trophic levels)
nFPar = 10;
fParAll = linspace(0,1,nFPar);
funnyVec = zeros(1,numel(fParAll));
percentChosen = zeros(1,numel(fParAll));
percentChosenO = zeros(1,numel(fParAll));
nWeb = 100;
para = false(S,1);
invert = false(S,1);
vert = false(S,1);
endo = false(S,1);
fPar = 0;
totAlive = zeros(nWeb,length(fParAll));
paraAlive = zeros(nWeb,length(fParAll));
count = 0;

for ii = 1:(nWeb);
    ii
    %Try to load the (nWeb)-th Niche Web:
    try
        LL= csvread(sprintf('%s/web%04u.csv',dataDirSC,ii-1));
        res = LL(:,1);
        con = LL(:,2);
        %otherwise, make one.
        
    catch MELL
        %Generate the nicheweb:
        initialWebBad = true;
        while initialWebBad
            [res, con,~,~,~] = NicheModel_nk(S,C);
            simMx = calculateSimilarity(res,con);
            
            initialWebBad = max(max(simMx))==1;
        end
        csvwrite(sprintf('%s/web%04u.csv',dataDirSC,ii-1),[res,con]);
    end
    %Calculate global and local properties of the web -- this is of trivial
    %complexity compared to the dynamic part.  Could **maybe** save some
    %time by saving these results?  I don't think it's necessarily worth
    %it.
    
    localProperties = calculateLocalProperties(res,con);
    swtl = localProperties(:,10);
    globalProperties = calculateGlobalProperties(res,con);
    %We ned to worry about a few things--
    % x_i - determined from allometric scaling (Williams 2007)
    % f_a - both fractions to take from  boit et al
    % f_m - both """"
    % r_i - follow distribution from Brose 2006
    % K - system wide carrying capacity from Brose 2006
    % e_ij - from Brose 2006
    % w_ij - strong generalist
    %
    
    basal = localProperties(:,2) == 0;
    nBasal = sum(basal);
    
    nFree = S-nBasal;
    %Pick the ectotherm Vertebrates, Invertebrates, and Endotherms
    %Ectotherm Vertebrates and Endotherms are some fraction of the high trophic
    %levels.
    %Q: How do we choose the appropriate a_x&y_ij for the metabolic scaling?
        %A: For the free livers, start with trying each(inconsistent data worries
        %me!)  Can also try a gradient of scaling constants.
        
        %Metabolic scalings from Yodzis & Innes 1992
        %Endo Ver Inv
        axsY = [54.9; 2.3; 0.5];
        
        %Metabolic scalings from Makarieva et al
        %Endo Ver(fish) Inv(aquatic,all); Inv Copepods and Krill; Percarids;
        %decapods
        % 38.9;.646;.813;.437;.617;.891
        %Ectotherms all about the same scale. O(.1~1) vs O(10)
        axsM = [10^(1.59); 10^(-0.19); 10^(-0.09); 10^(-.36); 10^(-.21); 10^(-.05)];
        
        %The variation in the above suggests that since there is variability
        %according to subdivisions, we could concievably have many different
        %constants.  Maybe the precies value doesn't matter and we can simply have
        %a distribution of these as well?
        
        %Metabolic scalings from Brose 2006
        %Endo (not given); Vert; Invert
        axsB = [.88; .314];
        
        %For the first set of simulations, take a_x = 0.5.
        
        
        
        axFree=0.314;
        %.314 for model 0.  .0001 otherwise
        axPara = 0.0001; %Coralie: .0001
        
         %Maximum ingestion rates for ectotherm invertebrates and vertebrates are
        %given in Brose 2006 as 8 and 4, respectively (when these are the only
        %species types in the web).  Doesn't give info for predation on plants.
        %What about endotherms?
        %Yodzis and Innes gives this:
        %(res,con); i = endo,vert,invert
        aJ = [89.2; 8.9; 9.7];
        aT = [54.9;2.3;0.5];
        fJ = [1;.2;.3];
        Y = fJ*(aJ./aT)';
        
        %use Y for more complicated arrangements;
        y = 8;
        yij = y * ones(size(res));
        %Assimilation efficiency is given in Brose et al as .45 and .85 for
        %herbivores and carnivores, respectively.  This would correspond to the
        %links we talkin bout.  Easy enough to figure:
        eij = zeros(size(res));
        
        eij(basal(res)) = .66;
        eij(~basal(res)) = .85;
        K = 450;
        
        %Using strong generalist model for preference.  (compare to weak
        %generalist?)
        wij = ones(size(res));
        
        h=1.2;
        
        halfSat = 80;
        %rand(S,1)*390+10;
    try
        nodeData = ...
            csvread(sprintf('%s/nodeData%04u.csv',dataDirSC,ii-1));
        ZFree = nodeData(:,1);
        ZPara = nodeData(:,2);
        mFree = nodeData(:,3);
        mPara = nodeData(:,4);
        r = nodeData(:,7);
        B0 = nodeData(:,8);
        
        %For testing MOdel 1
        xFree = axFree.*mFree.^(-.25);
        xPara = axPara.*mPara.^(-.25);
        
    catch MEmx
        %First approximation is free livers vs parasites.  This is average and
        %standard deviation of ln(con/res body mass).  These numbers are from the 3
        %webs published by Hechinger etal 2011.
        %The parasite number includes mosquitos since they represent that inverse
        %relationship we are trying to test - thus muParaMosq is the mean
        %of parasites, including mosquitos, same for stdmosq.
        
        %Saved as meanFree,MeanPA_Mosq,stdFree,&stdPA_Mosq if you run BigStats
        %first.
        %
        muFree = 6.89;
        stdFree = 4.09;
        
        muParaMosq = -12.59;
        %muParaMosq = -5.59; %Coralie = -12.89,3.90
        stdParaMosq = 3.98;
        
        ZFree = lognrnd(muFree,stdFree,S,1);
        %ZFree = 1000;
        ZPara = lognrnd(muParaMosq,stdParaMosq,S,1);
        %ZPara = .1;
        %These rules could be used later:
        %
        %Of the remaining species...
        %1. Pick 90% of the species below TL 3 to be invertebrates
        %2. Pick 25% of the species above TL 3 to be Ectotherm vertebrates
        %3. Pick 45% of the species above TL 3 to be Endotherms
        %4. Remainder of TL>3 are ectotherm invertebrates.
        %5. Remaidner of TL<3 are Endotherms.
        %{
%breaking apart the free liver metabolic classes
spTypeDraw = rand(S,1);

ectoInv1 = (swtl<3)&(~basal)&(spTypeDraw<0.9);
endo1 = (swtl<3)&(~basal)&(spTypeDraw>=0.9);
endo2 = (swtl>3)&(~basal)&(spTypeDraw>0.55);
ectoVert = (swtl>3)&(~basal)&(spTypeDraw<0.25);
ectoInv2 = (swtl>3)&(~basal)&(~endo2)&(~ectoVert);

ectoInv = ectoInv1|ectoInv2;
endo = endo1|endo2;

muEctoInv = 5.60;
stdEctoInv = 5.72;

muEctoVert = 6.51;
stdEctoVert = 4.01;

muEndo = 7.63;
stdEndo = 3.21;

muPara = -12.4;
stdPara = 5.02;

ZEI = lognrnd(muEctoInv,stdEctoInv,S,1);
ZEI = ZEI.*ectoInv;

ZEV = lognrnd(muEctoVert,stdEctoVert,S,1);
ZEV = ZEV.*ectoVert;

ZEN = lognrnd(muEndo,stdEndo,S,1);
ZEN = ZEN.*endo;
ZFree = ZEI + ZEV + ZEN;
axFree = zeros(S,1);
    axFree(ectoInv) = 10^(-0.19);
    axFree(endo) = 10^(-0.2);
    axFree(ectoVert) = 10^(-0.36);
        %}
        
        
        
        
        mFree = ZFree.^(swtl-1);
        mFree(ZFree==0) = 1;
        mPara = ZPara.^(swtl-1);
        
        
        
        xFree = axFree.*mFree.^(-.25);
        xPara = axPara.*mPara.^(-.25);
        
        
        
        %x=xFree;
        
       
        
        mu_r = .9;
        sd_r = 0.2;
        r = randn(S,1)*sd_r + mu_r;
        r(r<(mu_r-3*sd_r)) = mu_r-3*sd_r;
        r(r>(mu_r+3*sd_r)) = mu_r+3*sd_r;
        r(~basal) = 0;
        %r(basal) = 1;
        
        B0 = ones(S,1)*10;
        
        csvwrite(sprintf('%s/nodeData%04u.csv',dataDirSC,ii-1)...
            ,[ZFree,ZPara,mFree,mPara,xFree,xPara,r,B0]);
    end
    
    %Order of parasites should be the same.. could save this in the node
    %data?
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
    
    
    
    %Exclude in Models 0 and 1.  Include otherwise:
    %phi = .15*para;
    x=zeros(S,1);
    param = struct( 'S',S...
        ,'C',C...
        ,'res',res...
        ,'con',con...
        ,'K',K...
        ,'para',[]...
        ,'swtl',swtl...
        ,'eij',[]...
        ,'wij',wij...
        ,'yij',yij...
        ,'basal',basal...
        ,'r',r...
        ,'B0',B0...
        ,'h',h...
        ,'x',[]...
        ,'halfSat',halfSat...
        ,'Tf',8000 ...
        ,'phi',phi ...exclude in model 0
        );
            
    count = 0;
    
    yFinals = zeros(S,nFPar);
    yMeans = zeros(S,nFPar);
    yStds = zeros(S,nFPar);
    for fPar = fParAll
        
        para = false(S,1);
        count = count+1;
        %Choose parasites now.
        %Could also bias selection towards higher trophic levels & low/high mean
        %gen of prey (approximating a 2d RV?)
        nPar = round(fPar*nFree);
        
        
        
        para(idxPar(1:nPar)) = true;
        
        %potentialConcomittantLinks = ...
        %    identifyPotentialConcomittant(res,con,para);

        
        %previewing intelligent selection of links
%         funnyVec(ii,count) = length(potentialConcomittantLinks)/length(res);
%         m = zeros(S,1);
%         m = mFree.*(~para)+mPara.*(para);
%         data = [swtl(potentialConcomittantLinks) log(m(potentialConcomittantLinks))];
%         
%         [~,dataO1] = sort(data(:,1));
%         [~,dataO2] = sort(data(:,2));
%         [~,dataO3] = sort(data(:,3));
%         [~,dataO4] = sort(data(:,4));
%         
%         dataO = [dataO1 dataO2 dataO3 dataO4]/length(dataO1);
%         
%         concInd = predict(tree,data);
%         concIndO = predict(treeO,dataO);
%         percentChosen(ii,count) = sum(concInd)/length(potentialConcomittantLinks);
%         percentChosenO(ii,count) = sum(concIndO)/length(potentialConcomittantLinks);
        
        %Exclude in Model 0
        eij(para(con)) = .90;
        
        x(para) = xPara(para);
        x(~para) = xFree(~para);
        x(basal) = 0;
        
        param.x = x;
        param.para = para;
        param.eij = eij;
        
        sol = integrateATN(res,con,param);
        
        totAlive(ii,count) = sum(sol.(sprintf('part%u',sol.n)).y(:,end)>0);
        paraAlive(ii,count) = sum(sol.(sprintf('part%u',sol.n)).y(para,end)>0)/nPar;
        yFinals(:,count) = sol.yT;
        yMeans(:,count) = sol.mean;
        yStds(:,count) = sol.std;
    end
    if saveResults
        csvwrite(sprintf('%s/finals%04u.csv',resultsDir,ii-1),yFinals);
        csvwrite(sprintf('%s/means%04u.csv',resultsDir,ii-1),yMeans);
        csvwrite(sprintf('%s/stds%04u.csv',resultsDir,ii-1),yStds);
    end
end

% figure
% plot(sol.(sprintf('part%u',sol.n)).x,sol.(sprintf('part%u',sol.n)).y)
% [sol.(sprintf('part%u',sol.n)).y(:,end) basal r x]

