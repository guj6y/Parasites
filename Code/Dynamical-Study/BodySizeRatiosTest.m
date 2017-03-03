close all
webGeneration;

%Looking at effect of changing Z on mean body size ratios, maximum body
%sizes, and minimum metabolic rates.  Body size an metabolic generally have
%some predictive power for time of extinction.  Trying to see if using PATL
%has similar effect as higher Z for SWTL.  This could be important for the
%results of the model I'm working on now (changing portions of species with
%parasites vs. adding parasites to web nuggets).

nWebs = 1000;
S = 40;
C = 0.15;

fParMin = 0.05;
fParMax = .75;
nfPar = 8;
fParAll = linspace(fParMin,fParMax,nfPar);
fParAll = 0.3;

sList = 1:S;

dataDir = '~/Desktop/JulySimulations';
dataDirS = sprintf('%s/S%03u',dataDir,S);
dataDirSC = sprintf('%s/C%03u',dataDirS,C*100);
codeDir = pwd;

%%%%%%%%%%%% This controls the fraction of parasites - Important!  Huge!
fPar = 0.3;
%No parasites for now..

try
    cd(dataDirS)
    cd(codeDir)
catch
    mkdir(dataDirS);
end

try
    cd(dataDirSC)
    cd(codeDir)
catch
    mkdir(dataDirSC);
end
nZ = 8;
kFrees = 2;
kFrees = linspace(-2,5,nZ)';
kPara = -5;

kCount = 0;
bsEmpirical = zeros(0,3);
avgBSResEmp = [];
avgBSConEmpPara = [];
avgBSConEmpFree = [];
bsAllEmp = [];
paraAllEmp = [];
bsBasal = 1e-7;
nAll = [];

for jj = 1:6
    LL = linkListCell{jj};
    Sjj = max(max(LL));
    bsjj = propertiesCell{jj,1}(:,12);
    %bsjj(binAutoCell{jj}) = bsBasal;
    %bsjj = bsjj/bsBasal;
    %bsjj = bsjj/min(bsjj);
    resjj = LL(:,1);
    conjj = LL(:,2);
    
    parajj = binParaCell{jj};
    paraConjj = parajj(conjj);
    paraResjj = parajj(resjj);
    
    bsEmpirical = [bsEmpirical;...
        [bsjj(resjj), bsjj(conjj), 5*ones(size(resjj))+paraConjj+2*paraResjj]];
    
    bsResMx = full(sparse(resjj,conjj,bsjj(resjj),Sjj,Sjj));
    bsResMx(bsResMx==0) = nan;
    
    avgBSRes= mean(log10(bsResMx),'omitnan');
    
    bsConMx = full(sparse(resjj,conjj,bsjj(conjj),Sjj,Sjj));
    
    bsConMxFree = bsxfun(@times,bsConMx,~parajj);
    bsConMxFree(bsConMxFree==0) = nan;
    avgBSConFree = mean(log10(bsConMxFree),2,'omitnan');
    
    bsConMxPara = bsxfun(@times,bsConMx,parajj);
    bsConMxPara(bsConMxPara==0) = nan;
    avgBSConPara = mean(log10(bsConMxPara),2,'omitnan');
    
    avgBSResEmp = [avgBSResEmp;avgBSRes'];
    avgBSConEmpPara = [avgBSConEmpPara;avgBSConPara];
    avgBSConEmpFree = [avgBSConEmpFree;avgBSConFree];
    bsAllEmp = [bsAllEmp;bsjj];
    paraAllEmp = [paraAllEmp;parajj];
end
meanBSSWTL = zeros(0,8);
meanBSPATL = zeros(0,8);

bsrSWTLAll = zeros(0,8);
bsrPATLAll = zeros(0,8);

meanGenParaAll = zeros(0,16);
genAll = [];
vulAll = [];
for kFree = kFrees
    countfPar = 0;
    meanGenParaAll = zeros(0,16);
    for fPar = fParAll
        countfPar =  countfPar + 1;
        swtlAll = [];
        patlAll = [];
        MAll = zeros(0,2);
        paraAll = [];
        bsSWTL = zeros(0,3);
        bsPATL = zeros(0,3);
        
        avgBSResSimSWTL = [];
        avgBSConSimSWTLpara = [];
        avgBSConSimSWTLfree = [];
        
        avgBSResSimPATL = [];
        avgBSConSimPATLpara = [];
        avgBSConSimPATLfree = [];
        MAll = zeros(0,2*nZ);
        for ii = 1:nWebs
            
            if mod(ii,100) ==0
                clc
                fprintf('%.2f Percent Complete.\n',(ii+(kCount)*1000)/10000)
            end
            kCount = kCount+1;
            try
                IC0 = csvread(sprintf('%s/IC%04u.csv',dataDirSC,ii-1));
                
            catch
                IC0 = rand(S,1)*.95 + .05;
                csvwrite(sprintf('%s/IC%04u.csv',dataDirSC,ii-1),IC0);
            end
            
            try
                ncr = csvread(sprintf('%s/ncr%04u.csv',dataDirSC,ii-1));
                LL = csvread(sprintf('%s/web%04u.csv',dataDirSC,ii-1));
                res = LL(:,1);
                con = LL(:,2);
                n = ncr(:,1);
                c = ncr(:,2);
                r = ncr(:,3);
                
            catch
                webBad= true;
                
                while webBad
                    [res, con,n,c,r] = NicheModel_nk(S,C);
                    simMx = calculateSimilarity(res,con);
                    mx = sparse(res,con,1,S,S);
                    webBad = max(max(simMx))==1;
                end
                csvwrite(sprintf('%s/web%04u.csv',dataDirSC,ii-1),[res,con]);
                csvwrite(sprintf('%s/ncr%04u.csv',dataDirSC,ii-1),[n,c,r]);
            end
            
            nAll = [nAll;n];
            
            localProps = calculateLocalProperties(res,con);
            swtl = localProps(:,10);
            stl = localProps(:,8) + 1;
            gen = localProps(:,2);
            vul = localProps(:,3);
            vulAbs = vul*length(res)/S;
            genAbs = gen*length(res)/S;
            fpar = 0;
            basal = swtl==1;
            nFree = S-sum(basal);
            patl = 2*swtl - stl;

            try
                idxPar = csvread(sprintf('%s/idxPar%04u.csv',dataDirSC,ii-1));
                
            catch
                idxPar = datasample(sList(~basal),nFree,'Replace',false);
                csvwrite(sprintf('%s/idxPar%04u.csv',dataDirSC,ii-1),idxPar);
            end
            
            
            count = 0;
            
            %need to choose parasites now
            para = false(S,1);
            MSWTL = zeros(S,nZ);
            MPATL = zeros(S,nZ);
            
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
            para = para>0;
            free = ~para;
            
            ZFree = bsxfun(@power,(10*ones(S,1)),kFree');
            %ZFree = lognrnd(4.6,2.3,S,1);
            
            MSWTL(free,:) = bsxfun(@power,ZFree(free,:),(swtl(free)-1));
            MPATL(free,:) = bsxfun(@power,ZFree(free,:),(patl(free)-1));
            ZPara = (10*ones(S,1)).^kPara;
            %ZPara = lognrnd(-9.2,2.3,S,1);
            MSWTL(para>0,:) = bsxfun(@power,ZPara(para),((kPara+bsxfun(@times,kFree',(swtl(para>0)-2)))/kPara));
            MPATL(para>0,:) = bsxfun(@power,ZPara(para),((kPara+bsxfun(@times,kFree',(patl(para>0)-2)))/kPara));
            MSWTL(basal>0,:) = nan;
            MPATL(basal>0,:) = nan;
            MAll = [MAll;[MSWTL MPATL]];
            
            
            
             paraAll = [paraAll;para];
             
             bsrSWTL = MSWTL(con,:)./MSWTL(res,:);
             bsrPATL = MPATL(con,:)./MPATL(res,:);
             
             paraConLink = para(con);
             
             paraConParaLink = bsxfun(@and,bsrSWTL<1,paraConLink);
             paraConTropLink = bsxfun(@and,bsrSWTL>1,paraConLink);
             

             meanGenParaPara = sum(paraConParaLink)/length(con)*S/nPar;
             meanGenParaTrop = sum(paraConTropLink)/length(con)*S/nPar;
             
             meanGenParaAll = [meanGenParaAll; [meanGenParaPara meanGenParaTrop]];
             
%             %Bad way to do Parasites:
%             %M(para) = 10.^(kPara*(swtl(para)-1));
%              bsSWTL = [bsSWTL;[MSWTL(res,:),MSWTL(con,:)]];
%              bsPATL = [bsPATL;[MPATL(res,:),MPATL(con,:)]];
%             
%             meanBSSWTL = [meanBSSWTL;mean(log10(MSWTL(con,:)./MSWTL(res,:)))];
%             meanBSPATL = [meanBSPATL;mean(log10(MPATL(con,:)./MPATL(res,:)))];
%             
             swtlAll = [swtlAll;swtl];
             patlAll = [patlAll;patl];
            vulAll = [vulAll;vulAbs];
            genAll = [genAll;genAbs];
            
            
%             bsResMxSWTL = full(sparse(res,con,MSWTL(res),S,S));
%             bsResMxSWTL(bsResMxSWTL==0) = nan;
%             
%             avgBSResSWTL = mean(log10(bsResMxSWTL),'omitnan');
%             
%             bsConMxSWTL = full(sparse(res,con,MSWTL(con),S,S));
%             
%             bsConMxSWTLpara = bsxfun(@times,bsConMxSWTL,basal);
%             bsConMxSWTLpara(bsConMxSWTLpara==0) = nan;
%             avgBSConSWTLpara = mean(log10(bsConMxSWTLpara),2,'omitnan');
%             
%             bsConMxSWTLfree = bsxfun(@times,bsConMxSWTLpara,~basal);
%             bsConMxSWTLfree(bsConMxSWTLfree==0) = nan;
%             avgBSConSWTLfree = mean(log10(bsConMxSWTLfree),2,'omitnan');
%             
%             avgBSResSimSWTL = [avgBSResSimSWTL;avgBSResSWTL'];
%             avgBSConSimSWTLpara = [avgBSConSimSWTLpara;avgBSConSWTLpara];
%             avgBSConSimSWTLfree = [avgBSConSimSWTLfree;avgBSConSWTLfree];
%             
%             bsResMxPATL = full(sparse(res,con,MPATL(res),S,S));
%             bsResMxPATL(bsResMxPATL==0) = nan;
%             
%             avgBSResPATL = mean(log10(bsResMxPATL),'omitnan');
%             
%             bsConMxPATL = full(sparse(res,con,MPATL(con),S,S));
%             
%             bsConMxPATLpara = bsxfun(@times,bsConMxPATL,basal);
%             bsConMxPATLpara(bsConMxPATLpara==0) = nan;
%             avgBSConPATLpara = mean(log10(bsConMxPATLpara),2,'omitnan');
%             
%             bsConMxPATLfree = bsxfun(@times,bsConMxPATLpara,~basal);
%             bsConMxPATLfree(bsConMxPATLfree==0) = nan;
%             avgBSConPATLfree = mean(log10(bsConMxPATLfree),2,'omitnan');
%             
%             avgBSResSimPATL = [avgBSResSimPATL;avgBSResPATL'];
%             avgBSConSimPATLpara = [avgBSConSimPATLpara;avgBSConPATLpara];
%             avgBSConSimPATLfree = [avgBSConSimPATLfree;avgBSConPATLfree];
        end
        
        
%         minZ = 3;
%         maxZ = 5;
%         figure;
%         oneMx = ones(size(meanGenParaAll(:,minZ:maxZ)));
%         groupsMx = bsxfun(@times,oneMx,minZ:maxZ);
%         groupsVec = reshape(groupsMx,[],1);
%         meanGenParaAllTro = reshape(meanGenParaAll(:,minZ:maxZ),[],1);
%         meanGenParaAllPar = reshape(meanGenParaAll(:,8+(minZ:maxZ)),[],1);
%         gscatter(meanGenParaAllPar,meanGenParaAllTro,groupsVec,[],'.',15);
%         title(sprintf('fPar = %.2f',fPar))
%         xlabel('Mean Gen as parasite')
%         ylabel('Mean gen. as trophic')

%         bsPlot = [bsSWTL; bsEmpirical];
%         h0 = figure('Position', [0, 0, 1440, 900]);
%         h = gscatter(log10(bsPlot(:,1)),log10(bsPlot(:,2)),bsPlot(:,3));
%         set(h,'MarkerSize',10)
%         refline(1,kPara)
%         refline(1,kFree)
%         refline(1,2*kFree-kPara);
%         xlabel('log_{10}(BS_{res})')
%         ylabel('log_{10}(BS_{con})')
%         axis([min(log10(bs(:,1))) max(log10(bs(:,1))) min(log10(bs(:,2))) max(log10(bs(:,2)))]);
%         saveas(h0,'ZStudyScatter.jpg','jpg')
    end
end
% figure;
%         hist(swtlAll(swtlAll>1),100);
%         title('short-weighted')
%         figure;
%         hist(patlAll(patlAll>1),100)
%         title('prey-averaged')
omn = mod(swtlAll,1)==0;
omn2 = omn;
bas = swtlAll==1;

% figure;
%  hist(swtlAll(~(omn|bas)),100);
%         title('short-weighted')
xAll = MAll.^-.25;
plot(genAll,MAll)
        figure;
        binEdges = linspace(2,15,131);
        binCenters = linspace(2.05,14.95,130);
        
        [N,Edges] = histcounts(patlAll(~(omn2|bas)),binEdges,...
            'Normalization','pdf');
        bar(binCenters,N,1)
        [PHAT,PCI] = mle(patlAll(~(omn2|bas)),'distribution','gamma');
        pdfEst = gampdf(binCenters,PHAT(1),PHAT(2));
        line(binCenters,pdfEst,'LineWidth',2);
        title('prey-averaged')

        figure;
        hold on
        ax = histogram(patlAll(omn)-1,-0.5:16.5,'normalization','pdf');
       [PhatGeo,PCI2] = mle(patlAll(omn)-1,'distribution','geo');
       [PhatPoiss,PCI3] = mle(patlAll(omn)-1,'distribution','poiss');
       pdfEstGeo = geopdf(0:15,PhatGeo(1));
       pdfEstpoiss = poisspdf(0:15,Phat(1));
        errorbar(0:15,pdfEstGeo,'sk','MarkerSize',10);
        errorbar(0:15,pdfEstpoiss,'sr','MarkerSize',10);
        title('prey-averaged')
       
       
%Making and Saving a **TON** of histograms (okay like 10 but still it's
%ugly)
%{
bs = bsSWTL;
%%%%% Figure 1
h1 = figure('Position', [0, 0, 1440, 900]);

subplot(2,3,1)
histRatio = log10(bsEmpirical(bsEmpirical(:,3)==5|bsEmpirical(:,3)==7,2)./...
    bsEmpirical(bsEmpirical(:,3)==5|bsEmpirical(:,3)==7,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
ylabel('Count {SWTL}')
title('Emp. Distribution of consumer-resources bsr, non-parasite');
XL1 = xlim;

subplot(2,3,2)
histRatio = log10(bsEmpirical(bsEmpirical(:,3)==6|bsEmpirical(:,3)==8,2)./...
    bsEmpirical(bsEmpirical(:,3)==6|bsEmpirical(:,3)==8,1));
hist(histRatio ,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. Distribution of all consumer-resources bsr, parasite');
XL2 = xlim;

subplot(2,3,3)
histRatio = log10(bsEmpirical(:,2)./bsEmpirical(:,1));
hist(histRatio,20);
XL3 = xlim;
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25)))
title('Emp. Distribution of all consumer-resources bsr');


subplot(2,3,4)
histRatio = log10(bs(bs(:,3)==0|bs(:,3)==2,2)./bs(bs(:,3)==0|bs(:,3)==2,1));
hist(histRatio,100)
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title(sprintf('Sim. non-parasite consmer-resource bsr for 1000 S=%u,C=0.15 NM',S))
%xlim(XL1);

subplot(2,3,5)
histRatio = log10(bs(bs(:,3)==1|bs(:,3)==3,2)./bs(bs(:,3)==1|bs(:,3)==3,1));
hist(histRatio,100)
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title(sprintf('Sim. non-parasite consmer-resource bsr for 1000 S=%u,C=0.15 NM',S))
%xlim(XL2);

subplot(2,3,6);
histRatio = log10(bs(:,2)./bs(:,1));
hist(histRatio,100)
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
%xlim(XL3);
title(sprintf('Sim. non-parasite consmer-resource bsr for 1000 S=%u,C=0.15 NM',S))

subplot(2,3,1);
%xlim(XL3);
saveas(h1,'ZStudyFig1SWTL.jpg','jpg')

%%%%% Figure 2
h2 = figure('Position', [0, 0, 1440, 900]);
subplot(2,2,1)
histRatio= log10(bsEmpirical(bsEmpirical(:,3)==6,2)./...
    bsEmpirical(bsEmpirical(:,3)==6,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
ylabel('Count {SWTL}')
title('Emp. Distribution of all parasite-host body size ratio');
XL1 = xlim;

subplot(2,2,2)
histRatio = log10(bsEmpirical(bsEmpirical(:,3)==8,2)./...
    bsEmpirical(bsEmpirical(:,3)==8,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. Distribution of all parasite-parasite body size ratio');
XL2 = xlim;

subplot(2,2,3)
histRatio = log10(bs(bs(:,3)==1,2)./bs(bs(:,3)==1,1));
hist(histRatio,100);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Distribution of all parasite-host body size ratio');
%xlim(XL1);

subplot(2,2,4)
histRatio = log10(bs(bs(:,3)==3,2)./bs(bs(:,3)==3,1));
hist(histRatio,100);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Distribution of all parasite-parasite body size ratio');
%xlim(XL2);
saveas(h2,'ZStudyFig2SWTL.jpg','jpg')

%%%%% Figure 3
h3 = figure('Position', [0, 0, 1440, 900]);
subplot(2,2,1)
histRatio = log10(bsEmpirical(bsEmpirical(:,3)==5,2)./...
    bsEmpirical(bsEmpirical(:,3)==5,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
ylabel('Count {SWTL}')
title('Emp. Distribution of all free<-free body size ratio');
XL1 = xlim;

subplot(2,2,2)
histRatio = log10(bsEmpirical(bsEmpirical(:,3)==7,2)./...
    bsEmpirical(bsEmpirical(:,3)==7,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. Distribution of all free<-parasite body size ratio');
XL2 = xlim;

subplot(2,2,3)
histRatio = log10(bs(bs(:,3)==0,2)./bs(bs(:,3)==0,1));
hist(histRatio,100);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Distribution of all free<-free body size ratio');
%xlim(XL1);

subplot(2,2,4)
histRatio = log10(bs(bs(:,3)==2,2)./bs(bs(:,3)==2,1));
hist(histRatio,100);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Distribution of all free<-parasite body size ratio');
%xlim(XL2);
saveas(h3,'ZStudyFig3SWTL.jpg','jpg')

%%%%% Figure 4
h4 = figure('Position', [0, 0, 1440, 900]);
subplot(2,3,1)
histRatio = log10(bsAllEmp(~paraAllEmp))-avgBSResEmp(~paraAllEmp);
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
ylabel('Count {SWTL}')
title('Emp. Free consumer:mean Resource BSR');
XL1 = xlim;

subplot(2,3,2)
histRatio = log10(bsAllEmp(paraAllEmp>0))- avgBSResEmp(paraAllEmp>0);
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. Parasite consumer:mean Resource BSR');
XL2 = xlim;

subplot(2,3,3)
histRatio = log10(bsAllEmp)- avgBSResEmp;
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. All Consumer:Resource BSR');
XL3 = xlim;

subplot(2,3,4)
 histRatio = log10(MAll(~paraAll,1))-avgBSResSim(~paraAll);
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Free consumer:mean Resource BSR');
%xlim(XL1);

subplot(2,3,5)
 histRatio = log10(MAll(paraAll>0,1))-avgBSResSim(paraAll>0);
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Para consumer:mean Resource BSR');
%xlim(XL2);

subplot(2,3,6)
 histRatio = log10(MAll(:,1))-avgBSResSim;
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. All consumer:mean Resource BSR');
%xlim(XL3);
saveas(h4,'ZStudyFig4SWTL.jpg','jpg')


bs = bsPATL;
%%%%% Figure 1
h5 = figure('Position', [0, 0, 1440, 900]);

subplot(2,3,1)
histRatio = log10(bsEmpirical(bsEmpirical(:,3)==5|bsEmpirical(:,3)==7,2)./...
    bsEmpirical(bsEmpirical(:,3)==5|bsEmpirical(:,3)==7,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
ylabel('Count {PATL}')
title('Emp. Distribution of consumer-resources bsr, non-parasite');
XL1 = xlim;

subplot(2,3,2)
histRatio = log10(bsEmpirical(bsEmpirical(:,3)==6|bsEmpirical(:,3)==8,2)./...
    bsEmpirical(bsEmpirical(:,3)==6|bsEmpirical(:,3)==8,1));
hist(histRatio ,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. Distribution of all consumer-resources bsr, parasite');
XL2 = xlim;

subplot(2,3,3)
histRatio = log10(bsEmpirical(:,2)./bsEmpirical(:,1));
hist(histRatio,20);
XL3 = xlim;
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25)))
title('Emp. Distribution of all consumer-resources bsr');


subplot(2,3,4)
histRatio = log10(bs(bs(:,3)==0|bs(:,3)==2,2)./bs(bs(:,3)==0|bs(:,3)==2,1));
hist(histRatio,100)
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title(sprintf('Sim. non-parasite consmer-resource bsr for 1000 S=%u,C=0.15 NM',S))
%xlim(XL1);

subplot(2,3,5)
histRatio = log10(bs(bs(:,3)==1|bs(:,3)==3,2)./bs(bs(:,3)==1|bs(:,3)==3,1));
hist(histRatio,100)
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title(sprintf('Sim. non-parasite consmer-resource bsr for 1000 S=%u,C=0.15 NM',S))
%xlim(XL2);

subplot(2,3,6);
histRatio = log10(bs(:,2)./bs(:,1));
hist(histRatio,100)
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
%xlim(XL3);
title(sprintf('Sim. non-parasite consmer-resource bsr for 1000 S=%u,C=0.15 NM',S))

subplot(2,3,1);
%xlim(XL3);
saveas(h5,'ZStudyFig1PATL.jpg','jpg')

%%%%% Figure 2
h6 = figure('Position', [0, 0, 1440, 900]);
subplot(2,2,1)
histRatio= log10(bsEmpirical(bsEmpirical(:,3)==6,2)./...
    bsEmpirical(bsEmpirical(:,3)==6,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
ylabel('Count {PATL}')
title('Emp. Distribution of all parasite-host body size ratio');
XL1 = xlim;

subplot(2,2,2)
histRatio = log10(bsEmpirical(bsEmpirical(:,3)==8,2)./...
    bsEmpirical(bsEmpirical(:,3)==8,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. Distribution of all parasite-parasite body size ratio');
XL2 = xlim;

subplot(2,2,3)
histRatio = log10(bs(bs(:,3)==1,2)./bs(bs(:,3)==1,1));
hist(histRatio,100);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Distribution of all parasite-host body size ratio');
%xlim(XL1);

subplot(2,2,4)
histRatio = log10(bs(bs(:,3)==3,2)./bs(bs(:,3)==3,1));
hist(histRatio,100);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Distribution of all parasite-parasite body size ratio');
%xlim(XL2);
saveas(h6,'ZStudyFig2PATL.jpg','jpg')

%%%%% Figure 3
h7 = figure('Position', [0, 0, 1440, 900]);
subplot(2,2,1)
histRatio = log10(bsEmpirical(bsEmpirical(:,3)==5,2)./...
    bsEmpirical(bsEmpirical(:,3)==5,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
ylabel('Count {PATL}')
title('Emp. Distribution of all free<-free body size ratio');
XL1 = xlim;

subplot(2,2,2)
histRatio = log10(bsEmpirical(bsEmpirical(:,3)==7,2)./...
    bsEmpirical(bsEmpirical(:,3)==7,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. Distribution of all free<-parasite body size ratio');
XL2 = xlim;

subplot(2,2,3)
histRatio = log10(bs(bs(:,3)==0,2)./bs(bs(:,3)==0,1));
hist(histRatio,100);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Distribution of all free<-free body size ratio');
%xlim(XL1);

subplot(2,2,4)
histRatio = log10(bs(bs(:,3)==2,2)./bs(bs(:,3)==2,1));
hist(histRatio,100);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Distribution of all free<-parasite body size ratio');
%xlim(XL2);
saveas(h7,'ZStudyFig3PATL.jpg','jpg')

%%%%% Figure 4
h8 = figure('Position', [0, 0, 1440, 900]);
subplot(2,3,1)
histRatio = log10(bsAllEmp(~paraAllEmp))-avgBSResEmp(~paraAllEmp);
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
ylabel('Count {PATL}')
title('Emp. Free consumer:mean Resource BSR');
XL1 = xlim;

subplot(2,3,2)
histRatio = log10(bsAllEmp(paraAllEmp>0))- avgBSResEmp(paraAllEmp>0);
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. Parasite consumer:mean Resource BSR');
XL2 = xlim;

subplot(2,3,3)
histRatio = log10(bsAllEmp)- avgBSResEmp;
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. All Consumer:Resource BSR');
XL3 = xlim;

subplot(2,3,4)
 histRatio = log10(MAll(~paraAll,2))-avgBSResSim(~paraAll);
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Free consumer:mean Resource BSR');
%xlim(XL1);

subplot(2,3,5)
 histRatio = log10(MAll(paraAll>0,2))-avgBSResSim(paraAll>0);
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Para consumer:mean Resource BSR');
%xlim(XL2);

subplot(2,3,6)
 histRatio = log10(MAll(:,2))-avgBSResSim;
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. All consumer:mean Resource BSR');
%xlim(XL3);
saveas(h8,'ZStudyFig4PATL.jpg','jpg')

%%%%% Figure 9
h9 = figure('Position', [0, 0, 1440, 900]);
subplot(2,2,1)
histRatio = avgBSConEmpPara - log10(bsAllEmp) ;
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
ylabel('Count {SWTL}')
title('Emp. Free consumer:mean Resource BSR');
XL1 = xlim;

subplot(2,2,2)
histRatio =  avgBSConEmpFree - log10(bsAllEmp);
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. Parasite consumer:mean Resource BSR');
XL2 = xlim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
histRatio = avgBSConSimSWTLpara - log10(MAll(:,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. All Consumer:Resource BSR');
%xlim(XL1);

subplot(2,2,4)
 histRatio = avgBSConSimSWTLfree - log10(MAll(:,1));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Free consumer:mean Resource BSR');
%xlim(XL2);

saveas(h9,'ZStudyFig5SWTL.jpg','jpg')



%%%%% Figure 10
h10 = figure('Position', [0, 0, 1440, 900]);
subplot(2,2,1)
histRatio = avgBSConEmpPara - log10(bsAllEmp) ;
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
ylabel('Count {SWTL}')
title('Emp. Mean Parasite:mean Host BSR');
XL1 = xlim;

subplot(2,2,2)
histRatio =  avgBSConEmpFree - log10(bsAllEmp);
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Emp. Mean Free consumer: Resource BSR');
XL2 = xlim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
histRatio = avgBSConSimPATLpara - log10(MAll(:,2));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Mean Parasite:Host BSR');
%xlim(XL1);

subplot(2,2,4)
 histRatio = avgBSConSimPATLfree - log10(MAll(:,2));
hist(histRatio,20);
xlabel(sprintf('Mean = %.2f/Median = %.2f\nstd = %.2f/IQR = %.2f\nN=%u',...
    mean(histRatio,'omitnan'),...
    median(histRatio,'omitnan'),...
    std(histRatio,'omitnan'),...
    quantile(histRatio,.75)-quantile(histRatio,.25),...
    length(histRatio)))
title('Sim. Mean Free consumer:Resource BSR');
%xlim(XL2);

saveas(h10,'ZStudyFig5PATL.jpg','jpg')

%mean(log10(bs(bs(:,3)==1,2)./bs(bs(:,3)==1,1)),'omitnan')
%}


% title('All')
% figure;
% ratiosFree = ratiosAll(bs(:,3)==0);
% hist(ratiosFree)
% title('Free')
% figure;
% ratiosPara = ratiosAll(bs(:,3)==1);
% hist(ratiosPara)
% title('Para')