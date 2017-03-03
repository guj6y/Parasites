concomittantWebs = false;
%run('../Parasite-Classification/webGeneration');
run('../Parasite-Classification/BigStats');
close all
pSig = 0.05;
genWebsAll = zeros(2,6);
vulWebsAll = zeros(2,6);

genWebsAllC = zeros(2,6);
vulWebsAllC = zeros(2,6);
nWeb = 100;
for web = 1:6
S = totSpecies(web);
    C = connectances(web);
    Sp = totPara(web);
    gensFree = zeros(1,nWeb);
    vulsFree = zeros(1,nWeb);
    gensCFree = zeros(1,nWeb);
    vulsCFree = zeros(1,nWeb);
    
    gensPara = zeros(1,nWeb);
    vulsPara = zeros(1,nWeb);
    gensCPara = zeros(1,nWeb);
    vulsCPara = zeros(1,nWeb);
    for ii = 1:nWeb
    
    para = false(S,1);
    [res,con,~,~,~] = NicheModel_nk(totSpecies(web),connectances(web));
    localProp = calculateLocalProperties(res,con);
    basal = localProp(:,2) == 0;
paraNum = datasample(1:S,Sp,'replace',false);
    para(paraNum) = true;
    
    gensPara(:,ii) = mean(localProp(para,2));
    vulsPara(:,ii) = mean(localProp(para,3));
    
    gensFree(:,ii) = mean(localProp(~(para|basal),2));
    vulsFree(:,ii) = mean(localProp(~(para|basal),3));
    
    
    [conc] = identifyPotentialConcomittant(res,con,para);
    resC = [res;conc(:,1)];
    conC = [con;conc(:,2)];
    localPropC = calculateLocalProperties(resC,conC);
    
    gensCPara(:,ii) = mean(localPropC(para,2));
    vulsCPara(:,ii) = mean(localPropC(para,3));
    
    gensCFree(:,ii) = mean(localPropC(~(para|basal),2));
    vulsCFree(:,ii) = mean(localPropC(~(para|basal),3));
    end
genWebsAll(1,web) = mean(gensFree);
genWebsAll(2,web) = mean(gensPara);
vulWebsAll(1,web) = mean(vulsFree);
vulWebsAll(2,web) = mean(vulsPara);


genWebsAllC(1,web) = mean(gensCFree);
vulWebsAllC(1,web) = mean(vulsCFree);
genWebsAllC(2,web) = mean(gensCPara);
vulWebsAllC(2,web) = mean(vulsCPara);
end

fig= figure('Position',[0,0,1440,900]);
subplot(2,2,1);
gscatter(meanPropFreeCon(:,2),meanPropPara(:,2),pProp(:,2)<pSig,'kr','o')
line(genWebsAll(1,:),genWebsAll(2,:),'Marker','x','LineStyle','none');
line([meanPropFreeCon(:,2) genWebsAll(1,:)']',...
     [meanPropPara(:,2) genWebsAll(2,:)']');
axis([0,2,0,2])
refline(1,0)
xlabel('Mean Generality of Free-living consumers')
ylabel('Mean Generality of Parasites')
legend('Insignificant at P=0.05','Significant at P=0.05','NicheModels')
title('Generality of Empirical Webs (no concomittant)')

subplot(2,2,2);

gscatter(meanPropFreeCon(:,3),meanPropPara(:,3),pProp(:,3)<pSig,'kr','o')
line(vulWebsAll(1,:),vulWebsAll(2,:),'Marker','x','LineStyle','none');
line([meanPropFreeCon(:,3) vulWebsAll(1,:)']',...
     [meanPropPara(:,3) vulWebsAll(2,:)']');
axis([0,2.5,0,2.5])
refline(1,0)
xlabel('Mean Vulnerability of Free-living consumers')
ylabel('Mean Vulnerability of Parasites')
legend('Insignificant at P=0.05','Significant at P=0.05','NicheModels')
title('Vulnerability of Empirical Webs (no concomittant)')

concomittantWebs = true;
%run('../Parasite-Classification/webGeneration');
run('../Parasite-Classification/BigStats');

figure(fig);
subplot(2,2,3);

gscatter(meanPropFreeCon(:,2),meanPropPara(:,2),pProp(:,2)<pSig,'kr','o')
line(genWebsAllC(1,:),genWebsAllC(2,:),'Marker','x','LineStyle','none');
line([meanPropFreeCon(:,2) genWebsAllC(1,:)']',...
     [meanPropPara(:,2) genWebsAllC(2,:)']');
axis([0,2,0,2])
refline(1,0)
xlabel('Mean Generality of Free-living consumers')
ylabel('Mean Generality of Parasites')
legend('Insignificant at P=0.05','Significant at P=0.05','NicheModels')
title('Generality of Empirical Webs (concomittant)')

subplot(2,2,4);

gscatter(meanPropFreeCon(:,3),meanPropPara(:,3),pProp(:,3)<pSig,'kr','o')
line(vulWebsAllC(1,:),vulWebsAllC(2,:),'Marker','x','LineStyle','none');
line([meanPropFreeCon(:,3) vulWebsAllC(1,:)']',...
     [meanPropPara(:,3) vulWebsAllC(2,:)']');
axis([0,2.5,0,2.5])
refline(1,0)
xlabel('Mean Vulnerability of Free-living consumers')
ylabel('Mean Vulnerability of Parasites')
legend('Insignificant at P=0.05','Significant at P=0.05','NicheModels')
title('Vulnerability of Empirical Webs (concomittant)')