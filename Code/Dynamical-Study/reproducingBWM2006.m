%Reproducing Brose et al. 2006
%Neo wanted to make sure that my model can reproduce results similar to
%those seen in Brose 2006.  I will be using my 'integrateATN.m' code to do
%this, allowing for appropriate functional responses.  This script ensures
%that we have the correct webs since they imposed the (rather funny, but
%necessary because of what they were comparing?) requirement the webs must
%have 5 basal species.  Another odd thing about the model is that the
%bodymasses of all the basal species were the same and had the same
%intrinsic growth rate?

%Moving forward from, I think I will follow the general strategy of
%relaxing/modifying the dynamics present in Brose, Williams, Martinez 2006,
%but still using that as my baseline.  In this way, I can more explicitly
%test the changes that are being introduced (and neo can always compare
%back to the results that he knows from this paper).

%%%%%%% Description of methods in Brose2006:
%%%%Food Webs:
% S = 20,40,60
% C = 0.15
% All webs have 5 basal species *****
%
%%%%Dynamics:
% (basal) dB/dt = r_i(M_i)G_iB_i - loss from consumption
% (cnsmr) dB/dt = -x_i(M_i)B_i + increase from Cnsmptn - loss from cnsmptn
%
% Loss from consumption:
% sum_{j=consumers} x_jy_jB_jF_{ji}/(e_{ji}f_{ji}
%
% Gain from Consumption:
% sum_{j=resources} x_iB_iF_{ij}
%
% Functional Response (i eating j):
% F_{ij} = w_{ij}B_j^h/(B_0^h + c_iB_iB_0^h) + sum_{k=rsrcs_i}w_{ik}B_k^h)
%
% Extinctions were handled at the end; any species with body mass less than
% 10^{-30}(!) was considered extinct.
% They don't say how they ran these simulations.  They were probably (well)
% handled by Rich Williams in C or C++.
%%%%Parameters:
%%% Functional Response:
%w_ij = 1/n (weak generalist; preference for i by j is 1/(# in j's diet)
%h = 1,2
%c = 0,1
%B_0 = 0.5
%f_{ji} = 1
% (Functional responses tested: (h,c) = (1,0);(2,0);(1,1))
%
%%%% Metabolic Scaling:
%R_p = a_rM_P^{-0.25}
%X_C = a_XM_C^{-0.25}
%Y_C = a_YM_C^{-0.25}
%M_C = Z^T
%a_x = .314(inv.) | .88(ecto.vert.)
%a_r = 1
%a_y = 2.51(inv.) | 3.5(ecto.vert.)
%
%%% Dynamical Equaiton parameters
%y_j = 4(ecto. Vert.) | 8(inv.)
%r_i = 1;
%x_i = a_x/a_r(M_C/M_P)^{-0.25}
%y_i = a_y/a_x
%e_{ij} = .85 (carnivory) | .45 (herbivory)
%K = 1
%


saveResults = true;
%average of the 6 empirical webs: S = 142, C = .082
%looking at the levels:
richnessAll = [20 30 40];
%connectanceAll = 0.15;
C = .15;
nBasal = 5;
%funcitonalResponsesAll(1,:) = h
hAll = [1 ...
    1.2...
    2];
nZ = 8;
ZAll = logspace(-2,5,nZ);

axAll = [.314,.88];
yAll = [8 4];
%S%h%ax%wij
models = fullfact([3,3,2,2]);
%models = fullfact([1,1,1,1]);
nModels = length(models);
halfSat = 0.5;


nWeb = 100;

dataDir = sprintf('/Volumes/Buster/Research/Parasitism/Dynamics/Brose2006Webs');
whichModel = 0;
tic
for model = models'
    whichModel = whichModel +1;
    
    S = richnessAll(model(1));
    h = hAll(model(2));
    ax = axAll(model(3));
    y = yAll(model(3));
    %The abundance level web directory
    
    dataDirS = sprintf('%s/S%03u',dataDir,S);
    
    %The abundance,connectance level web directory
    dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));
    resultsDir = sprintf('%s/ResultsPATL%u%u%u%u',dataDir,model(1),model(2),model(3),model(4));
    
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
        ii
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
                webBasal = sum((sum(mx)==0));
                webBad = (webBasal~=nBasal)|max(max(simMx))==1;
            end
            csvwrite(sprintf('%s/web%04u.csv',dataDirSC,ii-1),[res,con]);
        end
        %Calculate global and local properties of the web -- this is of trivial
        %complexity compared to the dynamic part.  Could **maybe** save some
        %time by saving these results?  I don't think it's necessarily worth
        %it.
        
        localProperties = calculateLocalProperties(res,con);
        swtl = localProperties(:,10);
        sptb = localProperties(:,8);
        patl = 2*swtl - sptb;
        basal = swtl == 1;
        globalProperties = calculateGlobalProperties(res,con);
        %We still need to do a few things--
        % x_i - determined from allometric scaling (Brose 2006)
        % r_i - follow distribution from Brose 2006
        % K - system wide carrying capacity from Brose 2006
        % e_ij - from Brose 2006
        % w_ij - weak generalist
        
        nFree = S-nBasal;
        
        yij = y * ones(size(res));
        
        %Assimilation efficiency is given in Brose et al as .45 and .85 for
        %herbivores and carnivores, respectively.  This would correspond to the
        %links we talkin bout.  Easy enough to figure:
        eij = zeros(size(res));
        
        eij(basal(res)) = .45;
        eij(~basal(res)) = .85;
        K = 1;
        
        %Strong generalist model for preference:
        if model(4) == 2
        wij = ones(size(res));
        elseif model(4) ==1
        %Weak generalist model for preference:
        dietSizes = zeros(S,1);
        for jj = 1:S
            dietSizes(jj) = sum(con==jj);
        end
        wij = 1./dietSizes(con);
        end
        
        
        try
            B0 = ...
                csvread(sprintf('%s/IC%04u.csv',dataDirSC,ii-1));
            
            
        catch MEmx
            
            
            B0 = .95*rand(S,1)+.05;
            
            csvwrite(sprintf('%s/IC%04u.csv',dataDirSC,ii-1)...
                ,B0);
        end
        
                    r = zeros(S,1);
            r(basal) = 1;
        x=zeros(S,1);
        param = struct( 'S',S...
            ,'C',C...
            ,'res',res...
            ,'con',con...
            ,'K',K...
            ,'swtl',swtl...
            ,'eij',eij...h
            ,'wij',wij...
            ,'yij',yij...
            ,'basal',basal...
            ,'r',r...
            ,'B0',B0...
            ,'h',h...
            ,'x',x...
            ,'halfSat',halfSat...
            ,'Tf',10e3 ...
            ,'phi',[]...
            );
        
        count = 0;
        
        yFinals = zeros(S,nZ);
        yMeans = zeros(S,nZ);
        yStds = zeros(S,nZ);
        clc
        curTime = toc;
        fprintf('You have completed %u of 36 Dynamic models.\n',whichModel-1)
        fprintf('You are currently integrating web %u with different values of Z.\n',ii)
        fprintf('Approximately %.2f minutes have elapsed.\n',toc/60)
        numDone = (whichModel-1)*100 + ii-1;
        numLeft = 100*length(models)-numDone;
        rate = curTime/numDone;
        timeLeft = rate*numLeft;
        fprintf('Approximately %.2f minutes are left (%.2f s per web).\n',timeLeft/60,rate);
        fprintf('Z = ')
        for Z = ZAll
            fprintf('| %.1e ',Z)
            count = count+1;
            
            M = (ones(S,1)*Z).^(patl-1);
            x = ax.*M.^(-0.25);
            x(basal) = 0;

            
            param.x = x;
            
            
            sol = integrateATN(res,con,param);
            
            totAlive(ii,count) = sum(sol.(sprintf('part%u',sol.n)).y(:,end)>1e-30);
            yFinals(:,count) = sol.yT;
            yMeans(:,count) = sol.mean;
            yStds(:,count) = sol.std;
            save(sprintf('%s/sol%04uZ%u.mat',resultsDir,ii-1,count),'-struct','sol');
        end
        if saveResults
            csvwrite(sprintf('%s/finals%04u.csv',resultsDir,ii-1),yFinals);
            csvwrite(sprintf('%s/means%04u.csv',resultsDir,ii-1),yMeans);
            csvwrite(sprintf('%s/stds%04u.csv',resultsDir,ii-1),yStds);
           
        end
    end
end
% figure
% plot(sol.(sprintf('part%u',sol.n)).x,sol.(sprintf('part%u',sol.n)).y)
% [sol.(sprintf('part%u',sol.n)).y(:,end) basal r x]







