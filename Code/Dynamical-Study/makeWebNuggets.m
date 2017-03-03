%This code generates Web Nuggets.  I'll be adding parasites to these web
%nuggets.

%Starting with larger webs.  Hope is that we will land at around 40.  (Do I
%want to commit to 40?  Too late now...)  Could also test what changing C
%does, but to stay consistent I think keeping the same C is important.
%systematic changes in C tell us whether or more or less general species
%tend to go extinct.  Looking at extinction order is important too, in
%terms of generality/ Vulnerability.  (more or less general first, less or
%more later?  But at that point it can get so web dependent.)  Another
%thing to keep in mind is that the resulting webs are likely dependent on
%the choice of model.  SO, using different parameters *may* result in
%radically different web topologies.

%I also need to think about what I will get by having webs be ~40 and what
%the distribution of C values give.  Could possibly do regressions on that
%sort of thing since I'll have at least some variability there.  Eventually
%I will want to control for that by either factoring the effect out (via
%regression or something like that) or just using webs with S = 40 ( maybe
%+\- 1 to 4)
clear all
fclose('all');
S0 = 80;
C0 = 0.15;

%will add to this later; adding needs to be done carefully since I want to
%maintain some meta data files.  Will include room for up to 1,000,000 webs
%in my naming convention.


dataDir = '~/Desktop/JulySimulations/WebNuggets40';
OGWebDir = sprintf('%s/OGWebs',dataDir);
codeDir = pwd;

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



%Want to be careful about using this.  Worst goes to worst, the OG/WebIndex
%can be re-made.
try
    OGIndexFid = fopen(...
        sprintf('%s/webIndex.csv',OGWebDir)...
        ,'r');
    header = fgets(OGIndexFid);
    webIndex = fscanf(OGIndexFid,'%u,%u,%f,%f');
    webIndex = reshape(webIndex,4,[])';
    webStart = max(webIndex(:,1));
    fclose(OGIndexFid);
     OGIndexFid = fopen(...
        sprintf('%s/webIndex.csv',OGWebDir)...
        ,'a');
catch
     OGIndexFid = fopen(...
        sprintf('%s/webIndex.csv',OGWebDir)...
        ,'w');
    fprintf(OGIndexFid,'WebNumber,S0,C0,C\n');
    webStart = 0;
end


for ii = webStart:(webStart+nNewWebs-1);
    %Creating The intial web
    webBad = true;
    while webBad
        [res, con,n,c,r] = NicheModel_nk(S0,C0);
        simMx = calculateSimilarity(res,con);
        mx = sparse(res,con,1,S0,S0);
        webBad = max(max(simMx))==1;
    end
    OGNodeFilename = sprintf('%s/nodeData%06u.csv',OGWebDir,ii);
    OGLinkFilename = sprintf('%s/linkData%06u.csv',OGWebDir,ii);
    OGNodeFid = fopen(OGNodeFilename,'w');
    OGLinkFid = fopen(OGLinkFilename,'w');
    
    %Calculate global and local properties of the web -- this is of
    %trivial complexity compared to the dynamic part.  Could **maybe**
    %save some time by saving these results?  I don't think it's
    %necessarily worth it.
    
    localProperties = calculateLocalProperties(res,con);
    swtl = localProperties(:,10);
    sptb = localProperties(:,8)+1;
    patl = 2*swtl - sptb;
    
    basal = swtl == 1;
    nBasal = sum(basal);
    globalProperties = calculateGlobalProperties(res,con);
    C = length(res)/S0^2;
    
    %Now, need to setup the dynamic parameters.  Let's go ahead and just
    %hard code simulations for Z = 10,100,1000.
    
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
    B0 = .95*rand(S0,1)+.05;
    
    %Growth rate of basal species.  Not sure if this variability is
    %reasonable or not.
    gr = basal.*(randn(S0,1)*.1+1);
    
    %Saving the intial web & updating webIndex.
    fprintf(OGLinkFid,'%u,%u\n',[res,con]');
    fprintf(OGNodeFid,'%.9f,%.9f,%.9f,%.9f,%.9f\n',[n,c,r,gr,B0]');
    
    fprintf(OGIndexFid,'%06u,%03u,%.2f,%.9f\n',ii,S0,C0,C);
    
    fclose(OGLinkFid);
    fclose(OGNodeFid);
    
    
    for Z = [10,100]
        ZWebDir = sprintf('%s/Z%uNuggets',dataDir,Z);
        
       
        
        try
            cd(ZWebDir)
            cd(codeDir)
            
        catch
            mkdir(ZWebDir)
        end
        
        
         ZLinkFid = fopen(sprintf('%s/linkData%06u.csv',ZWebDir,ii),'w');
        ZNodeFid = fopen(sprintf('%s/nodeData%06u.csv',ZWebDir,ii),'w');
        
        try
            
            ZIndexFid = fopen(sprintf('%s/webIndex.csv',ZWebDir),'r');
            if ZIndexFid==-1
                error('webIndex for Z = %u does not exist.',Z)
            end
            fclose(ZIndexFid);
            ZIndexFid = fopen(sprintf('%s/webIndex.csv',ZWebDir),'a');
        catch
            ZIndexFid = fopen(sprintf('%s/webIndex.csv',ZWebDir),'w');
            fprintf(ZIndexFid,'WebNumber,S0,C0,C,Snew,Cnew,\n');

        end
        
        
        %Makes all basal species have a body size = 1.
        
        M = Z.^(patl - 1);
        
        %Prey-averaged trophic level does pretty well at reproducing
        %empirical body size ratio distributions.
        
        x = ax*M.^(-0.25);
        x(basal) = 0;
        y = 8*ones(S0,1);
        
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
            ,'B0',B0...
            ,'h',h...
            ,'x',x...
            ,'halfSat',halfSat...
            ,'Tf',Tf ...
            ,'extctThresh',1e-10...
            ,'AbsTol',1e-10...
            ,'RelTol',1e-7...
            ,'yij',yij...
            );
        
        sol = integrateToFilterWebs(res,con,param);
        
        fprintf(ZLinkFid,'%u,%u,%u\n',[res con sol.extctLink]');
        fclose(ZLinkFid);
        
        fprintf(ZNodeFid,'%.9f,%.9f,%.9f,%.9f,%u\n',...
            [n,c,r,sol.y(:,end),sol.extctOrder]');
        fclose(ZNodeFid);
        Snew = S0 - sum(sol.extctOrder<Inf);
        Lnew = sum(~sol.extctLink);
        Cnew = Lnew/(Snew^2);
        fprintf(ZIndexFid,'%u,%u,%.9f,%.9f,%u,%.9f\n',ii,S0,C0,C,Snew,Cnew);
        fclose(ZIndexFid);
    end
    
end

fclose(OGIndexFid);
