function [res,con,n,c,r,C_all,parasites]= ...
     NicheModelPara_test(N, C, Np, para_above)                
%------------------------------------------------------------------

%globalStream = RandStream.getGlobalStream;
%reset(globalStream);
Nf = N - Np;
tries=10000;                                                               
%error on the connectances
error_n=0.05;

%validation of the plain niche web:
ok_n=0; 

%We have to be careful about counting links.  What I have the same for
%parasites and free-livers in this code is the same average number of link
%per species.  Ensuring that parasites have a host would normally drive
%their number of links per species up

%Total links:
L = C*N^2;

%Expected links for parasites:
Lp = L*Np/N;

%Expected links for free-livers:
Lf = L*Nf/N;

%Since parasites are automatically given a host, this is the expected
%number of new links. hopefully it is a positive number... (-_-;) 
%back of the envelope says it should be except for very sparse webs.
Lp = Lp - Np; 

%Connectance of random free-liver links
Cf = Lf/(Nf^2 + Np*Nf);

%Connectace of random, additional parasitic links.
Cp = Lp/(Np^2 +Np*Nf);

%Randomly designate parasites (ensures parasites *can* always have a host
if para_above
    para = datasample(1:(N-1),Np,'Replace',false);
else
    para = datasample(2:N,Np,'Replace',false);
end

parasites = zeros(N,1);
parasites(para) = true;
parasites = parasites>0;
freeLivers = ~parasites; 


c = zeros(N,1);
r = zeros(N,1);
while ((N>0) && (tries>0) && (ok_n==0))
    tries=tries-1;
    %assign niche values from a uniform distribution
    n = rand(N,1);  
    n = sort(n);
    nf = n(freeLivers);
    np = n(parasites);
    cp = c(parasites);


%vector of ranges: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha = 1;
    if para_above
        Enp = (N/(N+1))/2;
        betap = alpha*(1-Enp)/Cp - alpha;
        rp = betarnd(alpha,betap,Np,1); 
        rp = (1-np).*rp;
    else
        Enp = (1+1/(N+1))/2;
        betap = alpha*Enp/Cp-alpha;
        rp = betarnd(alpha,betap,Np,1); 
        rp = rp.*np;
    end
    
    betaf = alpha/(2*Cf)-alpha; 
    
    rf = betarnd(alpha,betaf,Nf,1); 
    rf(1) = 0;
    rf = rf.*nf;
    r(freeLivers) = rf;
    r(parasites) = rp;
    % set center of range, uniformly distributed in [r_i/2,n_i]; 
    cf = min(1-rf./2,rand(Nf,1).*(nf-rf./2)+rf./2);
    
    if para_above ==1
        
        %For each parasite, pick a free-living niche value lower than
        %itself as the center of its range.  Makes sense since parasites 
        %should usually be pretty highly specialized to their hosts?  If
        %not, could also wiggle a bit to more closely mimic the niche model
        
        for jj = 1:Np
            cp(jj) = datasample(nf(nf>=np(jj)),1);
            
        end
        
    else


        for jj = 1:Np
            
            cp(jj) = datasample(nf(nf<=np(jj)),1);
        end
    end
        %Fixing if too high or too low center value - diet fully within niche
    %axis.
    over = (cp+rp/2)>1;
    cp(over) = 1-rp(over)/2;
    under  = (cp-rp/2)<0;
    cp(under) = rp(under)/2;
    %sort everything:
    c(parasites) = cp;
    c(freeLivers) = cf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    
    %lower border of niche range for every prey:
    preymins = c - r/2;
    
    %upper border of niche range for every predator:
    preymaxs = c + r/2;

    %fills the empty matrix with niche ranges:
    n_mx = n*ones(1,N); 
    %matrix with the lowest points of ranges in every column:
    preymins_mx = ones(N,1)*preymins'; 
    %same, with highest:
    preymaxs_mx = ones(N,1)*preymaxs';

    %Construct the web adjacency matrix;
    %if species in the row is in the diet range of the species in the 
    %column, it gets eaten (if species i is eaten by species j, set (i,j) =
    %1).
    
    web_mx=(n_mx>=preymins_mx)&(n_mx<=preymaxs_mx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  if there is an isolated species                                    %
    %  or something not connected to a basal                               %
    %  or a disconnected group                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %Make sure the graph is weakly connected (no isolated species)
    weak_comp = graphconncomp(sparse(web_mx),'Directed',1,'Weak',1);
    web_connected = (weak_comp == 1);
    
    %make sure that all species have a directed connection from a basal
    %species; only if connected.
   
     %make sure that all species have a directed connection from a basal
    %species; only if connected.
    if web_connected
        unconnected_species = (1:N)';
        basal = sum(web_mx)==0;
        basal_species = unconnected_species(basal);
        %if a parasite is a basal species, it has no host; no good
        if sum(basal(parasites))>0
            web_connected = 0;
        elseif sum(sum(web_mx(parasites,parasites))==sum(web_mx(:,parasites)))>0
            %Parasites have no free-living prey; they ain't parasites.
            web_connected = 0;
        else
            %All parasites have a (non-parasitic) host; check that all
            %energy flow starts at a basal
            connected_species = [];
            
            for kk = basal_species'
                connected_species = walk(kk,connected_species,web_mx);
            end
            
            unconnected_species(connected_species) = [];
            
            if isempty(unconnected_species)
                web_connected = 1;
            else
                web_connected = 0;
            end
        end
    end
    
    
    if web_connected
        % indices of links (which element in web_mx?)
        links = sum(sum(web_mx));  
        % Actual connectance
        C_web = links/(N^2);  

        if (abs(C_web-C)*1.0/C) > error_n
            ok_n=0;
        else
            ok_n=1;
        end
    end  

end


%Calculate Cff:
linksff = sum(sum(web_mx(freeLivers,freeLivers)));
Cff = linksff/(Nf^2);

%Calculate Cfp:
linksfp = sum(sum(web_mx(freeLivers,parasites)));
Cfp = linksfp/(Np*Nf);

%Calculate Cpf:
linkspf = sum(sum(web_mx(parasites,freeLivers)));
Cpf = linkspf/(Np*Nf);

%Calculate Cpp:
linkspp = sum(sum(web_mx(parasites,parasites)));
Cpp = linkspp/(Np*Np);

links = sum(sum(web_mx));
C_web = links/N^2;

C_all = [C_web,Cff,Cpf,Cfp,Cpp];

    [res,con] = find(web_mx);

    
end