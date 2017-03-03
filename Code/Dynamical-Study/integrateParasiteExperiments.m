function [solOut] = integrateParasiteExperiments(p)
%This program runs the dynamical simulations for ATN with parasites.
%{
This is the model: (Boit et al 2012 without detritus or extra mortality
and! Also incorporating distribtuions from Brose 2006 and scalings
from WIlliams 2007: Homage to yodzis and Innes...
dB_i/dt = -f_m x_i B_i  maintenance
            + f_a x_i B_i\sum_j[y_{ij}F_{ij}(B) Benefit from eating
            - \sum_j (x_j Y_ji B_j F_ji(B))/e_{ji} Loss to being eaten

For basal species the first two terms are replaced by
r_iB_iG_i(B) -- (exudation fraction?)

F_ij as in Boit 2012, no interference terms:
 
w_ij Bj^q/(B_0 + \sum_resources w_il B_l^q)

G_i as in Boit 2012:

1- sum_produceers B_j/K

p is an object with all the parameters needed (don't want to calculate
properties in this function; also, want to have access to these after code
runs.
p. ...
             'S',S...
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
            ,'extctThresh',1e-30...
            ,'AbsTol',1e-15...
            ,'RelTol',1e-5...
            ,'modelCode',modelCode...
            );


%}

%just use all of these and identify that as a potential drawback.
%Note that we don't need to identify concomittant links at all if we
includeFraction = p.modelCode(1)==2;
includeConcomittant = p.modelCode(2)==2;

res = p.res;
con = p.con;

phi = p.para.*p.phi + ~p.para;
B0h = p.halfSat^p.h;
one = ones(p.S,1);
extinctionThreshhold = p.extctThresh;


eij = sparse(res,con,1./p.eij,p.S,p.S);
eij2 = sparse(res,con,p.eij,p.S,p.S);
wij = sparse(res,con,p.wij,p.S,p.S);
yij = sparse(res,con,p.yij,p.S,p.S);
yeij = yij.*eij;

paraCon = p.para(con);

conM = p.M(con);
resM = p.M(res);

paraParaLink = paraCon & (conM<resM);

paraLinkPara = con(paraParaLink);
paraLinkHost = res(paraParaLink);

trophicLinkCon = con(~paraParaLink);
trophicLinkRes = res(~paraParaLink);

trophLink = ~paraParaLink;

yijParList = p.yij(paraParaLink);
yijTroList = p.yij(trophLink);

yijPar = sparse(paraLinkHost,paraLinkPara,yijParList,p.S,p.S);
yijTro = sparse(trophicLinkRes,trophicLinkCon,yijTroList,p.S,p.S);


eijParList = p.eij(paraParaLink);
eijTroList = p.eij(trophLink);

eijPar = sparse(paraLinkHost,paraLinkPara,1./eijParList,p.S,p.S);
eijTro = sparse(trophicLinkRes,trophicLinkCon,1./eijTroList,p.S,p.S);

yeijPar = yijPar.*eijPar;
yeijTro = yijTro.*eijTro;

wijParList = p.wij(paraParaLink);
wijTroList = p.wij(trophLink);

wijPar = sparse(paraLinkHost,paraLinkPara,wijParList,p.S,p.S);
wijTro = sparse(trophicLinkRes,trophicLinkCon,wijTroList,p.S,p.S);

    function [dB] = ATN_RHS(~,B)
        B(B<0) = 0;
        
        if includeFraction == 0
            Bh = B.^(1+p.h);
           
            %This is not as efficient is it strictly could be.
              yF = (Bh)*...
                 (1./(one*B0h + wij'*(Bh)))'.*yij;
             yFmx =yF;
             yFemx = yF.*eij2;
             
            %Just trying to avoid doing forming the full F matrix (no
            %sparsity pattern).  This doesn't really do much for S = 40.
            %Worth remembering if S is much larger; I think it saves time
            %there but for S = 40 this ends up taking longer.  Oh well.
%             den = (one*B0h + wij'*(Bh))';
%             yF = Bh(res)./den(con)'.*p.yij;
%             yFe = yF./p.eij;
%             yFmx = sparse(res,con,yF,p.S,p.S);
%             yFemx = sparse(res,con,yFe,p.S,p.S);
%             
%             
            xB = p.x.*B;
            
            dB = p.r.*B.*(1-p.basal'*B/p.K)... Basal logistic growth
                - xB...                        Metabolic losses
                +(yFemx'*one).*(xB)...            increase from consumption
                -(yFmx)*(xB);               %Loss due to predation.
            
            if includeConcomittant == 1
                
                %total losses due to trophic consumption of host(i.e.
                %non-parasitic consumer)
                totLossHost = yFe*(xB.*(~p.para));
                
                %Fraction of each p's total biomass in host h
                fhp = bsxfun(@times,yF,1./sum(yF));
                
                %I like this orientation better, for no reason ata ll
                fph = fhp';
                
                
                %Multiply each row by parasite population (that row's %B)
                fph = bsxfun(@times,fph,B.*(p.para));
                %Now fph is the amount of biomass of parasite p in host h.
                %THis zeros out any row tha tis not a parasite, which is
                %good.
                
                %Divide each column by biomass of that host (that col's B)
                fph = bsxfun(@times,fph,1./B');
                %now fph is the amount of biomass of parasite p per unit of
                %biomass of host h
                
                %Could be getting inf or nan; fixing those
                fph(~isfinite(fph)) = 0;
                
                %Should now be this simple; just making sure that we don't
                %change non-parasites..  Those rows hsould already be
                concomittantModifier = (fph*totLossHost);
                
                dB = dB - concomittantModifier.*p.para;
            end
            
        else
            %This model separates out trophic and parasitic links for
            %parasites. It probably takes about twice as long
            
            BPhi = B.*phi;
            BPhih = BPhi.^p.h;
            %need functional responses for each type of interaction:
            %trophic or parasitic.
            yFTro = (BPhih)*...
                (1./(one*B0h + wijTro'*(BPhih)))'.*yijTro;
            yFPar = (BPhih)*...
                (1./(one*B0h + wijPar'*(BPhih)))'.*yijPar;

            
            yFeTro = yFTro.*eijTro;
            yFePar = yFPar.*eijPar;
            
            xB = p.x.*B;
            
            dB = p.r.*B.*(1-p.basal'*B/p.K)...  Basal logistic growth
                - xB...                         Metabolic losses
                +(yFTro'*one).*(xB.*phi)...     trophic consumption
                +(yFPar'*one).*(xB.*(1-phi))... parasitic consumption
                -(yFeTro)*(xB.*phi)...          Loss due to predation.
                -(yFePar)*(xB.*(1-phi));       %Loss due to parasitism.
            
            
            if includeConcomittant == 1
                
                %Total losses due to trophic interactions for each host.
                totLossHost = yFeTro*(xB.*phi);
                
                %Fraction of each parasites p's total biomass in host h
                fhp = bsxfun(@times,yFPar,1./sum(yFPar));
                
                
                %Same, but makes more sense to me;
                %rows are parasites, columns are hosts
                fph = fhp';
                
                %Multiply each row by the biomass of the non free-living
                %portion of that parasite population (that row's B)
                fph = bsxfun(@times,fph,B.*(1-phi));
                %Now fph is the amount of biomass of parasite p in host h.
                %THis zeros out any row tha tis not a parasite, which is
                %good.
                
                %Divide each column by biomass of that host (that col's B)
                fph = bsxfun(@times,fph,1./B');
                %now fph is the amount of biomass of parasite p per unit of
                %biomass of host h
                
                %Could easily get division by zero; just set those to zero
                fph(~isfinite(fph)) = 0;
                
                %Multiply f_(p,h)by the total loss of each host, h and sum
                %the result: this is the matrix multiplication.
                concomittantModifier = (fph*(totLossHost));
                
                %Subtract it off.
                dB = dB - concomittantModifier;
            end
            
        end
        
        
        
        
        
        
    end

    function [value,isterminal,direction] = extinction(~,B)
        
        value = (log10(B)-log10(p.extctThresh));%*100; %/10
        value(B<0) = -1;
        isterminal = ones(p.S,1);
        direction = zeros(size(value));
    end


tnow = 0;
Tfinal = p.Tf;

try
    AbsTol = p.AbsTol;
    RelTol = p.RelTol;
catch
    AbsTol = extinctionThreshhold;
    RelTol = 1e-5;
end

options = odeset(p.options,'Events',@extinction,'AbsTol',AbsTol,'RelTol',RelTol);
Tinterval = [tnow Tfinal];

B0 = p.B0;

sol = struct;

solOut.extctOrder = zeros(p.S,1);
B0s = zeros(p.S,p.S);
B0s(:,1) = B0;
solOut.extctTime = zeros(1,p.S);
extinct_sp = false(p.S,1);
% opts = optimset('display','none');
%
% extct = inf(p.S,1);
% lsqnonlin(@ATN_RHS2,B0,zeros(p.S,1),extct);
% %[B0,~,~,~,~,~,~] = fmincon(@(x)0,B0,[],[],[],[],zeros(p.S,1),[],@ATN_RHS2,opts);
% extct((B0<1e-10)&(B0>0))=0;
% count = 0;
% while sum((B0<1e-10)&(B0>0))>0
%     %[B0,~,~,~,~,~,~] = fmincon(@(x)0,B0,[],[],[],[],zeros(p.S,1),[],@ATN_RHS2,opts);
%     lsqnonlin(@ATN_RHS2,B0,zeros(p.S,1),extct);
%     count = count+1;
% end
countExtct = 1;
countSol = 0;
while tnow<Tfinal
    countSol = countSol+1;
    %could just use odextend
    sol = p.odeSolver(@ATN_RHS,Tinterval,B0,options);
    tnow = sol.x(end);
    
    %Saving all solutions for diagnostic purposes.  Not necessary to save
    %the entire time series? The values at each extinction would be a good
    %compromise; can easily reconstruct particular parts of each sim. that
    %I want (but still, time-consuming to get the entire time series..).  I
    %just don't have the storage right now to make this work.
    
    solName = sprintf('sol%u',countSol);
    
    
    %Tfinal = tnow+1000; This *could* be a good idea
    if numel(sol.ie)>0
        Tfinal = max(Tfinal,ceil(tnow+2000));
    end
    %We need to catch species that get completely synchronized.  check:
    %correlation with extinct species?
    %relative error between extinct species?
    Tinterval = [tnow Tfinal];
    
    extctIndices = sol.ie;
    extinct_sp0 = extinct_sp;
    extinct_sp(extctIndices) = true;
    
    
    B0 = sol.y(:,end);
    
    
    bDead = B0(sol.ie);
    
    for ii = bDead'
        relErrori = abs(ii-B0)/ii;
        extinct_sp(relErrori<RelTol)=true;
    end
    
    extctIndices = find(extinct_sp - extinct_sp0);
    
    B0(extinct_sp) = 0;
    
    sol.ie = extctIndices;
    
    solOut.extctOrder(sol.ie) = countExtct;
    solOut.extctTime(countSol+1) = sol.x(end);
    solOut.(solName) = sol;
    
    for ii = 1:numel(extctIndices)
        countExtct = countExtct+1;
        B0s(:,countExtct+1) = B0;
    end
    
end
solOut.extctTime(countSol+2:end) = [];
B0s(:,countSol+1) = sol.y(:,end);
B0s(:,countSol+2:end) = [];
solOut.B0s = B0s;
solOut.n = countSol;

xMax = floor(sol.x(end));
xMin = xMax - 999;
xRange = linspace(xMin,xMax,1000);

y = deval(sol,xRange);
solOut.y  = y(:,end);
solOut.mean = mean(y,2);
solOut.std = std(y,0,2);


end


