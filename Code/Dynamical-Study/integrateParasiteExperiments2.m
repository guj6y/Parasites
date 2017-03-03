function [solOut] = integrateParasiteExperiments2(p)
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
wij = sparse(res,con,p.wij,p.S,p.S);
yij = sparse(res,con,p.yij,p.S,p.S);


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
yijTro = sparse(trophicLinkCon,trophicLinkRes,yijTroList,p.S,p.S);


eijParList = p.eij(paraParaLink);
eijTroList = p.eij(trophLink);

eijPar = sparse(paraLinkHost,paraLinkPara,1./eijParList,p.S,p.S);
eijTro = sparse(trophicLinkCon,trophicLinkRes,1./eijTroList,p.S,p.S);

wijParList = p.wij(paraParaLink);
wijTroList = p.wij(trophLink);

wijPar = sparse(paraLinkHost,paraLinkPara,wijParList,p.S,p.S);
wijTro = sparse(trophicLinkCon,trophicLinkRes,wijTroList,p.S,p.S);

    function [nadda, dB] = ATN_RHS(B)
        B(B<0) = 0;
        
        if includeFraction == 0
            F = (B.^p.h)*...
                (1./(one*B0h + wij'*(B.^p.h)))'*wij;
            
            yF = yij.*F;
            yFe = yF.*eij;
            xB = p.x.*B;
            
            dB = p.r.*B.*(1-p.basal'*B/p.K)... Basal logistic growth
                - xB...                        Metabolic losses
                +(yF'*one).*(xB)...            increase from consumption
                -(yFe)*(xB);               %Loss due to predation.
            
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
            FTro = (BPhih)*...
                (1./(one*B0h + wijTro'*(BPhih)))'*wijTro;
            FPar = (BPhih)*...
                (1./(one*B0h + wijPar'*(BPhih)))'*wijPar;
            
            yFTro = yijTro.*FTro;
            yFPar = yijPar.*FPar;
            
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
                
                %Multiply each row by the biomass of the free-living
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
        
        
                nadda = [];
        
    end


    function [value,isterminal,direction] = extinction(~,B)

            value = 1*(((abs(B)<extinctionThreshhold)&B>0)|B<0);

        isterminal = ones(p.S,1);
        direction = zeros(p.S,0);
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

options = odeset('Events',@extinction,'AbsTol',AbsTol,'RelTol',RelTol);
Tinterval = [tnow Tfinal];

B0 = p.B0;

sol = struct;
count = 1;
solOut.extctOrder = zeros(p.S,1);
B0s = zeros(p.S,p.S);

solOut = fmincon(@(x)0,B0,[],[],[],[],zeros(p.S,1),[],@ATN_RHS)


end


