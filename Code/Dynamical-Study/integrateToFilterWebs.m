function [sol] = integrateToFilterWebs(res,con,p)
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

B0h = p.halfSat^p.h;
one = ones(p.S,1);
extinctionThreshhold = p.extctThresh;


eij = sparse(res,con,1./p.eij,p.S,p.S);
wij = sparse(res,con,p.wij,p.S,p.S);
yij = sparse(res,con,p.yij,p.S,p.S);


    function [dB] = ATN_RHS(~,B)
        B(B<0) = 0;
        
        %There must be a more efficient way to do this.. 
        F = (B.^p.h)*...
            (1./(one*B0h + wij'*(B.^p.h)))'.*wij;
        
        yF = yij.*F;
        yFe = yF.*eij;
        xB = p.x.*B;
        
        dB = p.r.*B.*(1-p.basal'*B/p.K)... Basal logistic growth
            - xB...                        Metabolic losses
            + (yF'*one).*(xB)...            increase from consumption
            - (yFe)*(xB);               %Loss due to predation.
    end




    function [value,isterminal,direction] = extinction(~,B)
        if sum((abs(B)<extinctionThreshhold)&B>0)||sum(B<0);
            value = 0;
        else
            value = 1;
        end
        isterminal = 1;
        direction = 0;
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
count = 0;
extct = zeros(p.S,1);
extctOrder = zeros(p.S,1);

while tnow<Tfinal
    
    if numel(fieldnames(sol)) == 0
        sol = ode45(@ATN_RHS,Tinterval,B0,options);
    else
        sol = odextend(sol,@ATN_RHS,Tfinal,B0,options);
    end
    
    tnow = sol.x(end);
    
    
    %Tfinal = tnow+1000; This *could* be a good idea
    if numel(sol.xe)>0
        if (Tfinal - sol.xe(end))<2000
            Tfinal = max(Tfinal,ceil(tnow+2000));
        end    
    end
    
    extinct_sp = sol.y(:,end)<extinctionThreshhold;
    if numel(sol.xe)>0
        if tnow == sol.xe(end);
            count = count+1;
            newExtct = extinct_sp&~extct;
            extctOrder(newExtct) = count;
            extct = extinct_sp;
        end
    else
        extinct_sp = zeros(S,1);
    end
    B0 = sol.y(:,end);
    B0(extinct_sp) = 0;
    
end
extctOrder(~extct) = inf;
sol.extctOrder = extctOrder;
extctRes = extct(res);
extctCon = extct(con);
extctLink = extctRes|extctCon;
sol.extctLink = extctLink;
end