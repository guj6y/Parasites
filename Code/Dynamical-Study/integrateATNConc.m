function [sol] = integrateATNConc(res,con,p,concRes,concCon  )
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
           ,'numType',numType...
           ,'K',K...
           ,'para',para...
           ,'swtl',swtl...
           ,'eij',eij...
           ,'wij',wij...
           ,'yij',yij...
           ,'basal',basal...
           ,'x',x...
           ,'r',r...
           ,'halfSat',halfSat...
           ,'h',h...
           ,'Tf',Tf...
           ,'phi',phi...
            );

%}

extinctionThreshhold = 1e-9;

eij = sparse(res,con,1./p.eij,p.S,p.S);
wij = sparse(res,con,p.wij,p.S,p.S);
yij = sparse(res,con,p.yij,p.S,p.S);

B0h = p.halfSat^p.h;

one = ones(p.S,1);

    function [dB] = ATN_RHS(~,B)
        %trying for incremental increases in speed: (storing these in
        %memory vs. calculating and using on the fly?)
        B(B<0) = 0;
        
        F = (B.^p.h)*...
            (1./(one*B0h + wij'*(B.^p.h)))'.*wij;
%         F = ((~p.para+p.phi).*B.^p.h)*...
%         (1./(one*B0h + wij'*((~p.para+p.phi).*B.^p.h)))'.*wij;
        
        yF = yij.*F;
        xB = p.x.*B;
        
        dB = p.r.*B.*(1-p.basal'*B/p.K)...    %Basal logistic growth
            - xB...                                     %Metabolic losses
            +(yF'*one).*(xB)...     %increase from consumption
            -(yF.*eij)*(xB);     %Loss due to predation.
                


%         dB = p.r.*B.*(1-p.basal'*B/p.K)...    %Basal logistic growth
%         - xB...                                     %Metabolic losses
%         +(yF'*one).*xB...     %increase from consumption
%         -(yF.*eij)*(((~p.para)+p.phi).*xB);     %Loss due to predation.
        
        

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

options = odeset('Events',@extinction); ...,'nonNegative',one);
Tinterval = [tnow Tfinal];

B0 = p.B0;

sol = struct;
count = 0;
while tnow<Tfinal
    
    count = count+1;
    partName = sprintf('part%u',count);
    sol.(partName) = ode45(@ATN_RHS,Tinterval,B0,options);
    tnow = sol.(partName).x(end);
    
    
    Tinterval = [tnow Tfinal];
    extinct_sp = sol.(partName).y(:,end)<extinctionThreshhold;
    
    B0 = sol.(partName).y(:,end);
    B0(extinct_sp) = 0;
    
end
sol.n = count;


end


