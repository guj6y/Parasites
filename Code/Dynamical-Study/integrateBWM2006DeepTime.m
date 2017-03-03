function [sol] = integrateBWM2006DeepTime(res,con,p)
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

extinctionThreshhold = p.extctThresh;

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

        yF = yij.*F;
        xB = p.x.*B;
        
        dB = p.r.*B.*(1-p.basal.*B/p.K)...    %Basal logistic growth
            - xB...                                     %Metabolic losses
            +(yF'*one).*(xB)...     %increase from consumption
            -(yF.*eij)*(xB);     %Loss due to predation.
                



        
        

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


B0 = p.B0;

sol = struct;
count = 0;

snapshots = [5000,10000,25000,50000,100000,250000,500000,1e6];
TDeep = 1e6;
maxT = 1e4;

sol.snaphsots = zeros(p.S,length(snapshots));
sol.extinctions = zeros(0,3);
SList = 1:p.S;
snapshotsSaved = 1;

while tnow<Tfinal
    Tinterval = linspace(tnow,min(tnow+maxT,Tfinal),1000);
    
    partName = sprintf('part%u',count);
    sol.(partName) = ode45(@ATN_RHS,Tinterval,B0,options);
    tnow = sol.(partName).x(end);
    
    snapshotsInInterval = snapshots((snapshots>Tinterval(1))&(snapshots<=tnow));
    nSnapshotsInInterval = numel(snapshotsInInterval);
    
    if nSnapshotsInInterval>0
        sol.snapshots(:,snapshotsSaved:(snapshotsSaved+nSnapshotsInInterval-1)) = ...
            deval(sol.(partName),snapshotsInInterval);
        snapshotsSaved = snapshotsSaved+nSnapshotsInInterval;
    end
    %Tfinal = tnow+1000; This *could* be a good idea
    
    extinct_sp = (sol.(partName).ye<extinctionThreshhold)&...
        (sol.(partName).ye~=0);
    if sum(extinct_sp)>0
        count = count+1;
    fprintf('Extinction event %u at t = %.3f.  %u species remain.\n',count,tnow,p.S-count)
    extinctSpIdx = SList(extinct_sp)';
    extinctSpMetRate= p.x(extinct_sp);
    B0 = sol.(partName).y(:,end);
    B0(extinct_sp) = 0;

    try
    sol.extinctions = [sol.extinctions;
        sol.(partName).xe.*ones(size(extinctSpIdx)),extinctSpIdx,extinctSpMetRate];
    catch
        fprintf('WibblyWobbly goddammit!!\n')
    end
    else
        fprintf('The current time is t = %.4e\n',tnow)
        B0 = sol.(partName).y(:,end);
    end
    
end
sol.snapshots(:,end) = sol.(partName).y(:,end);
sol.n = count;

xMin = sol.(partName).x(1);
xMax = sol.(partName).x(end);
xRange = ceil(max(xMin,xMax-1000)):floor(xMax);

y = deval(sol.(partName),xRange);
sol.nx = numel(xRange);
sol.yT = y(:,end);
sol.mean = mean(y(:,floor(end/2+1):end),2);
sol.std = std(y(:,floor(end/2+1):end),0,2);


end


