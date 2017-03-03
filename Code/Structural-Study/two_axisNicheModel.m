function [res,con,n,c,r,failures,C_all] = ...
    two_axisNicheModel(Nf,Np,C4,para_above)

%TODO: (maybe?) only redraw subwebs; c, r, or n?

%This function creates a connected niche web that matches 4 connectances
%corresponding to all possible feeding relationships between 2 distinct
%classes of taxa.  In this case, one class of taxa uses reverse feeding rules;
%(Inverse niche model) that is, the (parasitic) taxa eats above it's niche
%range, and (parasitic) taxa lower on the niche axis are mroe generalist.

%defining error tolerances for each of the subwebs
err_ff = .15;
err_pf = .15;
err_fp = .15;
err_pp = .15;
err_all = .05;

%Total number of species:
N=Nf+Np;

%Connectance for f->f:
Cff = C4(1);

%Connectacnes for p->f:
Cpf = C4(2);

%Connectance for f->p:
Cfp = C4(3);

%Connectances for p->p:
Cpp = C4(4);

C = (Cff*Nf^2 + (Cpf + Cfp)*(Nf*Np) + Cpp*Np^2)/(N^2);



%parameters for beta random variables:
a = 1;
Bff = a/2/Cff -a;
Bpf = a/2/Cpf -a;
Bfp = a/2/Cfp -a;
Bpp = a/2/Cpp -a;


%Loop maintanence

tries = 10000;
ok_n=0;
connected = 0;

failures = zeros(6,1);
unconnected_species = (1:N)';
Cs = zeros(1001,5);
while (ok_n==0) && (tries >0) && (connected ==0)
    %these are for counting where webs fail.
    fail_con = 0;
    ff_fail = 0;
    fp_fail = 0;
    pf_fail = 0;
    pp_fail = 0;
    all_fail = 0;
    connected = 0;
    
    %Niche Values
    nf = rand(Nf,1);
    np = rand(Np,1);
    
    %beta scaling for ranges
    yff = betarnd(a,Bff,Nf,1);
    ypf = betarnd(a,Bpf,Nf,1);
    yfp = betarnd(a,Bfp,Np,1);
    ypp = betarnd(a,Bpp,Np,1);
    
    %yff = 2*Cff;
    %ypf = 2*Cpf;
    %yfp = 2*Cfp;
    %ypp = 2*Cpp;
    
    %Feeding ranges
    rff = nf.*yff;
    rpf = nf.*ypf;
    rfp = (1-np).*yfp;
    rpp = (1-np).*ypp;
    
    %Want to ensure all ranges are within niche space.
    rf = max(rff,rpf);
    rp = max(rfp,rpp);
    
    %diet centers
    cf = rand(Nf,1).*min(nf-rf./2,1-rf)+rf./2;
    if para_above
        cp = rand(Np,1).*min(1-np-rp./2,1-rp)+np;
    else
        cp = rand(Np,1).*(1-rp)+rp./2;
    end
    
    %sorting the vectors.
    [nf_sort, Indx_f] = sort(nf);
    [np_sort, Indx_p] = sort(np);
    
    rff_sort = rff(Indx_f);
    rpf_sort = rpf(Indx_f);
    rfp_sort = rfp(Indx_p);
    rpp_sort = rpp(Indx_p);
    
    cf_sort = cf(Indx_f);
    cp_sort = cp(Indx_p);
    
    %ensuring we have a basal species with lowest niche value
    rff_sort(1) = 0;
    rpf_sort(1) = 0;
    
    
    %constructing matrices to determine feeding relationships.
    onesf = ones(Nf,1);
    onesp = ones(Np,1);
    %ff submatrix
    Rff1 = onesf*((cf_sort-rff_sort/2)');
    Rff2 = onesf*((cf_sort+rff_sort/2)');
    
    %pf submatrix
    Rpf1 = onesp*((cf_sort-rpf_sort/2)');
    Rpf2 = onesp*((cf_sort+rpf_sort/2)');
    
    %fp submatrix
    Rfp1 = onesf*((cp_sort-rfp_sort/2)');
    Rfp2 = onesf*((cp_sort+rfp_sort/2)');
    
    %pp submatrix
    Rpp1 = onesp*((cp_sort-rpp_sort/2)');
    Rpp2 = onesp*((cp_sort+rpp_sort/2)');
    
    %upper and lower bounds matrices
    R1_mx = [Rff1 Rfp1;Rpf1 Rpp1];
    R2_mx = [Rff2 Rfp2;Rpf2 Rpp2];
    
    %nichevalue matrix
    N_mx = [nf_sort;np_sort]*ones(1,N);
    
    web_mx = (R1_mx<N_mx) & (N_mx < R2_mx);
    
    
    %Checking each of the subweb's connectances:
    
    %f->f connectance
    Lff = sum(sum(web_mx(1:Nf,1:Nf)));
    Cff_obs = Lff/Nf^2;
    Cs(tries,1) = Cff_obs;
    if abs((Cff_obs - Cff)/Cff)<err_ff
        ok_ff = 1;
    else
        ok_ff = 0;
        ff_fail = 1;
    end
    
    %p-f> connectance
    Lpf = sum(sum(web_mx(Nf+(1:Np),1:Nf)));
    Cpf_obs = Lpf/(Nf*Np);
    Cs(tries,2) = Cpf_obs;
    if abs((Cpf_obs - Cpf)/Cpf)<err_pf
        ok_pf = 1;
    else
        ok_pf = 0;
        pf_fail = 1;
    end
    
    %f_>p connectance
    Lfp = sum(sum(web_mx(1:Nf,Nf+(1:Np))));
    Cfp_obs = Lfp/(Nf*Np);
    Cs(tries,3) = Cfp_obs;
    if abs((Cfp_obs - Cfp)/Cfp)<err_fp
        ok_fp = 1;
    else
        ok_fp = 0;
        fp_fail = 1;
    end
    
    %f_>p connectance
    Lpp = sum(sum(web_mx(Nf+(1:Np),Nf+(1:Np))));
    Cpp_obs = Lpp/(Np*Np);
    Cs(tries,4) = Cpp_obs;
    if abs((Cpp_obs - Cpp)/Cpp)<err_pp
        ok_pp = 1;
    else
        ok_pp = 0;
        pp_fail = 1;
    end
    
    if (ok_ff)&&(ok_pf)&&(ok_fp)&&(ok_pp)
        ok_subs = 1;
    else
        ok_subs = 0;
        
    end
    
    %Checking overall connectances
    links = sum(sum(web_mx));
    % Actual connectance
    C_web = links/(N^2);
    Cs(tries,5) = C_web;
    if ((abs(C_web-C)*1.0/C) < err_all)
        ok_all=1;
    else
        ok_all=0;
        all_fail = 1;
    end
    
    ok_n = ok_subs&&ok_all;
    if ok_n
        %Make sure that all species are weakly connected to the specified 
        %basal species (food chains cross); but only if the connectances
        %are okay.  This is too time intensive to check every time.
        unconnected_species = (1:N)';
        
        connected_species = [];
        test_mx = (web_mx+web_mx')>0;
        
        connected_species = walk(1,connected_species,test_mx);
        unconnected_species(connected_species) = [];
        
        if isempty(unconnected_species)
            web_connected = 1;
        else
            web_connected = 0;
            fail_con = 1;
        end
        
        unconnected_species = (1:N)';
        basal_species = unconnected_species(sum(web_mx)==0)';
        connected_species = [];
        if sum(basal_species>Nf)
            web_connected = 0;
        elseif sum(...
                sum(web_mx(Nf+1:end,Nf+1:end))>0 &... %hyperparasites have
                (sum(web_mx(:,Nf+1:end)) ... equal diet size on parasites
                == sum(web_mx(Nf+1:end,Nf+1:end))))...  as total diet size
                %then...
            %Parasites have no free-living prey; they ain't parasites.
            web_connected = 0;
        else
            
            
            for kk = basal_species
                connected_species = walk(kk,connected_species,web_mx);
            end
            
            unconnected_species(connected_species) = [];
            
            if isempty(unconnected_species)
                web_connected = 1;
            else
                web_connected = 0;
                fail_con =1 ;
            end
        end
        
    end
    if ~ok_n
        failed= [fail_con;ff_fail;pf_fail;fp_fail;pp_fail;all_fail];
        failures = failed+failures;
    end
    ok_n = ok_n&&web_connected;
    tries = tries-1;
end
c = [cf_sort;cp_sort];
n = [nf_sort;np_sort];
r = [rff_sort rpf_sort;rfp_sort rpp_sort];
if ok_n ==1
%Calculate Cff:
linksff = sum(sum(web_mx(1:Nf,1:Nf)));
Cff = linksff/(Nf^2);

%Calculate Cfp:
linksfp = sum(sum(web_mx(1:Nf,(Nf+1):(Np+Nf))));
Cfp = linksfp/(Np*Nf);

%Calculate Cpf:
linkspf = sum(sum(web_mx((Nf+1):(Np+Nf),1:Nf)));
Cpf = linkspf/(Np*Nf);

%Calculate Cpp:
linkspp = sum(sum(web_mx((Nf+1):(Np+Nf),(Nf+1):(Np+Nf))));
Cpp = linkspf/(Np*Np);

C_all = [C_web, Cff, Cpf, Cfp, Cpp];    
    [res,con] = find(web_mx);
    10000-tries
    failures
else
    res = [];
    con = [];
    
    
end

end