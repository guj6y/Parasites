function [Res,Cons,n,c,r]= ...
     NicheModelParasites_(num_species2, C4,p_range)
%OUTPUT 
%Res = source half of link list
%Cons = target half of link list
%n = niche values for free livers and parasites.  The first n are free
%livers the last m are parasites.
%c = diet centers
%r = diet ranges

%INPUT
%num_species2 = ["#free livers" "#parasites"]
%C4 = [Cff Cpf Cfp Cpp]; (see below)
%p_range = range of parasites in niche space (they are restricted in this
%model)
 
%Needs:
%walk1.m
%checkConnected.m
%makeRawWeb_para.m

%Called by: 
%properties_experiments.m


%TODO (NK):  Need to streamline this.  Might be a good idea to make a separate
%function to generate the web w/out modified diet widths.  That chunk
%occurs three times in this code.  Also, a function that wraps up the
%connectance checks with walk1.


%globalStream = RandStream.getGlobalStream;
%reset(globalStream);

%third iteration; redraws all species if there is a problem.

np_min = p_range(1);
np_max = p_range(2);
Enp = (np_min+np_max)/2;

%Calculate expected width of diet 
Nf = num_species2(1);
Np = num_species2(2);

%Connectance for f->f:
Cff = C4(1);

%Connectacnes for p->f:
Cpf = C4(2);

%Connectance for f->p:
Cfp = C4(3);

%Connectances for p->p:
Cpp = C4(4);

%Overall connectance:
C = (Cff*Nf^2 + (Nf*Np)*(Cfp+Cpf) + Np^2*Cpp)/(Np+Nf)^2;


Cs = C4;
tries=10000;                                                               
%error on the connectance for free-livers
error_n=0.05;
err_sub = 0.075;
%error on the connectance for parasites

%validation of the plain niche web:
ok_n=0; 

num_species = Nf+Np;
N = num_species;



while ((num_species>0) && (tries>0) && (ok_n==0))
    tries=tries-1;
    
    
    %makeRawWeb_para generates a web that aims for all proper connectances:
    %f->f,f->p,p->f,p->p
    [web_mx,n,c,r,Cpf_new,Cpp_new] = makeRawWeb_para(Nf,Np,Cs,p_range);   
    Cs = [Cff Cpf_new Cfp Cpp_new];
    %Make sure that all species are connected properly
    [web_connected,unconnected_species] = checkConnected(web_mx,N,Nf,Np);
    
    %If web isn't connected, we start over- but with targeted connectances.
    
    if web_connected == 0
        
        continue
    end

    
    %need to make sure that the free-liver web has the right connectance.
    web_free = web_mx(1:Nf,1:Nf);
    
    
    if web_connected
        links = sum(sum(web_free));
        C_web_ff = links/Nf^2;
        errff = (abs(C_web_ff-Cff)*1.0/Cff);
        if  errff > err_sub
            sub_ok=0;
            web_connected = 0;
        else
            sub_ok=1;
        end
    end
    
    if web_connected
        links = sum(sum(web_mx(1:Nf,(Nf+1):(Np+Nf))));
        C_web_fp = links/Nf/Np;
        errfp = (abs(C_web_fp-Cfp)*1.0/Cfp);
        if errfp > err_sub
            sub_ok=0;
            web_connected = 0;
        end
    end
    
    if web_connected
        links = sum(sum(web_mx((Nf+1):(Np+Nf),1:Nf)));
        C_web_pf = links/Nf/Np;
        errpf = (abs(C_web_pf-Cpf)*1.0/Cpf);
        if errpf > err_sub
            sub_ok=0;
            web_connected = 0;
        end
    end
    
    if web_connected
        links = sum(sum(web_mx((Nf+1):(Np+Nf),(Nf+1):(Np+Nf))));
        C_web_pp = links/Np^2;
        errpp = (abs(C_web_pp-Cpp)*1.0/Cpp);
        if errpp > err_sub
            sub_ok=0;
            web_connected = 0;
        end
    end
    
    if web_connected && sub_ok
        % indices of links
        links = sum(sum(web_mx));  
        % Actual connectance
        C_web = links/(num_species^2);  
        err = (abs(C_web-C)*1.0/C);
        if  err > error_n
            ok_n=0;
            
        else
            ok_n=1;

        end
    end  
    
end
    [Res,Cons] = find(web_mx);
end