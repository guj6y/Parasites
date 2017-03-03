    %This script transforms data from http://esapubs.org/archive/ecol/... and
    %creates food webs (networks of trophic interactions between species).  I
    %want to consistenly aggregate and select the various species.  It saves
    %the created webs as link-lists:
    %
    % res,con,para
    %  i,j,b
    %
    %Here, i is the resource index, j is the consumer index, and b is a binary
    %variable that is 1 if the link is parasitic (i not a parasite, j is a
    %parasite).  I also save the web properties from calculateLocalProperites.

    %I had to do quite a bit of post-processing on the files inorder to make
    %them somewhat consistent and not a nightmare to import (the
    %post-processing was a nightmare in itself, but whatevs)


    %Detailed descriptions of the column headings and variables are available
    %in the metadata.

    %Requires a somewhat recent

    linkFiles = {'BSQweb_Links.txt',...
        'CSMweb_Links.txt',...
        'EPBweb_Links.txt',...
        'Flensburg_Data_Links.txt',...
        'Otago_Data_Links.txt',...
        'Sylt_Data_Links.txt'};

    numColLinks = [20,20,20,16,16,16];

    nodeFiles = {'BSQweb_Nodes.txt',... Has body sizes
        'CSMweb_Nodes.txt',...          Has body sizes
        'EPBweb_Nodes.txt',...          Has body sizes
        'Flensburg_Data_Nodes.txt',...  No body sizes
        'Otago_Data_Nodes.txt',...      NO body sizes
        'Sylt_Data_Nodes.txt'};        %No body sizes

    numColNodes = [46,46,46,41,41,41];

    n_webs = 6;
    
    if ~exist('concomittantWebs')
        concomittantWebs = false;
    end
    
    if concomittantWebs
        fprintf('Generating webs with concomittant links\n')
    else
        fprintf('Generating webs without concomittant links\n')
    end
    topoConcomittant = false;
    codeDir = pwd;

    shortWebNames = {'bahia','carp','punta','flens','otago','sylt'};

    %These are the links that I include in the food web.  We skip non-trophic
    %interactions, but it might be important to double-check which links I use.
    % Descriptions are found in the metadata.
    linkTypesTrophic = [1 ...   predation
        ,2 ...                  Social Predation
        ,3 ...                  MicroPredation
        ,18 ...                 Detritivory
        ,21 ...                 Facultative Micropredation
        ];

    linkTypesTrophicPara = [16 ...  Predation on free-living non-feeding stage
        ,19 ...                     Parasite Intraguild Antagonism
        ];

    linkTypesPara = [4 ...  Parasitic Castration
        ,5 ...              Pathogen Infection
        ,6 ...              Macroparasitism
        ,8 ...              Pollination
        ,10 ...             Trophically Transmitted Parasitic Castration
        ,11 ...             Trophically Transmitted pathogen infection
        ,12 ...             Trophically Transmitted Parasitism
        ];

    concLinks = [14 ...     Concurrent predation on symbionts
        ,15 ...             Trophic Transmission
        ];

    if concomittantWebs
        linkTypesIncluded = [linkTypesTrophic linkTypesTrophicPara...
            linkTypesPara concLinks];
    else
        linkTypesIncluded = [linkTypesTrophic linkTypesTrophicPara...
            linkTypesPara];
    end



    propertiesCell = cell(n_webs,2);
    oldNamesCell = cell(6,1);
    newNamesCell = cell(6,1);

    oldClassCell = cell(6,1);
    newClassCell = cell(6,1);

    nEndo = zeros(6,1);
    totSpecies = zeros(6,1);
    totPara = zeros(6,1);
    totBasal = zeros(6,1);
    totFreeCon = zeros(6,1);

    Ls = zeros(6,1);
    connectances = zeros(6,1);
    binEndoCell = cell(6,1);
    binEctoInvertCell = cell(6,1);
    binEctoVertCell = cell(6,1);
    binAutoCell = cell(6,1);
    binParaCell = cell(6,1);
    speciesTypeCell = cell(6,1);
    meanEndoPara = zeros(6,1);

    nodeDatas = cell(6,1);
    linkDatas = cell(6,1);
    linkListCell = cell(6,1);
    aggregationLevel = 1;
    jjLinksNInitial = zeros(6,length(linkTypesIncluded)+1);
    linkTypesCell = cell(6,1);
    linkParaCell = cell(6,1);

    meanGen = zeros(6,2);
    varGen = zeros(6,2);
    seGenPool = zeros(6,1);
    for ii = 1:6
        cd('../../Data/Raw/')
        cd('Links')
        linkData = readtable(linkFiles{ii},'Delimiter','tab');

        if ii>3
            linkData.Properties.VariableNames(7) = {'LinkTypeID'};
        end

        linkDatas{ii} = linkData;
        %linksIncluded is a binary array that indicates the links to include
        %in the initial food web.
        linksIncluded = zeros(size(linkData.LinkTypeID));

        %Detritius?

        %These Nodes are detritus - I need to be sure to take them out
        %completely since I was getting detritus as a basal species (!?) - some
        %links of snails eating detritus were labeled as predation; 4 webs have
        %detritus, 3 label 'detritus', the  other as 'Detritus'
        %detritusNodes = strcmp(nodeData.OrganismalGroup,'Detritus')|...
        %    strcmp(nodeData.OrganismalGroup,'detritus');

        % detritusNodeIDs = nodeData.NodeID(detritusNodes);
        % detritusSpeciesIDs = nodeData.SpeciesID(detritusNodes);

        %linksPara is a binary array that indicates whether the link type
        %implies a parasitic consuemr
        linksPara = linksIncluded;
        linksConc = linksIncluded;
        %Extract the links that I want; see metadata.
        count = 0;
        for jj = linkTypesIncluded
            count = count+1;
            jjLinks = linkData.LinkTypeID==jj;
            linksIncluded = (linksIncluded) + (jjLinks)*jj;
            jjLinksNInitial(ii,count) = sum(jjLinks);
            if sum(linkTypesPara==jj)
                linksPara = (linksPara) | (jjLinks);
            end
            if jj == 14
                linksConc = linksConc | jjLinks;
            end

        end
        jjLinksNInitial(ii,end) = sum(jjLinksNInitial(ii,1:(end-1)),2);
        %List of IDs of resources, consumers, and parasites as species; res0
        %and con0 define the foodweb; para0 will become a parasitic species
        %identifier.
        res0 = linkData.ResourceSpeciesID(linksIncluded>0);
        con0 = linkData.ConsumerSpeciesID(linksIncluded>0);
        para0 = linkData.ConsumerSpeciesID(linksPara);

        %Trying to add all conceivable concomittant links.  This tries to add
        %about 16K links for the first web... which is a *bit* high.  Don't do
        %this.
        if topoConcomittant
            resToAdd = [];
            conToAdd = [];
            for jj = unique(para0')
                jjsHosts = res0(con0==jj);
                if numel(jjsHosts >0)
                    for kk = jjsHosts'
                        kksPred = con0(res0==kk);
                        concAdded = [];
                        for ll = kksPred'
                            if sum(ll==jj)>0
                                continue
                            elseif sum(ll==jjsHosts) >0
                                continue
                            else
                                concAdded = [concAdded; ll];
                            end
                        end
                        resToAdd = [resToAdd; ones(size(concAdded))*jj];
                        conToAdd = [conToAdd; concAdded];
                    end
                end

            end
            res0 = [res0; resToAdd];
            con0 = [con0; conToAdd];
        end
        %Identifying Concomittant links; still not sure about intermediate host
        %at this point.  Before the following two lines, these are the same
        %size as the original link list from data.
        linksConc = linksConc(linksIncluded>0);
        linksPara = linksPara(linksIncluded>0);

        %Now, linksConc is a binary array of the same length as the initial
        %link list (i.e. it is L x 1 vector), =1 if the link is concomittant,
        %and zero otherwise.  Similarly for linksPara.

        %linkIDs Specifies the type of each link.  Could very well use the
        %original link types from the data, too!  Do this for now since I'm not
        %overly concerned about all those different types.
        linkTypes = linksIncluded(linksIncluded>0);
        linkIDs =  ones(size(linksConc)) + linksPara + 2*linksConc;

        %Actually Removing detrital links
        %for jj = detritusSpeciesIDs'
        %
        %    jjDetritalLinks = (res0==jj)|(con0==jj);
        %    res0(jjDetritalLinks) = [];
        %    con0(jjDetritalLinks) = [];
        %    para0(para0==jj) = [];

        %end


        %Unique IDs of resource, consumers, as stages:
        %resStage0 = linkData.ResourceNodeID(linksIncluded);
        %conStage0 = linkData.ConsumerNodeID(linksIncluded);
        %paraStage0 = linkData.ConsumerNodeID(linksPara);
        %
        %Actually Removing Detrital links ( and therefore species)
        %for jj = detritusNodeIDs'
        %    jjDetritalLinks = (res0==jj)|(con0==jj);
        %    resStage0(jjDetritalLinks) = [];
        %    conStage0(jjDetritalLinks) = [];
        %    paraStage0(paraStage0==jj) = [];
        %
        %end



        %Read node data
        cd('../Nodes')
        nodeData = readtable(nodeFiles{ii},'Delimiter','tab');
        cd('..');
        nodeDatas{ii} = nodeData;
        nodeData = sortrows(nodeData,'NodeID');

        %I hate using unique, but..
        idList = unique([res0;con0]);
        bodySizes = zeros(size(idList));
        biomasses = zeros(size(idList));

        %Doing the bodysizes here.
        for jj = 1:length(idList)
            %The body size of species j is the adult body size of species j.
            %The biomass of species j is the sum of the biomasses of all
            %stages.  NaN are omitted; biomass not that important, mainly as a
            %decent initial condition for simulations on these webs. If the
            %adult does not have a body size, take the body size of the most
            %abundant stage. If the most abundant stage does not have a body
            %size, take the lowest stage IDs body size.  Hopefully the only
            %thing that this leaves as NaN are the few messed up species, and
            %the plants.
            j = idList(jj);

            binj = nodeData.SpeciesID==j;
            minStage = min(nodeData.StageID(binj));
            binMin = nodeData.StageID==minStage;

            bodySizes(jj) = nodeData.BodySize_g_((binj)&(binMin));

            biomasses(jj) = sum(nodeData.Biomass_kg_ha_(binj),'omitnan');
            if sum(binj) > 1
                if isnan(bodySizes(jj))&&sum(~isnan(nodeData.Biomass_kg_ha_(binj))>0)
                    bodySizes(jj) = nodeData.BodySize_g_((binj)&...
                        (nodeData.Biomass_kg_ha_==...
                        max(...
                        nodeData.Biomass_kg_ha_(binj&(~isnan(nodeData.BodySize_g_)))...
                        )...
                        )&~isnan(nodeData.BodySize_g_));
                end


            end
        end


        %The total number of species in the Data set; the final number included
        %in my foodweb may be less than this!!  Actually it's not excactly
        %that for carpinteria, punta and bahia because of how their data is
        %saved.
        nAllSpecies = max(nodeData.SpeciesID);

        %SEtting up names
        oldNamesCell{ii} = cell(nAllSpecies,1);
        oldNamesCell{ii}(nodeData.SpeciesID) = nodeData.WorkingName;

        oldClassCell{ii} = cell(nAllSpecies,1);
        oldClassCell{ii}(nodeData.SpeciesID) = nodeData.Class;

        %oldNamesCell{ii}(~binSpeciesIncluded) = [];
        %List of all species Ids in the datase
        allSpeciesList = 1:nAllSpecies;

        %Binary array that indicates whether a particular species out of the
        %entire data set is included in my food web.  %WOrking with species Ids
        binSpeciesIncluded = false(nAllSpecies,1);
        binSpeciesIncluded([res0;con0]) = true;

        %List of only the old species Ids that are included in my food web.
        oldIdsInWeb = allSpeciesList(binSpeciesIncluded);

        %Number of species in new web:
        nSpecies = length(oldIdsInWeb);
        %List of new IDs in web
        newIds = 1:nSpecies;

        %New names in web
        newNamesCell{ii} = cell(nSpecies,1);
        newClassCell{ii} = cell(nSpecies,1);

        %The replacement vector works by placing the new ID at the old ID's
        %position.
        replacementVector = zeros(nAllSpecies,1);
        replacementVector(oldIdsInWeb) = newIds;

        %replacementVector('oldID of species j') = 'new ID of species j'
        %so,this re-does all the IDs correctly.  This doesn't change the length
        %of res or con, so now the link list has redundancies.  It doesn't
        %really matter that much since we can use 'sparse' to get the unique
        %link list.
        res = replacementVector(res0);
        con = replacementVector(con0);

        %Note that we don't need to worry link IDs since the lengths aren't
        %changing, so we still have the links correctly identified.  However,
        %there is still the question of whether or not concomittant links get
        %grouped with non-concomittant links.


        %Doesn't work so nicely for names,sizes, and biomasses.
        newNamesCell{ii} = oldNamesCell{ii}(oldIdsInWeb);
        newClassCell{ii} = oldClassCell{ii}(oldIdsInWeb);

        %Binary array that indicates whether a particular species is a
        %parasite; out of all species in the dataset
        binParaAll = false(nAllSpecies,1);
        binParaAll(para0) = true;

        %From that, only taking the species that are actually included.
        para = binParaAll(binSpeciesIncluded);

        %Is the final web a *trophic* web?
        trophicWeb = false;
        [LL, idxLL] = sortrows([con res]);
        con = LL(:,1);
        res = LL(:,2);
        linkTypes = linkTypes(idxLL);
        cd(codeDir)
        %Until the web is trophic
        while ~trophicWeb
            %CAUTION: this is gonna combine parasites with free-livers.
            %Probably not something we want to do.  It doesn't happen in the
            %six original ecosystems I looked at (punta, bahia, carpinteria,
            %sylt, flensburg, otago), but still a concern if this code gets
            %recycled.  I highly doubt it would happen, but you never know...

            %Is that ^ not an indication that they are a functionally distinct
            %class of species?>?

            %similarity Matrix
            simMx = calculateSimilarity(res,con);
            %Same matrix:  (sim(i,j) = 1, if i is the same as j.& =0 o.w.
            sameMx = triu(simMx>=aggregationLevel,1);
            if (sum(sum(sameMx)))==0
                trophicWeb = true;
                continue
            end
            %we want I -> J; I lists indices of superfluous species.
            [I,J] = find(sameMx);

            if sum(para(I) == ~para(J)) >0
                fprintf('You combined a parasite and a free-liver!!\n')
            end
            %number of species in the web.
            nSpecies = max([res;con]);

            %Species that we need to delete; the triu means that doing this
            %won't delete any species equivalence classes (trophic species)
            deleteVector = (1:nSpecies)';
            deleteVector(I) = 0;

            %Delete links that include a superfluous species; first, change
            %the species to zero; res0 and con0 are 1 if the link is to be
            %deleted, and zero otherweise.
            res0 = deleteVector(res)==0;
            con0 = deleteVector(con)==0;


            %Then, delete all links that contain a zero; must do this to both
            %the resource and consumer vector, as well as linkIDs; species are
            %aggregated by trimming links.
            res(res0|con0) = [];
            con(res0|con0) = [];
            linkTypes(res0|con0) = [];

            %My worry here is that a concomittant (or parsitic) link gets
            %aggregated with a non-concomittant (or parsitic) link.  That
            %manifests by deleting the particular superfluous link with the
            %link ID that I (may) want to keep.  Could go through the same
            %rigamarole as I do with the species names for the links, but that
            %would require another loop in here and honestly... I am too Lazy
            %to figure that out.
            %Adding concomittant links could cause more aggregation
            %(species that prey on parasites and hosts of parasites could look
            %like species that prey on hosts of parasites and concomittantly
            %prey on their parasites themselves.  But is that an important
            %distinction?)
            linkIDs(res0|con0) = [];

            %Delete the slots corresponding to superflous species.
            para(I) = [];

            %Renewing IDs - this is where the number of species gets trimmed
            oldIds = deleteVector(deleteVector>0);
            newIds = 1:length(oldIds);

            %The replacement vector works by placing the new ID at the old ID's
            %position.
            replacementVector = deleteVector;
            replacementVector(oldIds) = newIds;

            res = replacementVector(res);
            con = replacementVector(con);
            linkPara = (para(con))&(~para(res));

            %We also need to update the species names and body sizes.
            %This is not as straightforward as above, since all equivalent
            %species need to be combined in some way.

            %Combining the names
            %Averaging the body sizes.  Weight is the total biomass of each
            %species.
            I0 = I;
            J0 = J;
            for kk = I'
                %Equivalence class of the first species in array I:
                equivClasskk = [kk;J(I==kk)];


                %now, average the bodySizes.
                bodySizes(equivClasskk(end)) = sum(biomasses(equivClasskk)...
                    .*bodySizes(equivClasskk),'omitnan')/...
                    (sum(biomasses(equivClasskk).*...
                    (~isnan(bodySizes(equivClasskk))),'omitnan'));

                biomasses(equivClasskk(end)) = sum(biomasses(equivClasskk).*...
                    (~isnan(bodySizes(equivClasskk))),'omitnan');

                bodySizes(equivClasskk(1:end-1)) = 0;
                biomasses(equivClasskk(1:end-1)) = 0;
                for ll = equivClasskk(1:end-1)';
                    %Need the if statement in case equivClasskk is empty.
                    if numel(ll)>0
                        %Saving names
                        newNamesCell{ii}(equivClasskk(end)) = ...
                            strcat(newNamesCell{ii}(equivClasskk(end)),...
                            ';',newNamesCell{ii}(ll));
                        newClassCell{ii}(equivClasskk(end)) = ...
                            strcat(newClassCell{ii}(equivClasskk(end)),...
                            ';',newClassCell{ii}(ll));
                        J(I==ll) = [];
                        I(I==ll) = [];
                    end
                end


            end

            bodySizes(I0) = [];
            biomasses(I0) = [];
            newNamesCell{ii}(I0) = [];
            newClassCell{ii}(I0) = [];
        end
totSpecies(ii) = length(newClassCell{ii});

        MxTypes = sparse(res,con,linkTypes,totSpecies(ii),totSpecies(ii));
        [res,con,linkTypes] = find(MxTypes);
        linkPara = para(con)&(~para(res));
        linkParaCell{ii} = linkPara;
        linkListCell{ii} = [res con];
        
        totPara(ii) = sum(para);
        connectances(ii) = length(res)/totSpecies(ii)^2;
        Ls(ii) = length(res);

        properties = zeros(length(bodySizes),12);
        properties(:,1:11) = calculateLocalProperties(res,con);
        properties(:,12) = bodySizes;
        properties(:,13) = biomasses;
        propertiesCell{ii,1} = properties;
        propertiesCell{ii,2} = para;

        binParaCell{ii} = para;
        binEndoCell{ii} = zeros(totSpecies(ii),1);
        binEctoVertCell{ii} = zeros(totSpecies(ii),1);
        binEctoInvertCell{ii} = zeros(totSpecies(ii),1);
        binAutoCell{ii} = zeros(totSpecies(ii),1);


        %Ectotherm Vertebrates: one of these classes
        binAgna = ~cellfun('isempty',strfind(newClassCell{ii},'Agnatha'));
        binChon = ~cellfun('isempty',strfind(newClassCell{ii},'Chondrichthyes'));
        binOsth = ~cellfun('isempty',strfind(newClassCell{ii},'Ostheichthyes'));
        binOste = ~cellfun('isempty',strfind(newClassCell{ii},'Osteichthyes'));
        binAmph = ~cellfun('isempty',strfind(newClassCell{ii},'Amphibia'));
        binRept = ~cellfun('isempty',strfind(newClassCell{ii},'Reptilia'));
        binEctoVertCell{ii} = binAgna|binChon|binOste|binAmph|binRept|binOsth;

        %Endotherms: one of these classes
        binAvesCell = ~cellfun('isempty',strfind(newClassCell{ii},'Aves'));
        binMammCell = ~cellfun('isempty',strfind(newClassCell{ii},'Mammalia'));
        binEndoCell{ii} = binAvesCell|binMammCell;
        %sum(binEndoCell{ii}&binEctoVertCell{ii})
        %Autotrophs: generality = 0
        binAutoCell{ii} = properties(:,2)==0;
        totBasal(ii) = sum(binAutoCell{ii});
        totFreeCon(ii) = totSpecies(ii) -totBasal(ii) - totPara(ii);

        %sum(binEndoCell{ii}&binAutoCell{ii})
        %sum(binEctoVertCell{ii}&binAutoCell{ii})
        %Ectotherm invertebrates: the rest.  There are tiny little critters
        %that maybe shouldn't be classified as ectotherm invertebrates (i.e.
        %single cell organisms).. there aren't very many nodes, but might they
        %be important enough to delineate?
        binEctoInvertCell{ii} = ...
            ~(binEctoVertCell{ii}|binEndoCell{ii}|binAutoCell{ii});

        speciesTypeCell{ii} = nan(totSpecies(ii),1);
        speciesTypeCell{ii}(binEctoInvertCell{ii}) = 1;
        speciesTypeCell{ii}(binEctoVertCell{ii}) = 2;
        speciesTypeCell{ii}(binEndoCell{ii}) = 3;
        speciesTypeCell{ii}(binAutoCell{ii}) = 0;

        linkTypesCell{ii} = linkTypes;

        if concomittantWebs
            cd('../../Data/Processed')
            fid = fopen(sprintf('%sWebCon.csv',shortWebNames{ii}),'w');
            fprintf(fid,'res,con,para,con\n');
            fprintf(fid,'%u,%u,%u,%u\n',[res,con,linkTypes]');
            fclose(fid);
            cd(codeDir)

            fid = fopen(sprintf('%sPropertiesCon.csv',shortWebNames{ii}),'w');
            %                1         2    3         4             5           6               7           8               9           10      11
            %varNames = {'clustCoef','gen','vul','meanVulPrey','meanImpPrey','meanGenPred','meanImpPred','minSPToBasal','numConnBasal','SWTL','inLoop'};
            fprintf(fid,'\n');
            fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[properties para]');
            fclose(fid);

        elseif topoConcomittant
            continue
        else
            cd('../../Data/Processed')
            fid = fopen(sprintf('%sWeb.csv',shortWebNames{ii}),'w');
            fprintf(fid,'res,con,para\n');
            fprintf(fid,'%u,%u,%u\n',[res,con,linkPara]');
            fclose(fid);
            cd(codeDir)

            fid = fopen(sprintf('%sProperties.csv',shortWebNames{ii}),'w');
            %                1         2    3         4             5           6               7           8               9           10      11
            %varNames = {'clustCoef','gen','vul','meanVulPrey','meanImpPrey','meanGenPred','meanImpPred','minSPToBasal','numConnBasal','SWTL','inLoop'};
            fprintf(fid,'\n');
            fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[properties para]');
            fclose(fid);
        end
    end





