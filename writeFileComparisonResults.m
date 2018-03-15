function writeFileComparisonResults()
    varFnameVars1 = 'carloSimulatorVars';
    varFnameVars2 = 'lsaSimulatorVars';
    
    errTol = 0.0000001;
    
    varfname = 'comparisonResults.txt';
    fileId = fopen(varfname,'w');
    
    varsSimu1 = load(varFnameVars1);
    varsSimu2 = load(varFnameVars2);
    
    % Seed Number
    if varsSimu1.seedNumber == varsSimu2.seedNumber
        fprintf(fileId,'Seed Number: OK\r\n');
    else
        fprintf(fileId,'Seed Number: error\r\n');
    end
    
    fprintf(fileId,'\r\n');
    fclose(fileId);
    
    
    % cell positions
    
    nCells1 = length(varsSimu1.varcellpos_x);
    nCells2 = length(varsSimu2.varcellpos_x);
    
    if nCells1 == nCells2
        compareElemProperty(varfname,varsSimu1.varcellpos_x, ...
            varsSimu2.varcellpos_x,'Coordinate x of cells','Cell',errTol);

        compareElemProperty(varfname,varsSimu1.varcellpos_y, ...
            varsSimu2.varcellpos_y,'Coordinate y of cells','Cell',errTol);
    else
        fileId = fopen(varfname,'a');
        fprintf(fileId,'Number of Cells are not the same\r\n');
        fclose(fileId);
    end
    
    % ue positions
    
    nUes1 = length(varsSimu1.varuepos_x);
    nUes2 = length(varsSimu2.varuepos_x);
    
    if nUes1 ~= nUes2
        fileId = fopen(varfname,'a');
        fprintf(fileId,'Number of UEs are not the same\r\n');
        fclose(fileId);
    else
        compareElemProperty(varfname,varsSimu1.varuepos_x, ...
            varsSimu2.varuepos_x,'Coordinate x of UEs','UE',errTol);

        compareElemProperty(varfname,varsSimu1.varuepos_y, ...
            varsSimu2.varuepos_y,'Coordinate y of UEs','UE',errTol);
        
        compareElemProperty(varfname,varsSimu1.varcellue, ...
           varsSimu2.varcellue,'Serving cell of UEs','UE');
    end        
    
    % distance between cells and UEs
    
    compareVarsCellsByUsers(varfname,varsSimu1.vardistance,varsSimu2.vardistance,...
        'Distance comparison',nCells1,nUes1,errTol);

    compareVarsCellsByUsers(varfname,varsSimu1.varphiangle,varsSimu2.varphiangle,...
        'Phi angle comparison',nCells1,nUes1,errTol);

    compareVarsCellsByUsers(varfname,varsSimu1.varthetaangle,varsSimu2.varthetaangle,...
        'Theta angle comparison',nCells1,nUes1,errTol);
    
    compareVarsCellsByUsers(varfname,varsSimu1.varpathloss,varsSimu2.varpathloss,...
        'Pathloss comparison',nCells1,nUes1,errTol);

    compareVarsCellsByUsers(varfname,varsSimu1.varchcoeff,varsSimu2.varchcoeff,...
        'Ch Coeff comparison',nCells1,nUes1,errTol);
    
    compareVarsCellsByUsers(varfname,varsSimu1.varantgain,varsSimu2.varantgain,...
        'Antenna gain comparison',nCells1,nUes1,errTol);
    
    compareVarsCellsByUsers(varfname,varsSimu1.varchplusant,varsSimu2.varchplusant,...
        'Ch Coeff plus Antenna gain comparison',nCells1,nUes1,errTol);
    
    compareVarsCellsByUsers(varfname,varsSimu1.varisloss,varsSimu2.varisloss,...
        'LOS property comparison',nCells1,nUes1,errTol);
    
    compareVarsCellsByUsers(varfname,varsSimu1.varrxpw,varsSimu2.varrxpw,...
        'Received power cell by user',nCells1,nUes1,errTol);
   
    disp('Done!');
end

function compareVarsCellsByUsers(varfname,varData1,varData2,strTitle,...
    nCells,nUes, errTol)

    fileId = fopen(varfname,'a');
    varerror = 0;
    fprintf(fileId,'%s\r\n',strTitle);
    for ii = 1:nCells
        for j = 1:nUes
            vardiff = abs(varData1(ii,j) - varData2(ii,j));
            
            if vardiff > errTol
                fprintf(fileId,'Cell %d - UE %d\r\n',ii,j);
                fprintf(fileId,'%f    %f\r\n',varData1(ii,j),...
                    varData2(ii,j));
                fprintf(fileId,'\r\n');
                varerror = 1;
            end
        end
    end
    if varerror==0
        fprintf(fileId,'OK...\r\n');
    end
    fprintf(fileId,'\r\n');
    fclose(fileId);
end

function compareElemProperty(varfname,varData1, varData2,strTitle,...
    strNameElem,errTol)
    
    numElem1 = length(varData1);
    numElem2 = length(varData2);
    
    fileId = fopen(varfname,'a');
    if numElem1 == numElem2
        varerror = 0;
        fprintf(fileId,'%s\r\n',strTitle);
        for ii = 1:numElem1
            vardiff = abs(varData1(ii) - varData2(ii));
            if vardiff > errTol
                varerror = 1;
                fprintf(fileId,'%s %d: %f    %f\r\n',strNameElem,ii,...
                    varData1(ii),  varData2(ii));
            end
        end
        
        if varerror == 0
            fprintf(fileId,'OK...\r\n');
        end
        
        fprintf(fileId,'\r\n');
    
    end
    fclose(fileId);
end