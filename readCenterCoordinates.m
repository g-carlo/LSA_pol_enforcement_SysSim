function varHexCenterArray = readCenterCoordinates(varFname,varNumOfCells)
    
    varFile = fopen(varFname,'r');
    varHexCenterArray = (0 + 1i*0).*ones(1,varNumOfCells);
       
    for varCell=1:varNumOfCells
        varAux = fscanf(varFile,'Object ID: %d',1);
        varReal = fscanf(varFile,'%f',1);
        varImag = fscanf(varFile,'%f',1);
        varAux = fscanf(varFile,'%c',1);
        varAux = fscanf(varFile,'%c',1);
        varComplexVal = complex(varReal,varImag);
        varHexCenterArray(varCell)=varComplexVal;
    end        
    fclose(varFile);
end