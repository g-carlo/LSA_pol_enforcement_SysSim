function printPositionsFromArrayToFile(varfname,varPosArray, varIdArray)
    varFile = fopen(varfname,'w');
    varLen = length(varPosArray);
    
    for varInd=1:varLen
        fprintf(varFile,'Object ID: %d        %f %f\r\n',...
                varIdArray(varInd),real(varPosArray(varInd)),imag(varPosArray(varInd)));
    end
    fclose(varFile);
end