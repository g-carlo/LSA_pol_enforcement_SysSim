function cellsPerCarrier = readCoverageMap(numOfCarriers,totNumOfCells)
    
    varFname = 'coverageMap.txt';
    fileId = fopen(varFname,'r');
    
    cellsPerCarrier = zeros(numOfCarriers,totNumOfCells);
    carrierNum = 1;
    
    while ~feof(fileId) && carrierNum <= numOfCarriers
        varline = fgets(fileId);
        varlen = length(varline);
        j = 1;
        for ii=1:varlen
            varchar = varline(ii);
            if '0'<=varchar && varchar<='9'
                varnum(j)=varchar;
                j = j + 1;
            else
                if ~isempty(varnum)>0
                    varnum2 = str2num(varnum);
                    if varnum2<=totNumOfCells
                        cellsPerCarrier(carrierNum,varnum2) = 1;
                    end
                
                    j = 1;
                    varnum = '';
                end
            end
        end
        carrierNum = carrierNum + 1;
    end
    fclose(fileId);
end