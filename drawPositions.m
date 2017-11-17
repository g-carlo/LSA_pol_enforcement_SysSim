function drawPositions(objArray,numOfObj,cellTypeIndex,...
    firstIndex,lastIndex,varIntersiteDist,numOfCellsPerSite,varFname,...
    printCellPos)

    varObjPos = [objArray(firstIndex:lastIndex).pos];
    varComplexPos = [varObjPos(1:numOfObj).x] + 1i*[varObjPos(1:numOfObj).y];
    printPositionsFromArrayToFile(varFname,varComplexPos, firstIndex:lastIndex);
    
   
     if printCellPos
         drawRegularHexagon(varIntersiteDist,varFname,numOfObj,cellTypeIndex,...
             numOfCellsPerSite);
     else
         drawUsers(varFname,numOfObj);
     end
end