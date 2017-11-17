function drawRegularHexagon(varIntersiteDist,varFname,numOfCells,...
    varCellType,numCellsPerSite)
    
    global macroInd;
    global picoInd;
    global smallInd;
    
    varIndex = 1:numOfCells;
    
    varHalfIntersiteDist = varIntersiteDist/2;

    varSideLen = varHalfIntersiteDist/cos(pi/6);
    
    varHexCenterArray = readCenterCoordinates(varFname,numOfCells);
    
    varClr = cell(3,1);
    varMarker = cell(3,1);
    varColourIndex = ones(1,3);
    varMarkerIndex = ones(1,3);
    
    varClr{macroInd} = ['m' 'c' 'r' 'g' 'b'];
    varClr{picoInd} = ['m' 'c' 'r'];
    varClr{smallInd} = ['g' 'b'];
    
    varMarker{macroInd} = [];
    varMarker{picoInd} = ['+' 'x'];
    varMarker{smallInd} = '*';    
    
    switch varCellType
        case macroInd 
         
            %number of sites. Each hexagon represents one site.
            varNumOfHex = numOfCells/numCellsPerSite;
            cellCenterIndex = -numCellsPerSite+1;

            for varHex=1:varNumOfHex

                varX(1) = varSideLen/2;        
                varY(1) = varHalfIntersiteDist;

                varXpoint1 = varX(1);
                varYpoint1 = varY(1);

                varX(2) = 2*varXpoint1;
                varY(2) = 0;

                varX(3) = varXpoint1;
                varY(3) = -1*varYpoint1;

                varX(4) = -1*varXpoint1;        
                varY(4) = -1*varYpoint1;

                varX(5) = -2*varXpoint1;        
                varY(5) = 0;

                varX(6) = -1*varXpoint1;        
                varY(6) = varYpoint1;

                varX(7) = varXpoint1;        
                varY(7) = varYpoint1;

                cellCenterIndex = cellCenterIndex + numCellsPerSite;
                varCtrSiteX = real(varHexCenterArray(cellCenterIndex));
                varCtrSiteY = imag(varHexCenterArray(cellCenterIndex));

                varX = varX + varCtrSiteX;
                varY = varY + varCtrSiteY;

                plot(varX,varY);
                str1 =  num2str(varHex);
                text(varCtrSiteX,varCtrSiteY,str1)


                varLine1X = [0+varCtrSiteX 0+varCtrSiteX];
                varLine1Y = [0+varCtrSiteY varHalfIntersiteDist+varCtrSiteY];

                plot(varLine1X,varLine1Y,'k--');

                varLine1X = [0+varCtrSiteX varHalfIntersiteDist*cos(pi/6)+varCtrSiteX];
                varLine1Y = [0+varCtrSiteY -1*varHalfIntersiteDist*sin(pi/6)+varCtrSiteY];

                plot(varLine1X,varLine1Y,'k--');

                varLine1X = [0+varCtrSiteX -1*varHalfIntersiteDist*cos(pi/6)+varCtrSiteX];
                varLine1Y = [0+varCtrSiteY -1*varHalfIntersiteDist*sin(pi/6)+varCtrSiteY];

                plot(varLine1X,varLine1Y,'k--');
            end
        case {picoInd,smallInd}
            
            clrIndex = varColourIndex(varCellType);
            markIndex = varMarkerIndex(varCellType);
            
            varStyle = strcat(varClr{varCellType}(clrIndex),varMarker{varCellType}(markIndex));
            plot(real(varHexCenterArray(varIndex)),imag(varHexCenterArray(varIndex)),varStyle)
            
            if clrIndex+1>length(varClr{varCellType})
                markIndex = rem(markIndex,length(varMarker{varCellType}))+1;                
            end
            
            clrIndex = rem(clrIndex,length(varClr{varCellType}))+1;
            
            varColourIndex(varCellType) = clrIndex;
            varMarkerIndex(varCellType) = markIndex;
    end   
end

