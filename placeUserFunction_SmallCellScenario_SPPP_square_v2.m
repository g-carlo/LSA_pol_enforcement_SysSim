function placeUserFunction_SmallCellScenario_SPPP_square_v2(squareLen,...
    UEdistConfig,numOfUsers,minDistUEsmallCell,UEheight)

    posOfSmallCells = [cls(firstIndex:lastIndex).complexPosition];


    if UEdistConfig == 1
       for n = 1:numOfUsers

            xy =  (squareLen*rand(1)-squareLen/2) + 1i*(squareLen*rand(1)-squareLen/2);

            while  ( min(abs(xy - posOfSmallCells)) < minDistUEsmallCell)

                xy =  (squareLength*rand(1)-squareLength/2) + 1i*(squareLength*rand(1)-squareLength/2);

            end

            %%% THIS IS THE OUTPUT OF THE FUNCTION - BEGIN
            ues(n).setPosition( real(xy) ,imag(xy) , UEheight);
            %%% THIS IS THE OUTPUT OF THE FUNCTION - END


        end 
    end            
end