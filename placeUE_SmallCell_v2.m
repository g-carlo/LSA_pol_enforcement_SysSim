function placeUE_SmallCell_v2(posOfCells,ues,squareLength,minDistance,...
    firstIndex,lastIndex,UEheight)

    halfSqLen = squareLength/2;
    
    for n = firstIndex:lastIndex

        xy =  (squareLength*rand(1) - halfSqLen) + ...
            1i*(squareLength*rand(1) - halfSqLen);

        while  ( min(abs(xy - posOfCells)) < minDistance )
             xy =  (squareLength*rand(1) - halfSqLen) + ...
                1i*(squareLength*rand(1) - halfSqLen);
        end

        ues(n).setPosition( real(xy) ,imag(xy) , UEheight);
        
    end    
end        
