function drawUsers( varFname,numOfUsers )
    
    varStyle = strcat('k','o');
    varIndex = 1:numOfUsers-1;
    
    varUserPos = readCenterCoordinates(varFname,numOfUsers);

    plot(real(varUserPos(varIndex)),imag(varUserPos(varIndex)),varStyle)    
end

