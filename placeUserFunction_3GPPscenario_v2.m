 
function placeUserFunction_3GPPscenario_v2(ues,deploy,UEheight,cls,cellTypeArray)
    
    global macroInd;
    global picoInd;
    
    varIntersiteDist = deploy.intersiteDistance;
    UEdistribConfig = deploy.UEdistributionConfiguration;
        
    l = varIntersiteDist/2;
    b = l/(cos(pi/6));        

    initC = [b,b,2*l,l];
    
    switch UEdistribConfig    
        case 1 % configuration 1 of users distribution (see Req. Doc 3GPP TR 36.814)
            
            setUEposition(ues,UEheight,macroInd,cellTypeArray,initC,cls,...
                cellTypeArray(macroInd).firstUserId,macroInd,deploy);      

        case 2 % configuration 2 of users distribution (see Req. Doc 3GPP TR 36.814)
            
            rPico = cellTypeArray(picoInd).radius;

            setUEposition(ues,UEheight,[macroInd,picoInd],cellTypeArray,...
                initC,cls,cellTypeArray(macroInd).firstUserId,macroInd,deploy);

            initC  = [2*rPico,rPico,2*rPico,rPico];
            
            setUEposition(ues,UEheight,[macroInd,picoInd],cellTypeArray,initC,...
                cls,cellTypeArray(picoInd).firstUserId,picoInd,deploy);       
    end
end