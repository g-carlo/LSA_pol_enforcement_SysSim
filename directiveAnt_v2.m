function varAntGain = directiveAnt_v2(strAnt, cellOrient, phiAngle,...
                thetaAngle, elecDowntiltCarrier)

    gainMinTheta = strAnt.verticalMinGain;
    theta3db = strAnt.verticalTheta3dB;

    gainMinPhi = strAnt.horizontalMinGain;            
    phi3db = strAnt.horizontalPhi3dB;

    phi = angle(exp(1i .* (phiAngle - cellOrient)))./pi .* 180;

    theta = thetaAngle ./ pi .* 180;          
    gainH = -min(12 .* (( phi )./ phi3db ) .^ 2, gainMinPhi);
    gainV = -min(12 .* (( theta - elecDowntiltCarrier ) ./ theta3db ) .^ 2,...
        gainMinTheta);
    varAntGain = 10 .^ ( (-min(-(gainV+gainH),gainMinPhi) ./ 10 ));

end