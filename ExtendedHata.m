function [pathloss, FSpathloss] = ExtendedHata(varDist, freq, height_rx, height_tx, clutter)
%% function pathloss = ExtendedHata(distance, frequency, height_rx, height_tx, clutter)
%
% Calculates the pathloss between points based on the Extended Hata model
% 
% Martin Fenton 05/04/2011
%
% Inputs:
%
%	clutter   = Clutter category (of the receiver): 7 = Dense Urban;
%                                                   6 = Urban;
%                                                   5 = Open in urban;
%                                                   4 = Suburban;
%                                                   3 = Forest;
%                                                   2 = Rural; and
%                                                   1 = Water
%	height_tx = Transmitter height (m)
%	height_rx = Receiver height (m)
%	distance  = Matrix of horizontal distances (km) between transmitters 
%               and the receiver . 
%	frequency = frequency (MHz) of the transmissions
%
% Outputs:
%
%	pathloss  = Vector of pathlosses(dB)

%% Calculate pathloss
    
    [nCells,nUes] = size( varDist );
    pathloss = zeros(nCells,nUes);
    
    vartmp1 = 32.44 + 20 * log10(freq);
    vartmp2 = ((height_tx - height_rx) .^ 2) ./ 1000000;
        
    FSpathloss =  vartmp1 + 10 .* log10(  (varDist .^ 2) + vartmp2  );
    
    EHpathloss = ExtendedHata100(varDist, freq, nCells, nUes, height_rx, ...
        height_tx, clutter);
    
    varDistTmp = 0.04 .* ones(nCells,nUes);
    
    L40m = vartmp1 + 10 * log10((varDistTmp .^ 2) + vartmp2);
    
    varDistTmp = 0.1 .* ones(nCells,nUes);
    
    L100m = ExtendedHata100(varDistTmp, freq, nCells, nUes,height_rx, ...
        height_tx, clutter);
    
    vartmp = varDist<=0.04;
    pathloss(vartmp) = FSpathloss(vartmp);
    
    vartmp = varDist>0.1;
    pathloss(vartmp) = EHpathloss(vartmp);
    
    varA = varDist > 0.04;
    varB = varDist <= 0.1;
    
    vartmp =  logical(varA .* varB);
    
    vartmp2 = log10(0.1) - log10(0.04);
    vartmp3 = zeros(nCells,nUes);
    vartmp3(vartmp) = log10( varDist(vartmp) ) - log10(0.04);
    
    pathloss(vartmp) = L40m(vartmp) + ((L100m(vartmp) - L40m(vartmp)) .* (vartmp3(vartmp) ./ vartmp2));
end

function pathloss = ExtendedHata100(varDist, freq,nCells,nUes, height_rx, height_tx, clutter)
%% function pathloss = ExtendedHata100(distance, frequency, height_rx, height_tx, clutter)
% 
% Sub function called by the main ExtendedHata function to calculate the   
% pathloss between points based on the Extended Hata model(for distances
% greater than 100 metres)
%
% Martin Fenton 05/04/2011
%
% Inputs: - same as for the main ExtendedHata function above
%
% Outputs: pathloss = Vector of pathlosses (dB)

%% Calculate pathloss

    a = ((1.1 * log10(freq) - 0.7) * min(10, height_rx)) - (1.56 * log10(freq) - 0.8) + max(0, (20 * log10(height_rx / 10)));
    b = min(0, (20 * log10(height_tx / 30)));

    Alpha = ones(nCells,nUes);
    index = Alpha > 20;
    Alpha(index) = 1 + (0.14 + 0.000187 * freq + 0.00107 * height_tx) .* ((log10(varDist(index) ./ 20)) .^ 0.8);

    Beta = ((44.9 - (6.55 * log10(max(30, height_tx)))) .* (log10(varDist) .^ Alpha)) - (13.82 * log10(max(30, height_tx)));

    if freq <= 150
        Lurban = 69.6 + (26.2 .* log10(150)) - (20 .* log10(150 / freq)) ...
            + Beta - a - b;
    elseif freq <= 1500
        Lurban = 69.6 + (26.2 .* log10(freq)) + Beta - a - b;
    elseif freq <= 2000
        Lurban = 46.3 + (33.9 .* log10(freq)) + Beta - a - b;
    elseif freq <= 3000
        Lurban = 46.3 + (33.9 .* log10(freq)) + (10 .* log10(freq / 2000))...
            + Beta - a - b;
    end

    Gamma = min(max(150, freq), 2000);

    switch clutter
        case 7      % DenseUrban
            pathloss = Lurban + 3;
        case 6      % Urban
            pathloss = Lurban;
        case 4      % Suburban
            pathloss = Lurban - 5.4 - (2 .* ((log10(Gamma ./ 28)) .^ 2));
        otherwise   % Open/Rural
            pathloss = Lurban - 40.94 - (4.78 .* ((log10(Gamma)) .^ 2)) + (18.33 .* log10(Gamma));
    end
end
