%%% SON SIMULATOR: Final data analysis function class
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010
%   Last Modification: March, 2015
%

function dataAnalysisFnc_v2(objCarrier, varSimData,trafficGenType)

    %%% UE SINR and spectral efficiency 
    % size: 1 X total number of UEs, average spectral efficiency per UE
    ueSpectralEfficiency = objCarrier.ueCumSpectralEfficiency ./ ...
        objCarrier.ueSpectralEfficiency_cnt;

    ueSpectralEfficiency(objCarrier.ueSpectralEfficiency_cnt <= 0) = 0;
    UE_SINR_dB = 10*log10(2 .^ ueSpectralEfficiency-1);

    varSimData.UEspectralEfficiecy = ueSpectralEfficiency;

    varSimData.UE_SINR_dB = UE_SINR_dB;

    %%% network load
    if trafficGenType > 1       
        varSimData.resourceUsagePercentage = objCarrier.cumResourceUsage / objCarrier.cumResourceAvailability *100;
    end
end






