classdef carrierCellType_Cls<handle
    properties        
        typeIndex;
        carrierIndex;
        antennaGain;
        antennaType;
        maxTxPw;
        schedulingType;  
        pathLossModel;
        coeff;        
    end
    
    methods
        function obj = carrierCellType_Cls()
            obj.coeff = zeros(2,2);
        end
        
        function determineCoeffPathLoss(obj,varFreq,macroInd,picoInd,...
                smallInd)
            
            switch obj.pathLossModel
                case 1
                    switch obj.typeIndex
                        case macroInd
                             obj.coeff(1,1) = 128.1;
                             obj.coeff(1,2) = 37.6;
                        case picoInd
                            obj.coeff(1,1) = 140.7;
                            obj.coeff(1,2) = 36.7;
                        otherwise
                            varstr = ['Path loss model 1 is implemented for '...
                                      'macro and pico cells only'];
                            error(varstr);
                    end                    
                case 2
                    switch obj.typeIndex
                        case macroInd
                             % coefficients for channel with LOS property
                            obj.coeff(1,1) = 103.4;
                            obj.coeff(1,2) = 24.2;

                            % coefficients for channels with NLOS property
                            obj.coeff(2,1) = 131.1;
                            obj.coeff(2,2) = 42.8;
                        case picoInd
                            obj.coeff(1,1) = 103.8;
                            obj.coeff(1,2) = 20.9;

                            obj.coeff(2,1) = 145.4;
                            obj.coeff(2,2) = 37.5;
                        otherwise
                            varstr = ['Path loss model 2 is implemented for '...
                                      'macro and pico cells only'];
                            error(varstr);
                    end   
                    
                case 3
                    if obj.typeIndex==smallInd
                        obj.coeff(1,1) = 20*log10(4*pi*varFreq/3e8)+60;
                        obj.coeff(1,2) = 20;
                    else
                        varstr = ['Path loss model 3 is implemented for '...
                                  'small cells only'];
                        error(varstr);
                    end
                case 4
                    if obj.typeIndex == smallInd
                        obj.coeff(1,1) = 140.7;
                        obj.coeff(1,2) = 36.7;
                    else
                        varstr = ['Path loss model 4 is implemented for '...
                                  'small cells only'];
                        error(varstr);                    
                    end
            end
        end      
        
    end
    
end

