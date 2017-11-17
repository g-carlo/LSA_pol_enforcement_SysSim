function posIncumbents = compute_incunbent_positions(n_incumbents,areaLimVector)

    X_len = areaLimVector(2) - areaLimVector(1); 
    Y_len = areaLimVector(4) - areaLimVector(3); 

    %%%% This is the working version tested with long simulation
%     Define the location of the IU
% %     Unit = 500;
% %     N_IUs = 5; % new
% %     Center_IU = 10;
% %     IU_Loc = [Center_IU zeros(1,4)];
% %     Rad = [Unit zeros(1,3)];
% %     Cen = [Center_IU zeros(1,3)];
% %     
% %     sign_vec = [1 1; 1 -1; -1 1; -1 -1];
% %     
% %     for iu=2:N_IUs
% %         IU_Loc(iu) = Cen(iu-1) + [(1+rand)*Rad(iu-1) 1i*(1+rand)*Rad(iu-1)]*sign_vec(iu-1,:)' ;
% %         Cen(iu) = sum(IU_Loc)/iu;
% %         Rad(iu) =  abs(IU_Loc(iu)-Cen(iu)); 
% %     end
% %     
% %     posIncumbents = IU_Loc;
    %%%% END of "This is the working version tested with long simulation"
    
    %%%% This is the working version tested with long simulation
    % Define the location of the IU
    Radius = 1000;
    inter_Inc_distance = 200;
    
    N_IUs = 5; % new
    Center_IU = 800;
    IU_Loc = [Center_IU];
    for n = 2:N_IUs
        newPos = Center_IU + Radius + 2*Radius* ((rand(1)-0.5)+1i*(rand(1)-0.5));
    
        while abs( newPos - Center_IU )>Radius || any(abs( newPos - IU_Loc ) < inter_Inc_distance)
        
            newPos = Center_IU + Radius + 2*Radius* ((rand(1)-0.5)+1i*(rand(1)-0.5));
        
        end
        IU_Loc = [IU_Loc newPos];
    end
    
    posIncumbents = IU_Loc;
    %%%% END of "This is the working version tested with long simulation"  


%     posIncumbents = 10;
%     for n = 2:n_incumbents
%         candidatePos = (rand(1)*X_len  - X_len/2)  + 1i*( rand(1)*Y_len - Y_len/2) ;
%         dist_to_other_incumbents = abs(candidatePos-posIncumbents);
%         isTooClose = dist_to_other_incumbents <= 4*50;
%         while any(isTooClose)
%             candidatePos = (rand(1)*X_len  - X_len/2)  + 1i*( rand(1)*Y_len - Y_len/2) ;
%             dist_to_other_incumbents = abs(candidatePos-posIncumbents);
%             isTooClose = dist_to_other_incumbents <= 4*50;
%         end
%         posIncumbents = [posIncumbents candidatePos];
%     end
   

end