%%% SON SIMULATOR: antenna MkI class
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010
%
%   Description: this class implements the antenna class with the directive antenna
%   pattern given by 3GPP TR 36.814

function gainMtx = directiveAntExternalFnc(phiMtx,thetaMtx,boresightMtx)

global config;
global cls;

%%% Parameters for antenna gain

gainMinTheta    = 20;
gainMinPhi      = 25;
theta3db        = 10;
phi3db          = 70;
gainMinThetaMtx = gainMinTheta*ones(size(thetaMtx));
gainMinPhiMtx   = gainMinPhi*ones(size(phiMtx));
elet_downtilt   = 15;

phi = angle(exp(1i*(phiMtx-boresightMtx)))/pi*180;
theta = thetaMtx/pi*180;          
gainH = -min(12.*(( phi )./phi3db ).^2,gainMinPhiMtx);
% gainV = zeros(size(gainH));           % for calibration
gainV = -min(12.*(( theta - elet_downtilt )./theta3db ).^2,gainMinThetaMtx);
gainMtx = 10.^( (-min(-(gainV+gainH),gainMinPhiMtx) / 10 ));
