function [RFxyz_r,RFxyz_l] = rotateCrankForces(rTan,rRad,lTan,lRad,angleR,angleL)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
 
% Convert crank angle from CW to CCW
angleR_ccw = 2 * pi - angleR;
angleL_ccw = 2 * pi - angleL;

% Convert zero from TDC to global zero i.e CW pi/2 radians (90 deg)
thetaR = angleR_ccw + (pi / 2);
thetaL = angleL_ccw + (pi / 2);

% Pre-allocate output variable size
[RFxyz_r, RFxyz_l] = deal(NaN(3,numel(rTan)));

% Create a 2-D rotation matrix for each sample, which will transform 
% forces from the crank reference frame to the bicycle reference frame.
for iData = 1:numel(rTan)
    % 2-D rotation matrix = [cos() -sin(); sin() cos()]
    rotationMatrix2dRight = ...
        [cos(thetaR(iData)) -sin(thetaR(iData)); sin(thetaR(iData)) cos(thetaR(iData))];
    rotationMatrix2dLeft = ...
        [cos(thetaL(iData)) -sin(thetaL(iData)); sin(thetaL(iData)) cos(thetaL(iData))];
    
    % Set positive x-axis in crank reference frame as negative radial force.
    forceXaxisCrankRight = -rRad(iData); 
    forceXaxisCrankLeft = -lRad(iData);
    
    % Set positive y-axis in crank coordinate system as negative tangential force.
    forceYaxisCrankRight = -rTan(iData);
    forceYaxisCrankLeft = -lTan(iData); 
    
    % Multiply by 2D rotation matrix to transform to bicycle coordinate system.
    forceXyzBicycleRight = rotationMatrix2dRight * ...
        [forceXaxisCrankRight;forceYaxisCrankRight];
    forceXyzBicycleLeft = rotationMatrix2dLeft * ...
        [forceXaxisCrankLeft;forceYaxisCrankLeft];

    % Set z coordinate to zero. No forces in z-axis.
    forceXyzBicycleRight(3) = 0;
    forceXyzBicycleLeft(3) = 0;
    
    % Set input forces for OpenSim as reaction force.
    RFxyz_r(:,iData) = -forceXyzBicycleRight;
    RFxyz_l(:,iData) = -forceXyzBicycleLeft;
    
end

end

