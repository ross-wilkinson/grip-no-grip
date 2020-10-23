function [Pxyz_r, Pxyz_l] = calculateCrankForceOrigin(markerDataBicycle,markerLabels)
%UNTITLED Use foot markers to estimate origin of pedal force
%   Detailed explanation goes here

% Set marker variables
rtoe = markerDataBicycle(:,:,markerLabels == 'rtoe');
rmt5 = markerDataBicycle(:,:,markerLabels == 'rmt5');
rcal = markerDataBicycle(:,:,markerLabels == 'rcal');
ltoe = markerDataBicycle(:,:,markerLabels == 'ltoe');
lmt5 = markerDataBicycle(:,:,markerLabels == 'lmt5');
lcal = markerDataBicycle(:,:,markerLabels == 'lcal');

% Pre-allocate output variable size
[Pxyz_r, Pxyz_l] = deal(NaN(3,length(rtoe)));
    
for iData = 1:length(rtoe)
    % Get distance from calc to toe markers on each foot.
    dispCalcToToeRight = rtoe(:,iData) - rcal(:,iData);
    dispCalcToToeLeft = ltoe(:,iData) - lcal(:,iData);
    
    % Set virtual marker position as toe marker minus 1/3 of the distance to calc 
    % marker.
    virtualMarkerRight = rtoe(:,iData) - dispCalcToToeRight * (1/3);
    virtualMarkerLeft = ltoe(:,iData) - dispCalcToToeLeft * (1/3);
    
    % Use mt5 for Z coordinate of virtual marker.
    virtualMarkerRight(3) = rmt5(3,iData);
    virtualMarkerLeft(3) = lmt5(3,iData);
      
    % Swap axes from UQ motion capture lab to match OpenSim.
    % Xmodel = Xlab. Ymodel = Zlab. Zmodel = -Ylab.
    Pxyz_r(1,iData) = virtualMarkerRight(1);
    Pxyz_l(1,iData) = virtualMarkerLeft(1);
    Pxyz_r(2,iData) = virtualMarkerRight(3);
    Pxyz_l(2,iData) = virtualMarkerLeft(3);
    Pxyz_r(3,iData) = -virtualMarkerRight(2);
    Pxyz_l(3,iData) = -virtualMarkerLeft(2); 
end

% Convert marker coordinate units from millimeters (mm) to SI Units (m)
mmToM = 1000;
Pxyz_r = Pxyz_r / mmToM;
Pxyz_l = Pxyz_l / mmToM;

end

