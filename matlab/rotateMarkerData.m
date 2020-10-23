function [markerDataBicycle,markerLabels] = rotateMarkerData(path,filename)
%ROTATEMARKERDATA Rotates marker data from Qualisys reference frame to
%ergometer reference frame using markers placed perpendicular to the 
%heading angle of the ergometer.

% Load file
data = load([path '/' filename '.mat']);

% Set marker data, marker labels, and analog labels to new variables.
markerData = data.(filename).Trajectories.Labeled.Data(:,1:3,:);
markerLabels = categorical(data.(filename).Trajectories.Labeled.Labels);

% Permute marker data for writing .trc file. One marker per page.
markerData = permute(markerData,[2 3 1]);
        
% Create bicycle reference frame using markers on base of ergometer
B2 = markerData(:,1,markerLabels == 'LBIKE');
B3 = markerData(:,1,markerLabels == 'RBIKE');

% Use xyz distance between markers on rear legs of LODE to estimate heading 
% angle of ergometer.
dist = B2' - B3';

% Create new y-axis in bicycle coordinate system. Divide by Euclidean length 
% of the vector.
yAxisBicycle = dist / norm(dist);

% Use global z-axis to calculate x-axis in bicycle coordinate system.
zAxisBicycle = [0 0 1];

% Take cross product of y-axis unit vector (bicycle) and z-axis unit vector 
% (global) to find new x-axis unit vector for bicycle coordinate system.
xAxisBicycle = cross(yAxisBicycle,zAxisBicycle); 

% Now use x and y-axis in bicycle reference frame to calculate new z-axis.
zAxisBicycle = cross(xAxisBicycle,yAxisBicycle);
      
% Set bicycle unit vectors in rotation matrix (global -> bicycle).
rotationMatrix3d = [xAxisBicycle;yAxisBicycle;zAxisBicycle];

% Rotate marker data from global to bicycle reference frame
markerDataBicycle = zeros(size(markerData));
for i = 1:length(markerLabels)
    markerDataBicycle(:,:,i) = rotationMatrix3d * markerData(:,:,i);
end

end

