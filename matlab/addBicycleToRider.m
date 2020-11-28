function addBicycleToRider(modelPath,modelName)
% -------------------------------------------------------------------------- %
% OpenSim Moco: addBicycleToRider.m                                          %
% -------------------------------------------------------------------------- %
% Copyright (c) 2019 Stanford University and the Authors                     %
%                                                                            %
% Author(s): Ross Wilkinson                                                  %
%                                                                            %
% Licensed under the Apache License, Version 2.0 (the "License"); you may    %
% not use this file except in compliance with the License. You may obtain a  %
% copy of the License at http://www.apache.org/licenses/LICENSE-2.0          %
%                                                                            %
% Unless required by applicable law or agreed to in writing, software        %
% distributed under the License is distributed on an "AS IS" BASIS,          %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   %
% See the License for the specific language governing permissions and        %
% limitations under the License.                                             %
% -------------------------------------------------------------------------- %

% Import OpenSim Libraries into Matlab
% ====================================
import org.opensim.modeling.*

% Initialize the Model
% ====================
osimModel = Model( [modelPath '/' modelName] );

% Change Model name
% =================
newModelName = strrep(modelName,'.osim','WithBicycle.osim');
osimModel.setName(newModelName);

% Define bicycle geometries
% =========================
platformLength = 100; platformHeight = 0.100; platformWidth = 2;
wheelRadius = 0.35; wheelWidth = 0.023;
tubeRadius = 0.012;
crankLength = 0.175;
pedalLength = 0.060; pedalWidth = 0.040; pedalHeight = 0.010;
seatpostLength = 0.380;
saddleWidth = 0.133; saddleLength = 0.180; saddleHeight = 0.02;
handlebarWidth = 0.420; 
toptubeLength = 0.540;
wheelbaseLength = 0.978;
chainstayLength = 0.410;
seatstayLength = 0.450;
bbdropHeight = 0.072;
seattubeAngle = deg2rad(74);
stackHeight = 0.544;
reachLength = 0.384;
downtubeLength = 0.580;
downtubeAngle = 0.85;
frontcenterLength = sqrt(0.579^2 - 0.074^2);
rearcenterLength = sqrt(0.410^2 - 0.074^2);
headtubeAngle = deg2rad(73);
seattubeLength = 0.481;
rearstayLength = sqrt(seattubeLength^2+rearcenterLength^2);
headtubeLength = 0.143;
forkLength = 0.363;
stemLength = 0.100;
stemAngle = deg2rad(6);
bbWidth = 0.073; bbRadius = 0.021;
hubWidth = 0.130;
spindleWidth = 0.150;
steeringTubeLength = 0.183;
chainstayAngle = acos(bbdropHeight/chainstayLength);
hoodLength = 0.100;
hoodAngle = deg2rad(6);
dropLength = 0.075;

% Get a reference to the ground object
% ====================================
ground = osimModel.getGround();

% Define the acceleration of gravity
% ==================================
osimModel.setGravity( Vec3(0, -9.80665, 0) );

% Access existing Model bodies and joints
% =======================================
bodySet = osimModel.getBodySet();
jointSet = osimModel.getJointSet();

% Set reference to Model bodies of interest
% =========================================
leftToes = bodySet().get("toes_l");
rightToes = bodySet().get("toes_r");
leftHand = bodySet().get("hand_l");
rightHand = bodySet().get("hand_r");
pelvis = bodySet().get("pelvis");

% Set reference to Model joints of interest
% =========================================
groundToPelvis = jointSet().get("ground_pelvis");
back = jointSet().get("back");
leftHip = jointSet().get("hip_l");
rightHip = jointSet().get("hip_r");
leftKnee = jointSet().get("walker_knee_l");
rightKnee = jointSet().get("walker_knee_r");
leftAnkle = jointSet().get("ankle_l");
rightAnkle = jointSet().get("ankle_r");
leftMtp = jointSet().get("mtp_l");
rightMtp = jointSet().get("mtp_r");
leftShoulder = jointSet().get("acromial_l");
rightShoulder = jointSet().get("acromial_r");
leftElbow = jointSet().get("elbow_l");
rightElbow = jointSet().get("elbow_r");
leftRadioulnar = jointSet().get("radioulnar_l");
rightRadioulnar = jointSet().get("radioulnar_r");

% Update rider position
% =====================
% pelvis_tilt
pelvis_rx = groundToPelvis.upd_coordinates(0);
pelvis_rx.set_default_value( deg2rad(-20) );
% pelvis_tx
pelvis_tx = groundToPelvis.upd_coordinates(3);
pelvis_tx.set_default_value(0.67);
% pelvis_ty
pelvis_tz = groundToPelvis.upd_coordinates(4);
pelvis_tz.set_default_value(1.2);
% pelvis_tz
pelvis_tz = groundToPelvis.upd_coordinates(5);
pelvis_tz.set_default_value(0);
% lumbar_extension
torso_rz = back.upd_coordinates(0);
torso_rz.set_default_value( deg2rad(-80) );
% hip_flexion_l
leftHip_rz = leftHip.upd_coordinates(0);
leftHip_rz.set_default_value( deg2rad(8) );
% hip_adduction_l
leftHip_ry = leftHip.upd_coordinates(1);
leftHip_ry.set_default_value( deg2rad(-6) );
% hip_flexion_r
rightHip_rz = rightHip.upd_coordinates(0);
rightHip_rz.set_default_value( deg2rad(30) );
% hip_adduction_r
rightHip_ry = rightHip.upd_coordinates(1);
rightHip_ry.set_default_value( deg2rad(-6) );
% knee_angle_l
leftKnee_rz = leftKnee.upd_coordinates(0);
leftKnee_rz.set_default_value( deg2rad(30) );
% knee_angle_r
rightKnee_rz = rightKnee.upd_coordinates(0);
rightKnee_rz.set_default_value( deg2rad(30) );
% ankle_angle_l
leftAnkle_rz = leftAnkle.upd_coordinates(0);
leftAnkle_rz.set_default_value( deg2rad(20) );
% ankle_angle_r
rightAnkle_rz = rightAnkle.upd_coordinates(0);
rightAnkle_rz.set_default_value( deg2rad(10) );
% mtp_l
leftMtp_rz = leftMtp.upd_coordinates(0);
leftMtp_rz.set_default_value( deg2rad(-20) );
% mtp_r
rightMtp_rz = rightMtp.upd_coordinates(0);
rightMtp_rz.set_default_value( deg2rad(-20) );
% arm_flex_l
leftShoulder_rz = leftShoulder.upd_coordinates(0);
leftShoulder_rz.set_default_value( deg2rad(105) );
% arm_add_l
leftShoulder_ry = leftShoulder.upd_coordinates(1);
leftShoulder_ry.set_default_value( deg2rad(-6) );
% arm_flex_r
rightShoulder_rz = rightShoulder.upd_coordinates(0);
rightShoulder_rz.set_default_value( deg2rad(105) );
% arm_add_r
rightShoulder_ry = rightShoulder.upd_coordinates(1);
rightShoulder_ry.set_default_value( deg2rad(-6) );
% % elbow_flex_l
leftElbow_rz = leftElbow.upd_coordinates(0);
leftElbow_rz.set_default_value( deg2rad(5) );
% elbow_flex_r
rightElbow_rz = rightElbow.upd_coordinates(0);
rightElbow_rz.set_default_value( deg2rad(5) );
% radioulnar_l
leftRadioulnar_rz = leftRadioulnar.upd_coordinates(0);
leftRadioulnar_rz.set_default_value( deg2rad(90) );
% radioulnar_r
rightRadioulnar_rz = rightRadioulnar.upd_coordinates(0);
rightRadioulnar_rz.set_default_value( deg2rad(90) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Bodies and Joints %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Body: Platform
% ==============
% platform = Body();
% platform.setName('platform');
% platform.setMass(0.1);
% platform.setInertia( Inertia(0,0,0,0,0,0) );
% % Add geometry to the body
% platformGeometry = Brick( Vec3(platformLength/2,platformHeight/2,platformWidth/2) );
% platformGeometry.setColor( Vec3(1,0,0) );
% platform.attachGeometry(platformGeometry);
% % Add Body to the Model
% osimModel.addBody(platform);

% Joint: Platform {Pin} Ground 
% ============================
% platformToGround = PinJoint('platformToGround',... % Joint Name
%     ground,... % Parent Frame
%     Vec3(0,-1,0),... % Translation in Parent Frame
%     Vec3(0,0,0),... % Orientation in Parent Frame
%     platform,... % Child Frame
%     Vec3(0,0,0),... % Translation in Child Frame
%     Vec3(0,0,0) ); % Orientation in Child Frame
% % RotationZ
% platform_rz = platformToGround.upd_coordinates(0);
% platform_rz.setRange([deg2rad(-100), deg2rad(100)]);
% platform_rz.setName('platform_rz');
% platform_rz.setDefaultValue(deg2rad(0)); % 0 deg. slope
% platform_rz.setDefaultSpeedValue(0);
% platform_rz.setDefaultLocked(true)
% % Add the PlatformToGround Joint to the Model
% osimModel.addJoint(platformToGround);

% Body: Rear Wheel
% ================
rearWheel = Body();
rearWheel.setName('RearWheel');
rearWheel.setMass(0.1);
rearWheel.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry for display
rearWheelGeometry = Cylinder(wheelRadius,wheelWidth/2);
rearWheelGeometry.setColor( Vec3(0.1) );
rearWheel.attachGeometry(rearWheelGeometry);
% Add Body to the Model
osimModel.addBody(rearWheel);

% % Joint: Rear Wheel {Free} Platform 
% % =================================
% rearWheelToPlatform = FreeJoint('rearWheelToPlatform',...
%     platform,... % Parent Frame
%     Vec3(0,0,0),... % Translation in Parent Frame
%     Vec3(pi/2,0,0),... % Orientation in Parent Frame
%     rearWheel,... % Child Frame
%     Vec3(0,0,0),... % Translation in Child Frame
%     Vec3(0,0,0) ); % Orientation in Child Frame
% % RotationX
% rearWheel_rx = rearWheelToPlatform.upd_coordinates(0); 
% rearWheel_rx.setName('rearWheel_rx');
% rearWheel_rx.setDefaultValue(0);
% rearWheel_rx.set_locked(0);
% % RotationY
% rearWheel_ry = rearWheelToPlatform.upd_coordinates(1); 
% rearWheel_ry.setName('rearWheel_ry');
% rearWheel_ry.setRange([-pi,pi]);
% rearWheel_ry.setDefaultValue(0);
% % RotationZ
% rearWheel_rz = rearWheelToPlatform.upd_coordinates(2); 
% rearWheel_rz.setName('rearWheel_rz');
% rearWheel_rz.setDefaultValue(0);
% rearWheel_rz.set_locked(0);
% % TranslationX 
% rearWheel_tx = rearWheelToPlatform.upd_coordinates(3); 
% rearWheel_tx.setRange([0, 50]);
% rearWheel_tx.setName('rearWheel_tx');
% rearWheel_tx.setDefaultValue(0);
% rearWheel_tx.setDefaultSpeedValue(0)
% % TranslationY
% rearWheel_ty = rearWheelToPlatform.upd_coordinates(4); 
% rearWheel_ty.setName('rearWheel_ty');
% rearWheel_ty.setDefaultValue(0);
% rearWheel_ty.set_locked(0);
% % TranslationZ
% rearWheel_tz = rearWheelToPlatform.upd_coordinates(5);
% rearWheel_tz.setRange([-5, 5]);
% rearWheel_tz.setName('rearWheel_tz');
% rearWheel_tz.setDefaultValue(-wheelRadius);
% % Add Joint to model
% osimModel.addJoint(rearWheelToPlatform)

% Body: Right Chainstay
% =====================
rightChainstay = Body();
rightChainstay.setName('rightChainstay');
rightChainstay.setMass(0.1);
rightChainstay.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
rightChainstayGeometry = Cylinder(tubeRadius,chainstayLength/2);
rightChainstayGeometry.setColor( Vec3(1) );
rightChainstay.attachGeometry(rightChainstayGeometry);
% Add Body to the Model
osimModel.addBody(rightChainstay);

% Joint: Right Chainstay {Ball} Rear Wheel 
% =======================================
rightChainstayToRearWheel = BallJoint('rightChainstayToRearWheel',...
    rearWheel,... % Parent Frame
    Vec3(0,bbWidth/2,0),... % Translation in Parent Frame
    Vec3(0,0,pi/2),... % Orientation in Parent Frame
    rightChainstay,... % Child Frame
    Vec3(0,chainstayLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% RotationX
rightChainstay_rx = rightChainstayToRearWheel.upd_coordinates(0); 
rightChainstay_rx.setName('rightChainstay_rx');
rightChainstay_rx.setRange([-pi,pi]);
rightChainstay_rx.setDefaultValue(0);
% RotationY
rightChainstay_ry = rightChainstayToRearWheel.upd_coordinates(1); 
rightChainstay_ry.setName('rightChainstay_ry');
rightChainstay_ry.setDefaultValue(0);
rightChainstay_ry.set_locked(1);
% RotationZ
rightChainstay_rz = rightChainstayToRearWheel.upd_coordinates(2); 
rightChainstay_rz.setName('rightChainstay_rz');
rightChainstay_rz.setDefaultValue(0);
rightChainstay_rz.set_locked(1);
% Add Joint to model
osimModel.addJoint(rightChainstayToRearWheel)

% Body: Left Chainstay
% ====================
leftChainstay = Body();
leftChainstay.setName('leftChainstay');
leftChainstay.setMass(0.1);
leftChainstay.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
leftChainstayGeometry = Cylinder(tubeRadius,chainstayLength/2);
leftChainstayGeometry.setColor( Vec3(1) );
leftChainstay.attachGeometry(leftChainstayGeometry);
% Add Body to the Model
osimModel.addBody(leftChainstay);

% Joint: Left Chainstay {Ball} Rear Wheel 
% ======================================
leftChainstayToRearWheel = BallJoint('leftChainstayToRearWheel',...
    rearWheel,... % Parent Frame
    Vec3(0,-bbWidth/2,0),... % Translation in Parent Frame
    Vec3(0,0,pi/2),... % Orientation in Parent Frame
    leftChainstay,... % Child Frame
    Vec3(0,chainstayLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% RotationX
leftChainstay_rx = leftChainstayToRearWheel.upd_coordinates(0); 
leftChainstay_rx.setName('leftChainstay_rx');
leftChainstay_rx.setRange([-pi,pi]);
leftChainstay_rx.setDefaultValue(0);
% RotationY
leftChainstay_ry = leftChainstayToRearWheel.upd_coordinates(1); 
leftChainstay_ry.setName('leftChainstay_ry');
leftChainstay_ry.setDefaultValue(0);
leftChainstay_ry.set_locked(1);
% RotationZ
leftChainstay_rz = leftChainstayToRearWheel.upd_coordinates(2); 
leftChainstay_rz.setName('leftChainstay_rz');
leftChainstay_rz.setDefaultValue(0);
leftChainstay_rz.set_locked(1);
% Add Joint to model
osimModel.addJoint(leftChainstayToRearWheel)

% Body: Bottom Bracket
% ====================
bottomBracket = Body();
bottomBracket.setName('bottomBracket');
bottomBracket.setMass(0.1);
bottomBracket.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
bottomBracketGeometry = Cylinder(bbRadius,bbWidth/2);
bottomBracketGeometry.setColor( Vec3(1) );
bottomBracket.attachGeometry(bottomBracketGeometry);
% Add Body to the Model
osimModel.addBody(bottomBracket);

% Joint: Bottom Bracket {Weld} Right Chainstay 
% ============================================
bottomBracketToRightChainstay = WeldJoint('bottomBracketToRightChainstay',...
    rightChainstay,... % Parent Frame
    Vec3(0,-chainstayLength/2,0),... % Translation in Parent Frame
    Vec3(0,0,pi/2),... % Orientation in Parent Frame
    bottomBracket,... % Child Frame
    Vec3(0,-bbWidth/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(bottomBracketToRightChainstay)

% Joint: Bottom Bracket {Weld} Left Chainstay 
% ===========================================
bottomBracketToLeftChainstay = WeldJoint('bottomBracketToLeftChainstay',...
    leftChainstay,... % Parent Frame
    Vec3(0,-chainstayLength/2,0),... % Translation in Parent Frame
    Vec3(0,0,pi/2),... % Orientation in Parent Frame
    bottomBracket,... % Child Frame
    Vec3(0,bbWidth/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(bottomBracketToLeftChainstay)

% Body: Seat Tube
% ===============
seatTube = Body();
seatTube.setName('seatTube');
seatTube.setMass(0.1);
seatTube.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
seatTubeGeometry = Cylinder(tubeRadius,seattubeLength/2);
seatTubeGeometry.setColor( Vec3(1) );
seatTube.attachGeometry(seatTubeGeometry);
% Add Body to the Model
osimModel.addBody(seatTube);
% Add marker to the body
%seatTubeMarker1 = Marker("B1",seatTube,Vec3(tubeRadius*2,-0.120,0) );
%osimModel.addMarker(seatTubeMarker1)
% Add marker to the body
%seatTubeMarker2 = Marker("B2",seatTube,Vec3(tubeRadius*4,-0.060,-0.030) );
%osimModel.addMarker(seatTubeMarker2)
% Add marker to the body
%seatTubeMarker3 = Marker("B3",seatTube,Vec3(tubeRadius*4,-0.060,0.030) );
%osimModel.addMarker(seatTubeMarker3)

% Joint: Seat Tube {Weld} Bottom Bracket 
% ======================================
seatTubeToBottomBracket = WeldJoint('seatTubeToBottomBracket',...
    bottomBracket,... % Parent Frame
    Vec3(0,0,0),... % Translation in Parent Frame
    Vec3(pi/2,0,0),... % Orientation in Parent Frame
    seatTube,... % Child Frame
    Vec3(0,seattubeLength/2,0),... % Translation in Child Frame
    Vec3(0,0,(-pi/2)+seattubeAngle) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(seatTubeToBottomBracket)

% Body: Down Tube
% ===============
downTube = Body();
downTube.setName('downTube');
downTube.setMass(0.1);
downTube.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
downTubeGeometry = Cylinder(tubeRadius,downtubeLength/2);
downTubeGeometry.setColor( Vec3(1) );
downTube.attachGeometry(downTubeGeometry);
% Add Body to the Model
osimModel.addBody(downTube);

% Joint: Down Tube {Weld} Bottom Bracket 
% ======================================
downTubeToBottomBracket = WeldJoint('downTubeToBottomBracket',...
    bottomBracket,... % Parent Frame
    Vec3(0,0,0),... % Translation in Parent Frame
    Vec3(pi/2,0,0),... % Orientation in Parent Frame
    downTube,... % Child Frame
    Vec3(0,downtubeLength/2,0),... % Translation in Child Frame
    Vec3(0,0,downtubeAngle) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(downTubeToBottomBracket)

% Body: Top Tube
% ==============
topTube = Body();
topTube.setName('topTube');
topTube.setMass(0.1);
topTube.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
topTubeGeometry = Cylinder(tubeRadius,toptubeLength/2);
topTubeGeometry.setColor( Vec3(1) );
topTube.attachGeometry(topTubeGeometry);
% Add Body to the Model
osimModel.addBody(topTube);

% Joint: Top Tube {Weld} Seat Tube 
% ================================
topTubeToSeatTube = WeldJoint('topTubeToSeatTube',...
    seatTube,... % Parent Frame
    Vec3(0,-seattubeLength/2,0),... % Translation in Parent Frame
    Vec3(0,0,-pi+seattubeAngle),... % Orientation in Parent Frame
    topTube,... % Child Frame
    Vec3(0,toptubeLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(topTubeToSeatTube)

% Body: Head Tube
% ===============
headTube = Body();
headTube.setName('headTube');
headTube.setMass(0.1);
headTube.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
headTubeGeometry = Cylinder(bbRadius,headtubeLength/2);
headTubeGeometry.setColor( Vec3(1) );
headTube.attachGeometry(headTubeGeometry);
% Add Body to the Model
osimModel.addBody(headTube);

% Joint: Head Tube {Weld} Top Tube 
% ================================
headTubeToTopTube = WeldJoint('headTubeToTopTube',...
    topTube,... % Parent Frame
    Vec3(0,-toptubeLength/2,0),... % Translation in Parent Frame
    Vec3(0,0,-headtubeAngle),... % Orientation in Parent Frame
    headTube,... % Child Frame
    Vec3(0,headtubeLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(headTubeToTopTube)

% Body: Left Seatstay
% ===================
leftSeatstay = Body();
leftSeatstay.setName('leftSeatstay');
leftSeatstay.setMass(0.1);
leftSeatstay.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
leftSeatstayGeometry = Cylinder(tubeRadius,seatstayLength/2);
leftSeatstayGeometry.setColor( Vec3(1) );
leftSeatstay.attachGeometry(leftSeatstayGeometry);
% Add Body to the Model
osimModel.addBody(leftSeatstay);

% Joint: Left Seatstay {Weld} Left Chainstay 
% ==========================================
leftSeatStayToChainStay = WeldJoint('leftSeatStayToChainStay',...
    leftChainstay,... % Parent Frame
    Vec3(0,chainstayLength/2,0),... % Translation in Parent Frame
    Vec3(pi/2,0,0),... % Orientation in Parent Frame
    leftSeatstay,... % Child Frame
    Vec3(0,seatstayLength/2,0),... % Translation in Child Frame
    Vec3(pi/4,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(leftSeatStayToChainStay)

% Body: Right Seatstay
% ====================
rightSeatstay = Body();
rightSeatstay.setName('rightSeatstay');
rightSeatstay.setMass(0.1);
rightSeatstay.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
rightSeatstayGeometry = Cylinder(tubeRadius,seatstayLength/2);
rightSeatstayGeometry.setColor( Vec3(1) );
rightSeatstay.attachGeometry(rightSeatstayGeometry);
% Add Body to the Model
osimModel.addBody(rightSeatstay);

% Joint: Right Seatstay {Weld} Right Chainstay
% ============================================
rightSeatStayToChainStay = WeldJoint('rightSeatStayToChainStay',...
    rightChainstay,... % Parent Frame
    Vec3(0,chainstayLength/2,0),... % Translation in Parent Frame
    Vec3(pi/2,0,0),... % Orientation in Parent Frame
    rightSeatstay,... % Child Frame
    Vec3(0,seatstayLength/2,0),... % Translation in Child Frame
    Vec3(pi/4,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(rightSeatStayToChainStay)

% Body: Steering Tube
% ===================
steeringTube = Body();
steeringTube.setName('steeringTube');
steeringTube.setMass(0.1);
steeringTube.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
steeringTubeGeometry = Cylinder(tubeRadius,steeringTubeLength/2);
steeringTubeGeometry.setColor( Vec3(0.1) );
steeringTube.attachGeometry(steeringTubeGeometry);
% Add Body to the Model
osimModel.addBody(steeringTube);

% Joint: Steering Tube {Ball} Head Tube 
% =====================================
steeringTubeToHeadTube = BallJoint('steeringTubeToHeadTube',...
    headTube,... % Parent Frame
    Vec3(0,0,0),... % Translation in Parent Frame
    Vec3(0,0,0),... % Orientation in Parent Frame
    steeringTube,... % Child Frame
    Vec3(0,0,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% RotationX
steeringTube_rx = steeringTubeToHeadTube.upd_coordinates(0); 
steeringTube_rx.setName('steeringTube_rx');
steeringTube_rx.setDefaultValue(0);
steeringTube_rx.set_locked(1);
% RotationY
steeringTube_ry = steeringTubeToHeadTube.upd_coordinates(1); 
steeringTube_ry.setName('steeringTube_ry');
steeringTube_ry.setDefaultValue(0);
steeringTube_ry.setRange([-deg2rad(60),deg2rad(60)]);
% RotationZ
steeringTube_rz = steeringTubeToHeadTube.upd_coordinates(2); 
steeringTube_rz.setName('steeringTube_rz');
steeringTube_rz.setDefaultValue(0);
steeringTube_rz.set_locked(1);
% Add Joint to model
osimModel.addJoint(steeringTubeToHeadTube)

% Body: Stem
% ==========
stem = Body();
stem.setName('stem');
stem.setMass(0.1);
stem.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
stemGeometry = Cylinder(tubeRadius,stemLength/2);
stemGeometry.setColor( Vec3(0.1) );
stem.attachGeometry(stemGeometry);
% Add Body to the Model
osimModel.addBody(stem);

% Joint: Stem {Weld} Steering Tube 
% ================================
stemToSteeringTube = WeldJoint('stemToSteeringTube',...
    steeringTube,... % Parent Frame
    Vec3(0,steeringTubeLength/2,0),... % Translation in Parent Frame
    Vec3(0,0,headtubeAngle+stemAngle),... % Orientation in Parent Frame
    stem,... % Child Frame
    Vec3(0,stemLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(stemToSteeringTube)

% Body: Handlebar
% ===============
handlebar = Body();
handlebar.setName('handlebar');
handlebar.setMass(0.1);
handlebar.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
handlebarGeometry = Cylinder(tubeRadius,handlebarWidth/2);
handlebarGeometry.setColor( Vec3(0.1) );
handlebar.attachGeometry(handlebarGeometry);
% Add Body to the Model
osimModel.addBody(handlebar);

% Joint: Handlebar {Weld} Stem 
% ============================
handlebarToStem = WeldJoint('handlebarToStem',...
    stem,... % Parent Frame
    Vec3(0,-stemLength/2,0),... % Translation in Parent Frame
    Vec3(pi/2,0,0),... % Orientation in Parent Frame
    handlebar,... % Child Frame
    Vec3(0,0,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(handlebarToStem)

% Body: Left Hood
% ===============
leftHood = Body();
leftHood.setName('leftHood');
leftHood.setMass(0.1);
leftHood.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
leftHoodGeometry = Cylinder(tubeRadius,hoodLength/2);
leftHoodGeometry.setColor( Vec3(0.1) );
leftHood.attachGeometry(leftHoodGeometry);
% Add Body to the Model
osimModel.addBody(leftHood);

% Joint: Left Hood {Weld} Handlebar 
% =================================
leftHoodToHandlebar = WeldJoint('leftHoodToHandlebar',...
    handlebar,... % Parent Frame
    Vec3(0,-handlebarWidth/2,0),... % Translation in Parent Frame
    Vec3(pi/2,0,0),... % Orientation in Parent Frame
    leftHood,... % Child Frame
    Vec3(0,-hoodLength/2,0),... % Translation in Child Frame
    Vec3(0,0,-hoodAngle) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(leftHoodToHandlebar)

% Body: Right Hood
% ================
rightHood = Body();
rightHood.setName('rightHood');
rightHood.setMass(0.1);
rightHood.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
rightHoodGeometry = Cylinder(tubeRadius,hoodLength/2);
rightHoodGeometry.setColor( Vec3(0.1) );
rightHood.attachGeometry(rightHoodGeometry);
% Add Body to the Model
osimModel.addBody(rightHood);

% Joint: Right Hood {Weld} Handlebar 
% ==================================
rightHoodToHandlebar = WeldJoint('rightHoodToHandlebar',...
    handlebar,... % Parent Frame
    Vec3(0,handlebarWidth/2,0),... % Translation in Parent Frame
    Vec3(pi/2,0,0),... % Orientation in Parent Frame
    rightHood,... % Child Frame
    Vec3(0,-hoodLength/2,0),... % Translation in Child Frame
    Vec3(0,0,-hoodAngle) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(rightHoodToHandlebar)

% Body: Left Drop
% ===============
leftDrop = Body();
leftDrop.setName('leftDrop');
leftDrop.setMass(0.1);
leftDrop.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
leftDropGeometry = Cylinder(tubeRadius,dropLength/2);
leftDropGeometry.setColor( Vec3(0.1) );
leftDrop.attachGeometry(leftDropGeometry);
% Add Body to the Model
osimModel.addBody(leftDrop);

% Joint: Left Drop {Weld} Left Hood 
% =================================
leftDropToLeftHood = WeldJoint('leftDropToLeftHood',...
    leftHood,... % Parent Frame
    Vec3(0,0,0),... % Translation in Parent Frame
    Vec3(0,0,pi/2),... % Orientation in Parent Frame
    leftDrop,... % Child Frame
    Vec3(0,-dropLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(leftDropToLeftHood)

% Body: Left Drop 2
% =================
leftDrop2 = Body();
leftDrop2.setName('leftDrop2');
leftDrop2.setMass(0.1);
leftDrop2.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
leftDrop2Geometry = Cylinder(tubeRadius,dropLength/2);
leftDrop2Geometry.setColor( Vec3(0.1) );
leftDrop2.attachGeometry(leftDrop2Geometry);
% Add Body to the Model
osimModel.addBody(leftDrop2);

% Joint: Left Drop 2 {Weld} Left Drop 
% ===================================
leftDrop2ToLeftDrop = WeldJoint('leftDrop2ToLeftDrop',...
    leftDrop,... % Parent Frame
    Vec3(0,dropLength/2,0),... % Translation in Parent Frame
    Vec3(0,0,deg2rad(60)),... % Orientation in Parent Frame
    leftDrop2,... % Child Frame
    Vec3(0,-dropLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(leftDrop2ToLeftDrop)

% Body: Right Drop
% ================
rightDrop = Body();
rightDrop.setName('rightDrop');
rightDrop.setMass(0.1);
rightDrop.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
rightDropGeometry = Cylinder(tubeRadius,dropLength/2);
rightDropGeometry.setColor( Vec3(0.1) );
rightDrop.attachGeometry(rightDropGeometry);
% Add Body to the Model
osimModel.addBody(rightDrop);

% Joint: Right Drop {Weld} Right Hood 
% ===================================
rightDropToRightHood = WeldJoint('rightDropToRightHood',...
    rightHood,... % Parent Frame
    Vec3(0,0,0),... % Translation in Parent Frame
    Vec3(0,0,pi/2),... % Orientation in Parent Frame
    rightDrop,... % Child Frame
    Vec3(0,-dropLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(rightDropToRightHood)

% Body: Right Drop 2
% ==================
rightDrop2 = Body();
rightDrop2.setName('rightDrop2');
rightDrop2.setMass(0.1);
rightDrop2.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
rightDrop2Geometry = Cylinder(tubeRadius,dropLength/2);
rightDrop2Geometry.setColor( Vec3(0.1) );
rightDrop2.attachGeometry(rightDrop2Geometry);
% Add Body to the Model
osimModel.addBody(rightDrop2);

% Joint: Right Drop 2 {Weld} Right Drop 
% =====================================
rightDrop2ToRightDrop = WeldJoint('rightDrop2ToRightDrop',...
    rightDrop,... % Parent Frame
    Vec3(0,dropLength/2,0),... % Translation in Parent Frame
    Vec3(0,0,deg2rad(60)),... % Orientation in Parent Frame
    rightDrop2,... % Child Frame
    Vec3(0,-dropLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(rightDrop2ToRightDrop)

% Body: Left Fork
% ===============
leftFork = Body();
leftFork.setName('leftFork');
leftFork.setMass(0.1);
leftFork.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
leftForkGeometry = Cylinder(tubeRadius,forkLength/2);
leftForkGeometry.setColor( Vec3(1) );
leftFork.attachGeometry(leftForkGeometry);
% Add Body to the Model
osimModel.addBody(leftFork);

% Joint: Left Fork {Weld} Steering Tube 
% =====================================
leftForkToSteeringTube = WeldJoint('leftForkToSteeringTube',...
    steeringTube,... % Parent Frame
    Vec3(0,-steeringTubeLength/2,bbWidth/2),... % Translation in Parent Frame
    Vec3(0,0,0),... % Orientation in Parent Frame
    leftFork,... % Child Frame
    Vec3(0,forkLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(leftForkToSteeringTube)

% Body: Right Fork
% ================
rightFork = Body();
rightFork.setName('rightFork');
rightFork.setMass(0.1);
rightFork.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
rightForkGeometry = Cylinder(tubeRadius,forkLength/2);
rightForkGeometry.setColor( Vec3(1) );
rightFork.attachGeometry(rightForkGeometry);
% Add Body to the Model
osimModel.addBody(rightFork);

% Joint: Right Fork {Weld} Steering Tube 
% ======================================
rightForkToSteeringTube = WeldJoint('rightForkToSteeringTube',...
    steeringTube,... % Parent Frame
    Vec3(0,-steeringTubeLength/2,-bbWidth/2),... % Translation in Parent Frame
    Vec3(0,0,0),... % Orientation in Parent Frame
    rightFork,... % Child Frame
    Vec3(0,forkLength/2,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(rightForkToSteeringTube)

% Body: Front Wheel
% =================
frontWheel = Body();
frontWheel.setName('frontWheel');
frontWheel.setMass(0.1);
frontWheel.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
steeringTubeGeometry = Cylinder(wheelRadius,wheelWidth/2);
steeringTubeGeometry.setColor( Vec3(0.1) );
frontWheel.attachGeometry(steeringTubeGeometry);
% Add Body to the Model
osimModel.addBody(frontWheel);

% Joint: Front Wheel {Ball} Left Fork 
% ===================================
frontWheelToLeftFork = BallJoint('frontWheelToLeftFork',...
    leftFork,... % Parent Frame
    Vec3(0,-forkLength/2,-bbWidth/2),... % Translation in Parent Frame
    Vec3(-pi/2,0,0),... % Orientation in Parent Frame
    frontWheel,... % Child Frame
    Vec3(0,0,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% RotationX
frontWheel_rx = frontWheelToLeftFork.upd_coordinates(0); 
frontWheel_rx.setName('frontWheel_rx');
frontWheel_rx.setDefaultValue(0);
frontWheel_rx.set_locked(1);
% RotationY
frontWheel_ry = frontWheelToLeftFork.upd_coordinates(1); 
frontWheel_ry.setName('frontWheel_ry');
frontWheel_ry.setDefaultValue(0);
frontWheel_ry.setRange([-pi,pi]);
% RotationZ
frontWheel_rz = frontWheelToLeftFork.upd_coordinates(2); 
frontWheel_rz.setName('frontWheel_rz');
frontWheel_rz.setDefaultValue(0);
frontWheel_rz.set_locked(1);
% Add Joint to model
osimModel.addJoint(frontWheelToLeftFork)

% Body: Seatpost
% ==============
seatpost = Body();
seatpost.setName('seatpost');
seatpost.setMass(0.1);
seatpost.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
seatpostGeometry = Cylinder(tubeRadius*0.8,seatpostLength/2);
seatpostGeometry.setColor( Vec3(0.1) );
seatpost.attachGeometry(seatpostGeometry);
% Add Body to the Model
osimModel.addBody(seatpost);

% Joint: Seatpost {Sliding} Seat Tube 
% ===================================
seatpostToSeatTube = SliderJoint('seatpostToSeatTube',...
    seatTube,... % Parent Frame
    Vec3(0,-seattubeLength/2,0),... % Translation in Parent Frame
    Vec3(0,0,pi/2),... % Orientation in Parent Frame
    seatpost,... % Child Frame
    Vec3(0,0,0),... % Translation in Child Frame
    Vec3(0,0,pi/2) ); % Orientation in Child Frame
% TranslationX
seatpost_tx = seatpostToSeatTube.upd_coordinates(0); 
seatpost_tx.setName('seatpost_tx');
seatpost_tx.setRange([-seatpostLength,seatpostLength]);
seatpost_tx.setDefaultValue(0);
seatpost_tx.set_locked(1);
% Add Joint to model
osimModel.addJoint(seatpostToSeatTube)

% Body: Saddle
% ============
saddle = Body();
saddle.setName('saddle');
saddle.setMass(0.1);
saddle.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
saddleGeometry = Ellipsoid(saddleLength/2,saddleHeight/2,saddleWidth/2);
saddleGeometry.setColor( Vec3(0.1) );
saddle.attachGeometry(saddleGeometry);
% Add Body to the Model
osimModel.addBody(saddle);

% Joint: Saddle {Planar} Seatpost 
% ===============================
saddleToSeatpost = PlanarJoint('saddleToSeatpost',...
    seatpost,... % Parent Frame
    Vec3(0,-seatpostLength/2,0),... % Translation in Parent Frame
    Vec3(0,pi,(pi/2)-seattubeAngle),... % Orientation in Parent Frame
    saddle,... % Child Frame
    Vec3(0,0,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% RotationZ
saddle_rz = saddleToSeatpost.upd_coordinates(0); 
saddle_rz.setName('saddle_rz');
saddle_rz.setRange([-pi,pi]);
saddle_rz.setDefaultValue(0);
saddle_rz.set_locked(1);
% TranslationX
saddle_tx = saddleToSeatpost.upd_coordinates(1); 
saddle_tx.setName('saddle_tx');
saddle_tx.setRange([-seatpostLength,seatpostLength]);
saddle_tx.setDefaultValue(0);
saddle_tx.set_locked(1);
% TranslationY
saddle_ty = saddleToSeatpost.upd_coordinates(2); 
saddle_ty.setName('saddle_ty');
saddle_ty.setDefaultValue(0);
saddle_ty.set_locked(1);
% Add Joint to model
osimModel.addJoint(saddleToSeatpost)

% Body: Crank Spindle
% ===================
crankSpindle = Body();
crankSpindle.setName('crankSpindle');
crankSpindle.setMass(0.1);
crankSpindle.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
crankSpindleGeometry = Cylinder(tubeRadius,spindleWidth/2);
crankSpindleGeometry.setColor( Vec3(0.1) );
crankSpindle.attachGeometry(crankSpindleGeometry);
% Add Body to the Model
osimModel.addBody(crankSpindle);

% Joint: Crank Spindle {Ball} Bottom Bracket 
% ==========================================
crankSpindleToBottomBracket = BallJoint('crankSpindleToBottomBracket',...
    bottomBracket,... % Parent Frame
    Vec3(0,0,0),... % Translation in Parent Frame
    Vec3(0,0,0),... % Orientation in Parent Frame
    crankSpindle,... % Child Frame
    Vec3(0,0,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% RotationX
crankSpindle_rx = crankSpindleToBottomBracket.upd_coordinates(0); 
crankSpindle_rx.setName('crankSpindle_rx');
crankSpindle_rx.setDefaultValue(0);
crankSpindle_rx.set_locked(1);
% RotationY
crankSpindle_ry = crankSpindleToBottomBracket.upd_coordinates(1); 
crankSpindle_ry.setName('crankSpindle_ry');
crankSpindle_ry.setDefaultValue(0);
crankSpindle_ry.setRange([-pi,pi]);
% RotationZ
crankSpindle_rz = crankSpindleToBottomBracket.upd_coordinates(2); 
crankSpindle_rz.setName('crankSpindle_rz');
crankSpindle_rz.setDefaultValue(0);
crankSpindle_rz.set_locked(1);
% Add Joint to model
osimModel.addJoint(crankSpindleToBottomBracket)

% Body: Left Crank
% ================
leftCrank = Body();
leftCrank.setName('leftCrank');
leftCrank.setMass(0.1);
leftCrank.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
leftCrankGeometry = Brick( Vec3(crankLength/2,0.015,0.005) );
leftCrankGeometry.setColor( Vec3(0.1) );
leftCrank.attachGeometry(leftCrankGeometry);
% Add Body to the Model
osimModel.addBody(leftCrank);
% Add marker to the body
%leftCrankMarker = Marker("LC1",leftCrank,Vec3(-crankLength/2,0,-0.005) );
%osimModel.addMarker(leftCrankMarker)

% Joint: Left Crank {Weld} Crank Spindle 
% ======================================
leftCrankToSpindle = WeldJoint('leftCrankToSpindle',...
    crankSpindle,... % Parent Frame
    Vec3(0,spindleWidth/2,0),... % Translation in Parent Frame
    Vec3(pi/2,0,0),... % Orientation in Parent Frame
    leftCrank,... % Child Frame
    Vec3(-crankLength/2,0,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(leftCrankToSpindle)

% Body: Left Pedal
% ================
leftPedal = Body();
leftPedal.setName('leftPedal');
leftPedal.setMass(0.1);
leftPedal.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
leftPedalGeometry = Brick( Vec3(pedalLength/2,pedalHeight/2,pedalWidth/2) );
leftPedalGeometry.setColor( Vec3(0.1) );
leftPedal.attachGeometry(leftPedalGeometry);
% Add Body to the Model
osimModel.addBody(leftPedal);
% Add marker to the body
%leftPedalMarker = Marker("LP1",leftPedal,Vec3(0,0,-pedalWidth/2) );
%osimModel.addMarker(leftPedalMarker)

% Joint: Left Pedal {Pin} Left Crank
% ==================================
leftPedalToLeftCrank = PinJoint('leftPedalToLeftCrank',...
    leftCrank,... % Parent Frame
    Vec3(crankLength/2,0,-0.005),... % Translation in Parent Frame
    Vec3(0,0,0),... % Orientation in Parent Frame
    leftPedal,... % Child Frame
    Vec3(0,0,pedalWidth/2),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% RotationZ
leftPedal_rz = leftPedalToLeftCrank.upd_coordinates(0); 
leftPedal_rz.setName('leftPedal_rz');
leftPedal_rz.setRange([-pi,pi]);
leftPedal_rz.setDefaultValue(0);
% Add Joint to model
osimModel.addJoint(leftPedalToLeftCrank)

% Body: Right Crank
% =================
rightCrank = Body();
rightCrank.setName('rightCrank');
rightCrank.setMass(0.1);
rightCrank.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
rightCrankGeometry = Brick( Vec3(crankLength/2,0.015,0.005) );
rightCrankGeometry.setColor( Vec3(0.1) );
rightCrank.attachGeometry(rightCrankGeometry);
% Add Body to the Model
osimModel.addBody(rightCrank);
% Add marker to the body
%rightCrankMarker = Marker("RC1",rightCrank,Vec3(crankLength/2,0,0.005) );
%osimModel.addMarker(rightCrankMarker)

% Joint: Right Crank {Weld} Crank Spindle
% =======================================
rightCrankToSpindle = WeldJoint('rightCrankToSpindle',...
    crankSpindle,... % Parent Frame
    Vec3(0,-spindleWidth/2,0),... % Translation in Parent Frame
    Vec3(pi/2,0,0),... % Orientation in Parent Frame
    rightCrank,... % Child Frame
    Vec3(crankLength/2,0,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(rightCrankToSpindle)

% Body: Right Pedal
% =================
rightPedal = Body();
rightPedal.setName('rightPedal');
rightPedal.setMass(0.1);
rightPedal.setInertia( Inertia(0,0,0,0,0,0) );
% Add geometry to the body
rightPedalGeometry = Brick( Vec3(pedalLength/2,pedalHeight/2,pedalWidth/2) );
rightPedalGeometry.setColor( Vec3(0.1) );
rightPedal.attachGeometry(rightPedalGeometry);
% Add Body to the Model
osimModel.addBody(rightPedal);
% Add marker to the body
%rightPedalMarker = Marker("RP1",rightPedal,Vec3(0,0,pedalWidth/2) );
%osimModel.addMarker(rightPedalMarker)

% Joint: Right Pedal {Pin} Right Crank
% ====================================
rightPedalToRightCrank = PinJoint('rightPedalToRightCrank',...
    rightCrank,... % Parent Frame
    Vec3(-crankLength/2,0,0.005),... % Translation in Parent Frame
    Vec3(0,0,0),... % Orientation in Parent Frame
    rightPedal,... % Child Frame
    Vec3(0,0,-pedalWidth/2),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% RotationZ
rightPedal_rz = rightPedalToRightCrank.upd_coordinates(0); 
rightPedal_rz.setName('rightPedal_rz');
rightPedal_rz.setRange([-pi,pi]);
rightPedal_rz.setDefaultValue(0);
% Add Joint to model
osimModel.addJoint(rightPedalToRightCrank)

% Joint: Left Toes {Weld} Left Pedal
% ====================================
leftToesToLeftPedal = WeldJoint('leftToesToLeftPedal',...
    leftPedal,... % Parent Frame
    Vec3(0,0,0),... % Translation in Parent Frame
    Vec3(0,0,pi),... % Orientation in Parent Frame
    leftToes,... % Child Frame
    Vec3(0,0,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(leftToesToLeftPedal)

% Joint: Right Toes {Weld} Right Pedal
% ====================================
rightToesToRightPedal = WeldJoint('rightToesToRightPedal',...
    rightPedal,... % Parent Frame
    Vec3(0,0,0),... % Translation in Parent Frame
    Vec3(0,0,pi),... % Orientation in Parent Frame
    rightToes,... % Child Frame
    Vec3(0,0,0),... % Translation in Child Frame
    Vec3(0,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(rightToesToRightPedal)

% Joint: Left Hand {Weld} Left Drop 2
% ====================================
leftHandToLeftDrop2 = WeldJoint('leftHandToLeftDrop2',...
    leftDrop2,... % Parent Frame
    Vec3(0,0,tubeRadius),... % Translation in Parent Frame
    Vec3(0,pi/2,0),... % Orientation in Parent Frame
    leftHand,... % Child Frame
    Vec3(0,-0.100,0),... % Translation in Child Frame
    Vec3(pi/2,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(leftHandToLeftDrop2)

% Joint: Right Hand {Weld} Right Drop 2
% ====================================
rightHandToRightDrop2 = WeldJoint('rightHandToRightDrop2',...
    rightDrop2,... % Parent Frame
    Vec3(0,0,-tubeRadius),... % Translation in Parent Frame
    Vec3(0,-pi/2,0),... % Orientation in Parent Frame
    rightHand,... % Child Frame
    Vec3(0,-0.100,0),... % Translation in Child Frame
    Vec3(-pi/2,0,0) ); % Orientation in Child Frame
% Add Joint to model
osimModel.addJoint(rightHandToRightDrop2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Construct Contact Geometries %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Contact: Platform {Half Space}
% ==========================
% platformContact = ContactHalfSpace();
% platformContact.setLocation( Vec3(0,platformHeight/2,0) );
% platformContact.setOrientation( Vec3(0,0,-pi/2) );
% platformContact.setFrame(platform)
% platformContact.setName('platformContact');
% osimModel.addContactGeometry(platformContact);

% Contact: Rear Wheel {Sphere}
% ============================
% rearWheelContact = ContactSphere();
% rearWheelContact.setRadius(wheelRadius);
% rearWheelContact.setLocation( Vec3(0) );
% rearWheelContact.setFrame(rearWheel)
% rearWheelContact.setName('rearWheelContact');
% osimModel.addContactGeometry(rearWheelContact);

% Contact: Front Wheel {Sphere}
% =============================
% frontWheelContact = ContactSphere();
% frontWheelContact.setRadius(wheelRadius);
% frontWheelContact.setLocation( Vec3(0) );
% frontWheelContact.setFrame(frontWheel)
% frontWheelContact.setName('frontWheelContact');
% osimModel.addContactGeometry(frontWheelContact);

% Contact: Left Pedal {Sphere}
% ============================
% leftPedalContact = ContactSphere();
% leftPedalContact.setRadius(pedalHeight/2);
% leftPedalContact.setLocation( Vec3(0) );
% leftPedalContact.setFrame(leftPedal);
% leftPedalContact.setName('leftPedalContact');
% osimModel.addContactGeometry(leftPedalContact);

% Contact: Right Pedal {Sphere}
% =============================
% rightPedalContact = ContactSphere();
% rightPedalContact.setRadius(pedalHeight/2);
% rightPedalContact.setLocation( Vec3(0) );
% rightPedalContact.setFrame(rightPedal);
% rightPedalContact.setName('rightPedalContact');
% osimModel.addContactGeometry(rightPedalContact);

% Contact: Left Drop 2 {Sphere}
% =============================
% leftDrop2Contact = ContactSphere();
% leftDrop2Contact.setRadius(tubeRadius);
% leftDrop2Contact.setLocation( Vec3(0) );
% leftDrop2Contact.setFrame(leftDrop2);
% leftDrop2Contact.setName('leftDrop2Contact');
% osimModel.addContactGeometry(leftDrop2Contact);

% Contact: Right Drop 2 {Sphere}
% ==============================
% rightDrop2Contact = ContactSphere();
% rightDrop2Contact.setRadius(tubeRadius);
% rightDrop2Contact.setLocation( Vec3(0) );
% rightDrop2Contact.setFrame(rightDrop2);
% rightDrop2Contact.setName('rightDrop2Contact');
% osimModel.addContactGeometry(rightDrop2Contact);

% Contact: Saddle {Sphere}
% ========================
% saddleContact = ContactSphere();
% saddleContact.setRadius(saddleWidth/4);
% saddleContact.setLocation( Vec3(0,saddleWidth/4,0) );
% saddleContact.setFrame(saddle);
% saddleContact.setName('saddleContact');
% osimModel.addContactGeometry(saddleContact);

% Contact: Pelvis {Sphere}
% ========================
% pelvisContact = ContactSphere();
% pelvisContact.setRadius(0.08);
% pelvisContact.setLocation( Vec3(-0.07,-0.05,0) );
% pelvisContact.setFrame(pelvis);
% pelvisContact.setName('pelvisContact');
% osimModel.addContactGeometry(pelvisContact);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Hunt-Crossley Forces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set properties
% stiffness = 1000000;
% dissipation = 2.0;
% staticFriction = 0.75;
% dynamicFriction = 0.004;
% viscousFriction = 0.47;
% transitionVelocity = 0.001;

% Contact Force: Rear Wheel
% =========================
% rearWheelHuntCrossley = HuntCrossleyForce();
% rearWheelHuntCrossley.setName('rearWheelForce');
% rearWheelHuntCrossley.addGeometry('rearWheelContact');
% rearWheelHuntCrossley.addGeometry('platformContact');
% rearWheelHuntCrossley.setStiffness(stiffness);
% rearWheelHuntCrossley.setDissipation(dissipation);
% rearWheelHuntCrossley.setStaticFriction(staticFriction);
% rearWheelHuntCrossley.setDynamicFriction(dynamicFriction);
% rearWheelHuntCrossley.setViscousFriction(viscousFriction);
% rearWheelHuntCrossley.setTransitionVelocity(transitionVelocity);
% % Add force to model
% osimModel.addForce(rearWheelHuntCrossley);

% Contact Force: Front Wheel
% ==========================
% frontWheelHuntCrossley = HuntCrossleyForce();
% frontWheelHuntCrossley.setName('frontWheelForce');
% frontWheelHuntCrossley.addGeometry('frontWheelContact');
% frontWheelHuntCrossley.addGeometry('platformContact');
% frontWheelHuntCrossley.setStiffness(stiffness);
% frontWheelHuntCrossley.setDissipation(dissipation);
% frontWheelHuntCrossley.setStaticFriction(staticFriction);
% frontWheelHuntCrossley.setDynamicFriction(dynamicFriction);
% frontWheelHuntCrossley.setViscousFriction(viscousFriction);
% frontWheelHuntCrossley.setTransitionVelocity(transitionVelocity);
% % Add force to model
% osimModel.addForce(frontWheelHuntCrossley);

% Contact Force: Saddle
% ==========================
% Set properties
% stiffness = 1000000;
% dissipation = 20.0;
% staticFriction = 1;
% dynamicFriction = 0.75;
% viscousFriction = 1;
% transitionVelocity = 1;
% 
% saddleHuntCrossley = HuntCrossleyForce();
% saddleHuntCrossley.setName('saddleForce');
% saddleHuntCrossley.addGeometry('saddleContact');
% saddleHuntCrossley.addGeometry('pelvisContact');
% saddleHuntCrossley.setStiffness(stiffness);
% saddleHuntCrossley.setDissipation(dissipation);
% saddleHuntCrossley.setStaticFriction(staticFriction);
% saddleHuntCrossley.setDynamicFriction(dynamicFriction);
% saddleHuntCrossley.setViscousFriction(viscousFriction);
% saddleHuntCrossley.setTransitionVelocity(transitionVelocity);
% % Add force to model
% osimModel.addForce(saddleHuntCrossley);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize and Save Model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check model for consistency
osimModel.initSystem();
% Save the model to a file
cd(modelPath)
osimModel.print(newModelName);
disp('Full body model with bicycle printed!');

end
