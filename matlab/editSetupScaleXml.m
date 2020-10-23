function editSetupScaleXml(datDir,modDir,resDir,setDir,name,mass,ht,age,filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Load structure of generic scale file
load([modDir '/structureScale.mat'],'Tree')

% Edit ScaleTool -> subject info
Tree.ScaleTool.ATTRIBUTE.name = name;
Tree.ScaleTool.mass = mass;
Tree.ScaleTool.height = ht;
Tree.ScaleTool.age = age;

% Edit ScaleTool -> GenericModelMaker
Tree.ScaleTool.GenericModelMaker.ATTRIBUTE.name = name;

% Edit ScaleTool -> ModelScaler
Tree.ScaleTool.ModelScaler.ATTRIBUTE.name = name;
Tree.ScaleTool.ModelScaler.marker_file = [datDir '/' filename '.trc'];

% ScaleTool -> MarkerPlacer
Tree.ScaleTool.MarkerPlacer.ATTRIBUTE.name = name;
Tree.ScaleTool.MarkerPlacer.marker_file = [datDir '/' filename '.trc'];
Tree.ScaleTool.MarkerPlacer.output_model_file = ...
    [resDir '/' name 'modelScaled.osim'];
Tree.ScaleTool.MarkerPlacer.output_motion_file = ...
    [resDir '/' filename '.mot'];
Tree.ScaleTool.MarkerPlacer.output_marker_file = ...
    [resDir '/' name 'markersScaled.xml'];

% Set inputs for xml_write
rootName = 'OpenSimDocument';                                               
Pref.StructItem = false;

% Write .xml file
cd(setDir)
xml_write([name 'setupScale.xml'],Tree,rootName,Pref);                           

% Save structure
save([name 'structureScale.mat'],'Tree')

end

