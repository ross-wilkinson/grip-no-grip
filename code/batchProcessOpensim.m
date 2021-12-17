%% Batch process OpenSim
% Run after creation of subject-specific scaled model

%% Initialize workspace
clear; clc; close all

%% Initialize OpenSim API
import org.opensim.modeling.*

%% Set project folder paths
expDir = '/Users/rosswilkinson/Google Drive/projects/grip-no-grip';
datDir = [expDir '/data'];
codDir = [expDir '/matlab'];
modDir = [expDir '/model'];
resDir = [expDir '/results'];
setDir = [expDir '/setup'];

%% Specify subject to analyze
prompt = 'Which subject would you like to process? s01, s02, etc.: ';
subName = input(prompt,'s');

%% Get subject filenames
cd(datDir)
subFileList = dir([subName '*.c3d']);
nFiles = size(subFileList,1);

%% Batch process all trials
for iFiles = 2:nFiles
    trialFileName = strrep(subFileList(iFiles).name,'.c3d','');
        
    % Setup and run the OpenSim IK Tool for this trial
    ikSetupFileName = [setDir '/' trialFileName 'setupInverseKinematics.xml'];
    ikTool = InverseKinematicsTool(ikSetupFileName);
    ikTool.run();
    
    % Setup and run the OpenSim ID Tool for this trial
    idSetupFileName = [setDir '/' trialFileName 'setupInverseDynamics.xml'];
    idTool = InverseDynamicsTool(idSetupFileName);
    idTool.run();
    
    % Setup and run the OpenSim Anlayze Tool for this trial
    analyzeSetupFileName = [setDir '/' trialFileName 'setupBodyKinematics.xml'];
    analyzeTool = AnalyzeTool(analyzeSetupFileName);
    analyzeTool.run();
    
end
