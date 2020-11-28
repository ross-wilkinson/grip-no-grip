function mat2trc(filename)
%MAT2TRC Converts motion capture data from MatLab structure to trc file
%   Input:
%           - .mat file containing motion capture data
%
%   Output:
%           - .trc file containing motion capture data
%

%% Load .mat file
data = load([filename '.mat']);

%% Set labels and trajectories to vars
traj_labels = data.(filename).Trajectories.Labeled.Labels;
traj_data = data.(filename).Trajectories.Labeled.Data;

%% Reshape marker data
traj_data_trc = permute(traj_data,[2 3 1]);

%% Transform data from UQ Qualisys XYZ to OpenSim XYZ
% UQ X = OS X
% UQ Y = OS -Z
% UQ Z = OS Y

temp = traj_data_trc;

traj_data_trc(2,:,:) = temp(3,:,:); % UQ Z to OS Y
traj_data_trc(3,:,:) = -temp(2,:,:); % UQ -Y to OS Z

%% Check if analog data was collected
analogCheck = isfield(data.(filename),'Analog');

%% Create new .trc file
newfilename = [filename '.trc']; %add .trc file ext to filename
pathName = cd;
fid = fopen([pathName '/' newfilename],'w'); %open the file

if analogCheck == 1
    %% Create header information for .trc file
    fprintf(fid,'PathFileType\t4\t(X/Y/Z)\t %s\n',newfilename);
    fprintf(fid,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
    fprintf(fid,'%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n',...
    data.(filename).Analog.Frequency,...
    data.(filename).FrameRate,...
    data.(filename).Frames,...
    data.(filename).Trajectories.Labeled.Count,...
    'mm',...
    data.(filename).Analog.Frequency,...
    data.(filename).StartFrame,...
    data.(filename).Frames);
else
    %% Create header information for .trc file
    fprintf(fid,'PathFileType\t4\t(X/Y/Z)\t %s\n',newfilename);
    fprintf(fid,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
    fprintf(fid,'%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n',...
    data.(filename).FrameRate,...
    data.(filename).FrameRate,...
    data.(filename).Frames,...
    data.(filename).Trajectories.Labeled.Count,...
    'mm',...
    data.(filename).FrameRate,...
    data.(filename).StartFrame,...
    data.(filename).Frames);
end

%% Create data frame for .trc file
header4 = 'Frame#\tTime\t'; %create header4 using marker labels
header5 = '\t\t'; %create header5 using XYZ and marker numbers
format_text = '%i\t%2.4f\t';
start_frame = data.(filename).StartFrame;
end_frame = data.(filename).Frames + start_frame - 1;
nframes = start_frame:end_frame;
frate = data.(filename).FrameRate;
time = nframes/frate;
data_out = [nframes; time];

for mm = 1:length(traj_labels)
    header4 = [header4 traj_labels{mm} '\t\t\t'];
    header5 = [header5 'X' num2str(mm) '\t' 'Y' num2str(mm) '\t' 'Z' num2str(mm) '\t'];
    format_text = [format_text '%f\t%f\t%f\t'];
    data_out = [data_out;traj_data_trc(1:3,:,mm)];
end

header4 = [header4 '\n'];
header5 = [header5 '\n'];
format_text = [format_text '\n'];

%% Write header and data to .trc file
fprintf(fid,header4);
fprintf(fid,header5);
fprintf(fid,format_text,data_out);
fclose(fid);

end

