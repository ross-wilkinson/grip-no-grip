function [t2, S] = sliceTablePmax(path,flag)
% Analyze grip vs no grip data

cd(path)

dataFileList = dir('*.mat');
nFiles = size(dataFileList,1);

for iFiles = 1:nFiles
    dataFileName = dataFileList(iFiles).name;
    trialName = strrep(dataFileName,'dataTable.mat','');
    subject = trialName(1:3);
    condition = trialName(4:8);

    data = load(dataFileName);

    % find index of Pmax
    [~,I] = max(data.t_.crankPowerT);

    % find index of each tdc
    [~, locs] = findpeaks(data.t_.crankAngleR,'minPeakDistance',50,'minPeakHeight',5);
    tdcStart = max(locs(locs < I));
    tdcEnd = min(locs(locs > I));

    % slice table at cycle containing index of Pmax
    v = table2array(data.t_(tdcStart:tdcEnd,:));

    % interpolate slice to 361 data points
    x = tdcStart:tdcEnd; % sample points
    xq = tdcStart:(length(v)-1)/360:tdcEnd; % new query points
    method = 'pchip';
    vq = interp1(x,v,xq,method);

    varNames = genvarname(data.t_.Properties.VariableNames);
    t2 = array2table(vq,'VariableNames',varNames);

    for iVars = 1:length(varNames)
        S.(varNames{iVars}).(subject).(condition) = t2.(varNames{iVars});
    end
    
    if flag
        % write excel spreadsheet
        writetable(t2,[path '/' trialName 'dataPmax.xlsx'])
    end
end

end

