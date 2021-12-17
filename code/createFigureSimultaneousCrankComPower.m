%% Create figure - Simultaneous crank and CoM power

%% Initialize
clear;clc;close all

%% Set directories
expDir = '/Users/rosswilkinson/Google Drive/projects/grip-no-grip';
datDir = [expDir '/data'];
docDir = [expDir '/docs'];
codDir = [expDir '/matlab'];
resDir = [expDir '/results'];

%% Load data 
[~,S] = sliceTablePmax(resDir,0);

%% Calculate group mean total crank power
subjectList = fieldnames(S.time);
nSubjects = numel(subjectList);
nConditions = 4;
aCrankPowerT = NaN(361,nConditions);
[a1, a2, a3, a4] = deal(NaN(361,nSubjects));

measure = 'crankPowerT';

for iSub = 1:nSubjects
    t1 = S.(measure).(subjectList{iSub}).c0101;  
    t2 = S.(measure).(subjectList{iSub}).c0102;
    t3 = S.(measure).(subjectList{iSub}).c0103;
    t4 = S.(measure).(subjectList{iSub}).c0201;  
    t5 = S.(measure).(subjectList{iSub}).c0202;
    try
        t6 = S.(measure).(subjectList{iSub}).c0203;
    catch
        t6 = NaN(361,1);
    end
    t7 = S.(measure).(subjectList{iSub}).c0301;  
    t8 = S.(measure).(subjectList{iSub}).c0302;
    try
        t9 = S.(measure).(subjectList{iSub}).c0303;
    catch
        t9 = NaN(361,1);
    end
    t10 = S.(measure).(subjectList{iSub}).c0401;  
    t11 = S.(measure).(subjectList{iSub}).c0402;
    t12 = S.(measure).(subjectList{iSub}).c0403;
    
    a1(:,iSub) = mean([t1 t2 t3],2,'omitnan');
    a2(:,iSub) = mean([t4 t5 t6],2,'omitnan');
    a3(:,iSub) = mean([t7 t8 t9],2,'omitnan');
    a4(:,iSub) = mean([t10 t11 t12],2,'omitnan');
end

aCrankPowerT(:,1) = mean(a1,2);
aCrankPowerT(:,2) = mean(a2,2);
aCrankPowerT(:,3) = mean(a3,2);
aCrankPowerT(:,4) = mean(a4,2);

%% Calculate group mean CoM power
aComPower = NaN(361,nConditions);
[a1, a2, a3, a4] = deal(NaN(361,nSubjects));

measure = 'comPower';

for iSub = 1:nSubjects
    t1 = S.(measure).(subjectList{iSub}).c0101 - mean(S.(measure).(subjectList{iSub}).c0101);  
    t2 = S.(measure).(subjectList{iSub}).c0102 - mean(S.(measure).(subjectList{iSub}).c0102);
    t3 = S.(measure).(subjectList{iSub}).c0103 - mean(S.(measure).(subjectList{iSub}).c0103);
    t4 = S.(measure).(subjectList{iSub}).c0201 - mean(S.(measure).(subjectList{iSub}).c0201);  
    t5 = S.(measure).(subjectList{iSub}).c0202 - mean(S.(measure).(subjectList{iSub}).c0202);
    try
        t6 = S.(measure).(subjectList{iSub}).c0203 - mean(S.(measure).(subjectList{iSub}).c0203);
    catch
        t6 = NaN(361,1);
    end
    t7 = S.(measure).(subjectList{iSub}).c0301 - mean(S.(measure).(subjectList{iSub}).c0301);  
    t8 = S.(measure).(subjectList{iSub}).c0302 - mean(S.(measure).(subjectList{iSub}).c0302);
    try
        t9 = S.(measure).(subjectList{iSub}).c0303 - mean(S.(measure).(subjectList{iSub}).c0303);
    catch
        t9 = NaN(361,1);
    end
    t10 = S.(measure).(subjectList{iSub}).c0401 - mean(S.(measure).(subjectList{iSub}).c0401);  
    t11 = S.(measure).(subjectList{iSub}).c0402 - mean(S.(measure).(subjectList{iSub}).c0402);
    t12 = S.(measure).(subjectList{iSub}).c0403 - mean(S.(measure).(subjectList{iSub}).c0403);
    
    a1(:,iSub) = mean([t1 t2 t3],2,'omitnan');
    a2(:,iSub) = mean([t4 t5 t6],2,'omitnan');
    a3(:,iSub) = mean([t7 t8 t9],2,'omitnan');
    a4(:,iSub) = mean([t10 t11 t12],2,'omitnan');
end

aComPower(:,1) = mean(a1,2);
aComPower(:,2) = mean(a2,2);
aComPower(:,3) = mean(a3,2);
aComPower(:,4) = mean(a4,2);

%% Create figure
fig = figure('color','w','position',[50 65 250 600]);

%% Plot crank and com power - seated, grip
ax1 = subplot(411);

plot(aCrankPowerT(:,1))
hold on
plot(aComPower(:,1))

% Edit axes
box off
xlim([0 380])
ylim([-300 1700])
ylabel('Power (W)')
ax1.XAxisLocation = 'origin';
ax1.XTick = 0:90:360;
ax1.YTick = -300:300:1500;
ax1.YMinorTick = 'on';
ax1.YAxis.MinorTickValues = 10:4:22;
title('A. \rmSeated, Grip')
ax1.TitleHorizontalAlignment = 'left';

%% Plot crank and com power - seated, no grip
ax2 = subplot(412);

plot(aCrankPowerT(:,2))
hold on
plot(aComPower(:,2))

% Edit axes
box off
xlim([0 380])
ylim([-300 1700])
ylabel('Power (W)')
ax2.XAxisLocation = 'origin';
ax2.XTick = 0:90:360;
ax2.YTick = -300:300:1500;
ax2.YMinorTick = 'on';
ax2.YAxis.MinorTickValues = 10:4:22;
title('B. \rmSeated, No Grip')
ax.TitleHorizontalAlignment = 'left';

