%% Create figure - Upper-body power

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

%% Set n for subjects and conditions
subjectList = fieldnames(S.time);
nSubjects = numel(subjectList);
nConditions = 4;

%% Upper-body power
aResPwr = NaN(361,nConditions);
[a1, a2, a3, a4] = deal(NaN(361,nSubjects));

measure = 'residualPower';

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

aResPwr(:,1) = mean(a1,2);
aResPwr(:,2) = mean(a2,2);
aResPwr(:,3) = mean(a3,2);
aResPwr(:,4) = mean(a4,2);

%% Create figure
fig = figure('color','w','position',[50 65 250 290]);

%% Subplot 1 - U-B Power Seated
ax1 = subplot(211);

plot(aResPwr(:,1),'k','linestyle','-','linewidth',1)
hold on
plot(aResPwr(:,2),'k','linestyle','-','linewidth',0.5)

% Edit axes
box off
xlim([0 380])
ylim([-500 500])
ylabel({'Upper-Body','Mechanical Power (W)'})
% ax1.YLabel.Position = [-68 0 0];
ax1.XAxisLocation = 'origin';
ax1.XTick = 0:90:360;
ax1.XMinorTick = 'on';
ax1.XAxis.MinorTickValues = 45:90:315;
ax1.YTick = -500:250:500;
ax1.YMinorTick = 'on';
ax1.YAxis.MinorTickValues = -250:500:250;
title('A. \rmSeated')
ax1.TitleHorizontalAlignment = 'left';
leg = legend('Grip','No Grip');
leg.Box = 'off';
leg.Position = [0.66 0.91 0 0];

%% Subplot 1 - U-B Power Non-Seated
ax2 = subplot(212);

plot(aResPwr(:,3),'k','linestyle','-','linewidth',1)
hold on
plot(aResPwr(:,4),'k','linestyle','-','linewidth',0.5)

% Edit axes
box off
xlim([0 380])
ylim([-500 500])
ylabel({'Upper-Body','Mechanical Power (W)'})
% ax1.YLabel.Position = [-68 0 0];
ax2.XAxisLocation = 'origin';
ax2.XTick = 0:90:360;
ax2.XMinorTick = 'on';
ax2.XAxis.MinorTickValues = 45:90:315;
ax2.YTick = -500:250:500;
ax2.YMinorTick = 'on';
ax2.YAxis.MinorTickValues = -250:500:250;
title('B. \rmNon-Seated')
ax2.TitleHorizontalAlignment = 'left';
leg = legend('Grip','No Grip');
leg.Box = 'off';
leg.Position = [0.64 0.92 0 0];

%% Save figure
cd(docDir)
export_fig(fig,'fig_UbPwr_grey','-eps','-png','-grey','-r900')
export_fig(fig,'fig_UbPwr_color','-eps','-png','-cmyk','-r900')
