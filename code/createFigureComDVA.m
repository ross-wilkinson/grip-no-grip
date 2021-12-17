%% Create figure - CoM displacement, velocity, acceleration

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

%% Calculate group mean CoM displacement
aComPosY = NaN(361,nConditions);
[a1, a2, a3, a4] = deal(NaN(361,nSubjects));

measure = 'comPosY';

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

aComPosY(:,1) = mean(a1,2);
aComPosY(:,2) = mean(a2,2);
aComPosY(:,3) = mean(a3,2);
aComPosY(:,4) = mean(a4,2);

%% CoM velocity
aComVelY = NaN(361,nConditions);
[a1, a2, a3, a4] = deal(NaN(361,nSubjects));

measure = 'comVelY';

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

aComVelY(:,1) = mean(a1,2);
aComVelY(:,2) = mean(a2,2);
aComVelY(:,3) = mean(a3,2);
aComVelY(:,4) = mean(a4,2);

%% CoM acceleration
aComAccY = NaN(361,nConditions);
[a1, a2, a3, a4] = deal(NaN(361,nSubjects));

measure = 'comAccY';

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

aComAccY(:,1) = mean(a1,2);
aComAccY(:,2) = mean(a2,2);
aComAccY(:,3) = mean(a3,2);
aComAccY(:,4) = mean(a4,2);

%% Create figure
fig = figure('color','w','position',[50 65 500 580]);

%% Plot CoM Pos Y - seated
ax1 = subplot(421);

plot(aComPosY(:,1),'k','linestyle','-','linewidth',1)
hold on
plot(aComPosY(:,2),'k','linestyle','-','linewidth',0.5)

% Edit axes
box off
xlim([0 380])
ylim([-0.02 0.02])
ylabel({'Vertical CoM','Displacement (m)'})
ax1.YLabel.Position = [-68 0 0];
ax1.XAxisLocation = 'origin';
ax1.XTick = 0:90:360;
ax1.XMinorTick = 'on';
ax1.XAxis.MinorTickValues = 45:90:315;
ax1.YTick = -0.02:0.01:0.02;
ax1.YMinorTick = 'on';
ax1.YAxis.MinorTickValues = -0.015:0.01:0.015;
title('A                  \rmSeated')
ax1.TitleHorizontalAlignment = 'left';
leg = legend('Grip','No Grip');
leg.Box = 'off';
leg.Position = [0.405 0.905 0 0];

%% Plot CoM Pos Y - non-seated
ax2 = subplot(422);

plot(aComPosY(:,3),'k','linestyle','-','linewidth',1)
hold on
plot(aComPosY(:,4),'k','linestyle','-','linewidth',0.5)

% Edit axes
box off
xlim([0 380])
ylim([-0.02 0.02])
% ylabel({'Vertical CoM','Displacement (cm)'})
ax2.XAxisLocation = 'origin';
ax2.XTick = 0:90:360;
ax2.XMinorTick = 'on';
ax2.XAxis.MinorTickValues = 45:90:315;
ax2.YTick = -0.02:0.01:0.02;
ax2.YMinorTick = 'on';
ax2.YAxis.MinorTickValues = -0.015:0.01:0.015;
title('B              \rmNon-Seated')
ax2.TitleHorizontalAlignment = 'left';
leg = legend('Grip','No Grip');
leg.Box = 'off';
leg.Position = [0.855 0.905 0 0];

%% Plot CoM Vel Y - seated
ax3 = subplot(423);

plot(aComVelY(:,1),'k','linestyle','-','linewidth',1)
hold on
plot(aComVelY(:,2),'k','linestyle','-','linewidth',0.5)

% Edit axes
box off
xlim([0 380])
ylim([-0.4 0.4])
ylabel({'Vertical CoM','Velocity (m s^{–1})'})
ax3.YLabel.Position = [-64 0 0];
ax3.XAxisLocation = 'origin';
ax3.XTick = 0:90:360;
ax3.XMinorTick = 'on';
ax3.XAxis.MinorTickValues = 45:90:315;
ax3.YTick = -0.4:0.2:0.4;
ax3.YMinorTick = 'on';
ax3.YAxis.MinorTickValues = -0.3:0.2:0.3;
title('C')
ax3.TitleHorizontalAlignment = 'left';

%% Plot CoM Vel Y - non-seated
ax4 = subplot(424);

plot(aComVelY(:,3),'k','linestyle','-','linewidth',1)
hold on
plot(aComVelY(:,4),'k','linestyle','-','linewidth',0.5)

% Edit axes
box off
xlim([0 380])
ylim([-0.4 0.4])
% ylabel({'Vertical CoM','Velocity (m s^{–1})'})
ax4.XAxisLocation = 'origin';
ax4.XTick = 0:90:360;
ax4.XMinorTick = 'on';
ax4.XAxis.MinorTickValues = 45:90:315;
ax4.YTick = -0.4:0.2:0.4;
ax4.YMinorTick = 'on';
ax4.YAxis.MinorTickValues = -0.3:0.2:0.3;
title('D')
ax4.TitleHorizontalAlignment = 'left';

%% Plot CoM Acc Y - seated
ax5 = subplot(425);

plot(aComAccY(:,1),'k','linestyle','-','linewidth',1)
hold on
plot(aComAccY(:,2),'k','linestyle','-','linewidth',0.5)

% Edit axes
box off
xlim([0 380])
ylim([-10 10])
ylabel({'Vertical CoM','Acceleration (m s^{–2})'})
ax5.YLabel.Position = [-64 0 0];
ax5.XAxisLocation = 'origin';
ax5.XTick = 0:90:360;
ax5.XMinorTick = 'on';
ax5.XAxis.MinorTickValues = 45:90:315;
ax5.YTick = -10:5:10;
ax5.YMinorTick = 'on';
ax5.YAxis.MinorTickValues = -7.5:5:7.5;
title('E')
ax5.TitleHorizontalAlignment = 'left';

%% Plot CoM Acc Y - non-seated
ax6 = subplot(426);

plot(aComAccY(:,3),'k','linestyle','-','linewidth',1)
hold on
plot(aComAccY(:,4),'k','linestyle','-','linewidth',0.5)

% Edit axes
box off
xlim([0 380])
ylim([-10 10])
% ylabel({'Vertical CoM','Acceleration (m s^{–2})'})
% ax6.YLabel.Position = [-64 0 0];
ax6.XAxisLocation = 'origin';
ax6.XTick = 0:90:360;
ax6.XMinorTick = 'on';
ax6.XAxis.MinorTickValues = 45:90:315;
ax6.YTick = -10:5:10;
ax6.YMinorTick = 'on';
ax6.YAxis.MinorTickValues = -7.5:5:7.5;
title('F')
ax6.TitleHorizontalAlignment = 'left';

%% Run SPM analysis on CoM Acc Y - seated and non-seated
t21        = spm1d.stats.ttest_paired(a2', a1');
t43        = spm1d.stats.ttest_paired(a4', a3');

alpha      = 0.05;
nTests     = 1;
p_critical = spm1d.util.p_critical_bonf(alpha, nTests);
t21i       = t21.inference(p_critical, 'two_tailed', true, 'interp',true);
t43i       = t43.inference(p_critical, 'two_tailed', true, 'interp',true);
% t21i       = t21.inference(p_critical, 'two_tailed', false);
% t43i       = t43.inference(p_critical, 'two_tailed', false);
 
%% Plot SPM - seated
ax7 = subplot(427);
t21i.plot();

% Edit axes
box off
xlim([0 380])
ylim([-10 10])
ylabel({'SPM (t)','Vert. CoM Acc.'})
ax7.YLabel.Position = [-64 0 0];
ax7.XAxisLocation = 'origin';
ax7.XTick = 0:90:360;
ax7.XMinorTick = 'on';
ax7.XAxis.MinorTickValues = 45:90:315;
% ax7.YTick = -10:5:10;
% ax7.YMinorTick = 'on';
% ax7.YAxis.MinorTickValues = -7.5:5:7.5;
title('G')
ax7.TitleHorizontalAlignment = 'left';

% Annotate
text(180,-12,'Crank angle (\circ)','HorizontalAlignment','center')

%% Plot SPM - non-seated
ax8 = subplot(428);
t43i.plot();

% Edit axes
box off
xlim([0 380])
ylim([-10 10])
ylabel('')
ax8.YLabel.Position = [-64 0 0];
ax8.XAxisLocation = 'origin';
ax8.XTick = 0:90:360;
ax8.XMinorTick = 'on';
ax8.XAxis.MinorTickValues = 45:90:315;
% ax7.YTick = -10:5:10;
% ax7.YMinorTick = 'on';
% ax7.YAxis.MinorTickValues = -7.5:5:7.5;
title('H')
ax8.TitleHorizontalAlignment = 'left';

% Annotate
text(180,-12,'Crank angle (\circ)','HorizontalAlignment','center')

str1 = 'p = .01';
x1 = t43i.clusters{1}.xy(1);
y1 = t43i.clusters{1}.xy(2)+2.5;
text(x1,y1,str1,'HorizontalAlignment','center')

str2 = 'p < .001';
x2 = t43i.clusters{2}.xy(1);
y2 = t43i.clusters{2}.xy(2)-2.5;
text(x2,y2,str2,'HorizontalAlignment','center')

str3 = 'p = .04';
x3 = t43i.clusters{3}.xy(1);
y3 = t43i.clusters{3}.xy(2)+2.5;
text(x3,y3,str3,'HorizontalAlignment','center')

str4 = 'p < .001';
x4 = t43i.clusters{5}.xy(1);
y4 = t43i.clusters{5}.xy(2)-2.5;
text(x4,y4,str4,'HorizontalAlignment','center')

%% Save figure
cd(docDir)
% export_fig(fig,'fig_ComDVA_grey','-eps','-png','-grey','-r900')
% export_fig(fig,'fig_ComDVA_color','-eps','-png','-cmyk','-r900')
