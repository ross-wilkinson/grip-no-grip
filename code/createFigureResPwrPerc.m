%% Create figure - % contribution of lower- and upper-body power to crank power

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
aResPwr = NaN(nSubjects,nConditions);
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

aResPwr(:,1) = mean(a1);
aResPwr(:,2) = mean(a2);
aResPwr(:,3) = mean(a3);
aResPwr(:,4) = mean(a4);

%% Leg power
aLegPwr = NaN(nSubjects,nConditions);
[a1, a2, a3, a4] = deal(NaN(361,nSubjects));

measure = 'legPowerT';

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

aLegPwr(:,1) = mean(a1);
aLegPwr(:,2) = mean(a2);
aLegPwr(:,3) = mean(a3);
aLegPwr(:,4) = mean(a4);

%% Crank power
aCrankPwr = NaN(nSubjects,nConditions);
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

aCrankPwr(:,1) = mean(a1);
aCrankPwr(:,2) = mean(a2);
aCrankPwr(:,3) = mean(a3);
aCrankPwr(:,4) = mean(a4);

%% Calculate % contributions
aLegPwrPerc = aLegPwr ./ aCrankPwr * 100;
aResPwrPerc = 100 - aLegPwrPerc;

%% Create figure
fig = figure('color','w','position',[50 65 400 500]);

%% Absolute difference SEGR vs. SENG - Upper-body
ax1 = subplot(221);

plot(aResPwr(:,1:2)')
hold on
set(ax1,'ColorOrderIndex',1)
for i = 1:nSubjects
    scatter([1 2],aResPwr(i,1:2),20,'filled','o')
end
plot(mean(aResPwr(:,1:2),'omitnan'),'k-','LineWidth',2)
scatter([1 2],mean(aResPwr(:,1:2),'omitnan'),60,'k','filled','o')

xlim([0.5 2.5])
ylim([-100 300])
box off
ylabel({'Cycle Max.','Upper-Body Power (W)'})
ax1.XAxisLocation = 'origin';
ax1.XTick = 1:2;
sp = '\newline\newline';
ax1.XTickLabel = {[sp 'Grip'],[sp 'No Grip']};
ax1.YTick = -100:100:300;
ax1.YMinorTick = 'on';
ax1.YAxis.MinorTickValues = -50:100:250;
title('A. \rmSeated')
ax1.TitleHorizontalAlignment = 'left';

%% Inset: % crank power u-b seated
ax5 = axes('Position',[0.37 0.82 0.1 0.1]);

y = aResPwrPerc(:,1:2);

plot(mean(y,'omitnan'),'k-','LineWidth',2)
hold on
scatter([1 2],mean(y,'omitnan'),60,'k','filled','o')

xlim([0.5 2.5])
ylim([3 7])
box off
ylabel('% Crank Power')
ax5.XTick = 1:2;
sp = '';
ax5.XTickLabel = {[sp 'G'],[sp 'NG']};
ax5.YTick = 3:7;
ax5.YMinorTick = 'on';
ax5.YAxis.MinorTickValues = 3.5:6.5;
ax5.XAxis.FontSize = 8;
ax5.YAxis.FontSize = 8;

%% Absolute difference NSGR vs. NSNG - Upper-body
ax2 = subplot(222);

plot(aResPwr(:,3:4)')
hold on
set(ax2,'ColorOrderIndex',1)
for i = 1:nSubjects
    scatter([1 2],aResPwr(i,3:4),20,'filled','s')
end
plot(mean(aResPwr(:,3:4),'omitnan'),'k-','LineWidth',2)
scatter([1 2],mean(aResPwr(:,3:4),'omitnan'),60,'k','filled','s')

xlim([0.5 2.5])
ylim([-100 300])
box off
% ylabel({'Cycle Max.','Upper-Body Power (W)'})
ax2.XAxisLocation = 'origin';
ax2.XTick = 1:2;
sp = '';
ax2.XTickLabel = {[sp 'Grip'],[sp 'No Grip']};
ax2.YTick = -100:100:300;
ax2.YMinorTick = 'on';
ax2.YAxis.MinorTickValues = -50:100:250;
title('B. \rmNon-Seated')
ax2.TitleHorizontalAlignment = 'left';

%% Inset: % crank power u-b non-seated
ax6 = axes('Position',[0.82 0.82 0.1 0.1]);

y = aResPwrPerc(:,3:4);

plot(mean(y,'omitnan'),'k-','LineWidth',2)
hold on
scatter([1 2],mean(y,'omitnan'),60,'k','filled','o')

xlim([0.5 2.5])
ylim([5 9])
box off
ylabel('% Crank Power')
ax6.XTick = 1:2;
sp = '';
ax6.XTickLabel = {[sp 'G'],[sp 'NG']};
ax6.YTick = 5:9;
ax6.YMinorTick = 'on';
ax6.YAxis.MinorTickValues = 5.5:8.5;
ax6.XAxis.FontSize = 8;
ax6.YAxis.FontSize = 8;

%% Absolute difference SEGR vs. SENG - Leg
ax3 = subplot(223);

y = aLegPwr(:,1:2);

plot(y')
hold on
set(ax3,'ColorOrderIndex',1)
for i = 1:nSubjects
    scatter([1 2],y(i,:),20,'filled','o')
end
plot(mean(y,'omitnan'),'k-','LineWidth',2)
scatter([1 2],mean(y,'omitnan'),60,'k','filled','o')

xlim([0.5 2.5])
ylim([400 2000])
box off
ylabel('Cycle Max. Leg Power (W)')
ax3.XAxisLocation = 'origin';
ax3.XTick = 1:2;
sp = '';
ax3.XTickLabel = {[sp 'Grip'],[sp 'No Grip']};
ax3.YTick = 400:400:2000;
ax3.YMinorTick = 'on';
ax3.YAxis.MinorTickValues = 600:400:1800;
title('C. \rmSeated')
ax3.TitleHorizontalAlignment = 'left';

%% Inset: % crank power leg seated
ax7 = axes('Position',[0.37 0.34 0.1 0.1]);

y = aLegPwrPerc(:,1:2);

plot(mean(y,'omitnan'),'k-','LineWidth',2)
hold on
scatter([1 2],mean(y,'omitnan'),60,'k','filled','o')

xlim([0.5 2.5])
ylim([93 97])
box off
ylabel('% Crank Power')
ax7.XTick = 1:2;
sp = '';
ax7.XTickLabel = {[sp 'G'],[sp 'NG']};
ax7.YTick = 93:97;
ax7.YMinorTick = 'on';
ax7.YAxis.MinorTickValues = 93.5:96.5;
ax7.XAxis.FontSize = 8;
ax7.YAxis.FontSize = 8;

%% Absolute difference NSGR vs. NSNG - Leg
ax4 = subplot(224);

y = aLegPwr(:,3:4);

plot(y')
hold on
set(ax4,'ColorOrderIndex',1)
for i = 1:nSubjects
    scatter([1 2],y(i,:),20,'filled','s')
end
plot(mean(y,'omitnan'),'k-','LineWidth',2)
scatter([1 2],mean(y,'omitnan'),60,'k','filled','s')

xlim([0.5 2.5])
ylim([400 2000])
box off
% ylabel({'Cycle Max.','Leg Power (W)'})
ax4.XAxisLocation = 'origin';
ax4.XTick = 1:2;
sp = '';
ax4.XTickLabel = {[sp 'Grip'],[sp 'No Grip']};
ax4.YTick = 400:400:2000;
ax4.YMinorTick = 'on';
ax4.YAxis.MinorTickValues = 600:400:1800;
title('D. \rmNon-Seated')
ax4.TitleHorizontalAlignment = 'left';

%% Inset: % crank power u-b non-seated
ax8 = axes('Position',[0.82 0.34 0.1 0.1]);

y = aLegPwrPerc(:,3:4);

plot(mean(y,'omitnan'),'k-','LineWidth',2)
hold on
scatter([1 2],mean(y,'omitnan'),60,'k','filled','o')

xlim([0.5 2.5])
ylim([91 95])
box off
ylabel('% Crank Power')
ax8.XTick = 1:2;
sp = '';
ax8.XTickLabel = {[sp 'G'],[sp 'NG']};
ax8.YTick = 91:95;
ax8.YMinorTick = 'on';
ax8.YAxis.MinorTickValues = 91.5:94.5;
ax8.XAxis.FontSize = 8;
ax8.YAxis.FontSize = 8;

%% Save figure
cd(docDir)
% export_fig(fig,'fig_UbPwrPerc_grey','-eps','-png','-grey','-r900')
% export_fig(fig,'fig_UbPwrPerc_color','-eps','-png','-cmyk','-r900')
