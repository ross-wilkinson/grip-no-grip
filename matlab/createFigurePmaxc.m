%% Create figure - Pmax.i

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

%% Run stats
run statisticalAnalysis_Pmax_i.m

%% Create a table storing the responses by column
subjectList = fieldnames(S.time);
nSubjects = numel(subjectList);
nConditions = 4;
a = NaN(nSubjects,nConditions);

measure = 'crankPowerTkg';
yLabel = {'Cycle Max.','Crank Power (W kg^{–1})'};

for iSub = 1:nSubjects
    t1 = mean(S.(measure).(subjectList{iSub}).c0101);  
    t2 = mean(S.(measure).(subjectList{iSub}).c0102);
    t3 = mean(S.(measure).(subjectList{iSub}).c0103);
    t4 = mean(S.(measure).(subjectList{iSub}).c0201);  
    t5 = mean(S.(measure).(subjectList{iSub}).c0202);
    try
        t6 = mean(S.(measure).(subjectList{iSub}).c0203);
    catch
        t6 = NaN;
    end
    t7 = mean(S.(measure).(subjectList{iSub}).c0301);  
    t8 = mean(S.(measure).(subjectList{iSub}).c0302);
    try
        t9 = mean(S.(measure).(subjectList{iSub}).c0303);
    catch
        t9 = NaN;
    end
    t10 = mean(S.(measure).(subjectList{iSub}).c0401);  
    t11 = mean(S.(measure).(subjectList{iSub}).c0402);
    t12 = mean(S.(measure).(subjectList{iSub}).c0403);
    
    a(iSub,1) = mean([t1 t2 t3],'omitnan');
    a(iSub,2) = mean([t4 t5 t6],'omitnan');
    a(iSub,3) = mean([t7 t8 t9],'omitnan');
    a(iSub,4) = mean([t10 t11 t12],'omitnan');
end

t = array2table(a,'VariableNames',{'SEGR','SENG','NSGR','NSNG'});

%% Calculate percentage change between Factor 2 within each level of Factor 1 = No Grip vs. Grip
a01 = a;
a01(:,1) = a(:,1) - a(:,1);
a01(:,2) = a(:,2) - a(:,1);
a01(:,3) = a(:,3) - a(:,3);
a01(:,4) = a(:,4) - a(:,3);

aPerc1 = a01;
aPerc1(:,2) = a01(:,2) ./ a(:,1) * 100;
aPerc1(:,4) = a01(:,4) ./ a(:,3) * 100;

%% Calculate percentage change between Factor 1 within each level of Factor 2 = Non-Seated vs. Seated
a02 = a;
a02(:,1) = a(:,1) - a(:,1);
a02(:,2) = a(:,2) - a(:,2);
a02(:,3) = a(:,3) - a(:,1);
a02(:,4) = a(:,4) - a(:,2);

aPerc2 = a02;
aPerc2(:,3) = a02(:,3) ./ a(:,1) * 100;
aPerc2(:,4) = a02(:,4) ./ a(:,2) * 100;

%% Create figure
fig = figure('color','w','position',[50 65 400 500]);

%% Subplot 1 - absolute difference SEGR vs. SENG
ax1 = subplot(221);

plot(a(:,1:2)')
hold on
set(ax1,'ColorOrderIndex',1)
for i = 1:nSubjects
    scatter([1 2],a(i,1:2),20,'filled','o')
end
plot(mean(a(:,1:2),'omitnan'),'k-','LineWidth',2)
scatter([1 2],mean(a(:,1:2),'omitnan'),60,'k','filled','o')

xlim([0.5 2.5])
ylim([8 24])
box off
ylabel(yLabel)
ax1.XTick = 1:2;
ax1.XTickLabel = {'Grip','No Grip'};
ax1.YTick = 8:4:24;
ax1.YMinorTick = 'on';
ax1.YAxis.MinorTickValues = 10:4:22;
title('A. \rmSeated')
ax1.TitleHorizontalAlignment = 'left';

%% Subplot 2 - absolute difference NSGR vs. NSNG
ax2 = subplot(222);

plot(a(:,3:4)')
hold on
set(ax2,'ColorOrderIndex',1)
for i = 1:nSubjects
    scatter([1 2],a(i,3:4),20,'filled','s')
end
plot(mean(a(:,3:4),'omitnan'),'k-','LineWidth',2)
scatter([1 2],mean(a(:,3:4),'omitnan'),60,'k','filled','s')

xlim([0.5 2.5])
ylim([8 24])
box off
% ylabel(yLabel)
ax2.XTick = 1:2;
ax2.XTickLabel = {'Grip','No Grip'};
ax2.YTick = 8:4:24;
ax2.YMinorTick = 'on';
ax2.YAxis.MinorTickValues = 10:4:22;
title('B. \rmNon-Seated')
ax2.TitleHorizontalAlignment = 'left';

%% CalculatePd within Seated
meas = -aPerc1(:,2); % data
[meas_pdf, x_pdf, CI95, mu] = calculatePd(meas);

%% Subplot 3 - % change and pd Seated
ax3 = subplot(223);

% Violin plot
violinPlot(meas, meas_pdf, x_pdf, CI95, mu, 'o')

% Edit axes
xlim([0.9 1.1])
ylim([-10 30])
ax3.XTick = [];
ax3.XAxisLocation = 'origin';
ylabel({'% Change','Grip \itvs. \rmNo Grip'})
title('C. \rmSeated')
ax3.TitleHorizontalAlignment = 'left';

% Annotate
text(1.03,32,{'ES = 1.4','p = .002'})

%% CalculatePd within Non-Seated
meas = -aPerc1(:,4); % data
[meas_pdf, x_pdf, CI95, mu] = calculatePd(meas);

%% Subplot 4 - % change and pd Non-Seated
ax4 = subplot(224);

% Violin plot
violinPlot(meas, meas_pdf, x_pdf, CI95, mu, 's')

% Edit axes
xlim([0.9 1.1])
ylim([-10 30])
ax4.XTick = [];
ax4.XAxisLocation = 'origin';
% ylabel({'% Change','No Grip \itvs. \rmGrip'})
title('D. \rmNon-Seated')
ax4.TitleHorizontalAlignment = 'left';

% Annotate
text(1.03,32,{'ES = 2.7','p < .001'})

%% Save figure
cd(docDir)
export_fig(fig,'fig_Pmaxc_grey','-eps','-png','-grey','-r900')
export_fig(fig,'fig_Pmaxc_color','-eps','-png','-cmyk','-r900')
