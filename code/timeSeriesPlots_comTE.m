% Statistical analysis Grip vs No Grip - Pmax cycles

%% Load data structure 
docDir = '/Users/rosswilkinson/Google Drive/projects/grip-no-grip/docs';
resDir = '/Users/rosswilkinson/Google Drive/projects/grip-no-grip/results';
[~,S] = sliceTablePmax(resDir,0);

%% Create a table storing the responses by column
a = NaN(10,4);
subjectList = fieldnames(S.time);

measure1 = 'comTE';
measure2 = '';
yLabel = '\Delta whole-body mechanical energy (J)';
yLim = [-20 20];

for iSub = 1:11
    t1 = S.(measure1).(subjectList{iSub}).c0101 - mean(S.(measure1).(subjectList{iSub}).c0101);  
    t2 = S.(measure1).(subjectList{iSub}).c0102 - mean(S.(measure1).(subjectList{iSub}).c0102);
    t3 = S.(measure1).(subjectList{iSub}).c0103 - mean(S.(measure1).(subjectList{iSub}).c0103);
    t4 = S.(measure1).(subjectList{iSub}).c0201 - mean(S.(measure1).(subjectList{iSub}).c0201);  
    t5 = S.(measure1).(subjectList{iSub}).c0202 - mean(S.(measure1).(subjectList{iSub}).c0202);
    try
        t6 = S.(measure1).(subjectList{iSub}).c0203 - mean(S.(measure1).(subjectList{iSub}).c0203);
    catch
        t6 = NaN(361,1);
    end
    t7 = S.(measure1).(subjectList{iSub}).c0301 - mean(S.(measure1).(subjectList{iSub}).c0301);  
    t8 = S.(measure1).(subjectList{iSub}).c0302 - mean(S.(measure1).(subjectList{iSub}).c0302);
    try
        t9 = S.(measure1).(subjectList{iSub}).c0303 - mean(S.(measure1).(subjectList{iSub}).c0303);
    catch
        t9 = NaN(361,1);
    end
    t10 = S.(measure1).(subjectList{iSub}).c0401 - mean(S.(measure1).(subjectList{iSub}).c0401);  
    t11 = S.(measure1).(subjectList{iSub}).c0402 - mean(S.(measure1).(subjectList{iSub}).c0402);
    t12 = S.(measure1).(subjectList{iSub}).c0403 - mean(S.(measure1).(subjectList{iSub}).c0403);
    
    a1(:,iSub) = -mean([t1 t2 t3],2,'omitnan');
    a2(:,iSub) = -mean([t4 t5 t6],2,'omitnan');
    a3(:,iSub) = -mean([t7 t8 t9],2,'omitnan');
    a4(:,iSub) = -mean([t10 t11 t12],2,'omitnan');
end

%% Figures
figure('color','w')
%plot seated, grip
subplot(411)
x = 0:360;
plot(x,a1)
hold on
line([180 180],yLim,'LineStyle',':')

%edit appearance
xlim([0 380])
ylim(yLim)
% ylim auto
box off
ylabel(yLabel)
ax = gca;
ax.XTick = 0:45:380;
xlabel('Crank angle (\circ)')
ax.XAxisLocation = 'origin';
title('Seated, Grip')

%plot seated, no grip
subplot(412)
x = 0:360;
plot(x,a2)
hold on
line([180 180],yLim,'LineStyle',':')

%edit appearance
xlim([0 380])
ylim(yLim)
% ylim auto
box off
ylabel(yLabel)
ax = gca;
ax.XTick = 0:45:380;
xlabel('Crank angle (\circ)')
ax.XAxisLocation = 'origin';
title('Seated, Grip')

%plot non-seated, grip
subplot(413)
x = 0:360;
plot(x,a3)
hold on
line([180 180],yLim,'LineStyle',':')

%edit appearance
xlim([0 380])
ylim(yLim)
% ylim auto
box off
ylabel(yLabel)
ax = gca;
ax.XTick = 0:45:380;
xlabel('Crank angle (\circ)')
ax.XAxisLocation = 'origin';
title('Non-Seated, Grip')

%plot non-seated, grip
subplot(414)
x = 0:360;
plot(x,a4)
hold on
line([180 180],yLim,'LineStyle',':')

%edit appearance
xlim([0 380])
ylim(yLim)
% ylim auto
box off
ylabel(yLabel)
ax = gca;
ax.XTick = 0:45:380;
xlabel('Crank angle (\circ)')
ax.XAxisLocation = 'origin';
title('Non-Seated, No Grip')

%% Export figure
cd(docDir)
% export_fig(gcf,['fig_' measure1 '_color'],'-eps','-cmyk','-r1200')
