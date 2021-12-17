% Statistical analysis Grip vs No Grip - Pmax cycles

%% Load data structure 
docDir = '/Users/rosswilkinson/Google Drive/projects/grip-no-grip/docs';
resDir = '/Users/rosswilkinson/Google Drive/projects/grip-no-grip/results';
[~,S] = sliceTablePmax(resDir,0);

%% Create a table storing the responses by column
a = NaN(10,4);
subjectList = fieldnames(S.time);

measure1 = 'anklePowerR';
measure2 = 'crankPwrR';
figureName = [measure1 'perc'];
yLabel = {'Ankle power','contribution to crank power (%)'};
yLim = [0 100];

for iSub = 1:10
    t1 = mean(S.(measure1).(subjectList{iSub}).c0101) / mean(S.(measure2).(subjectList{iSub}).c0101);  
    t2 = mean(S.(measure1).(subjectList{iSub}).c0102) / mean(S.(measure2).(subjectList{iSub}).c0102) ;
    t3 = mean(S.(measure1).(subjectList{iSub}).c0103) / mean(S.(measure2).(subjectList{iSub}).c0103);
    t4 = mean(S.(measure1).(subjectList{iSub}).c0201) / mean(S.(measure2).(subjectList{iSub}).c0201);  
    t5 = mean(S.(measure1).(subjectList{iSub}).c0202) / mean(S.(measure2).(subjectList{iSub}).c0202);
    try
        t6 = mean(S.(measure1).(subjectList{iSub}).c0203) / mean(S.(measure2).(subjectList{iSub}).c0203);
    catch
        t6 = NaN;
    end
    t7 = mean(S.(measure1).(subjectList{iSub}).c0301) / mean(S.(measure2).(subjectList{iSub}).c0301);  
    t8 = mean(S.(measure1).(subjectList{iSub}).c0302) / mean(S.(measure2).(subjectList{iSub}).c0302);
    try
        t9 = mean(S.(measure1).(subjectList{iSub}).c0303) / mean(S.(measure2).(subjectList{iSub}).c0303);
    catch
        t9 = NaN;
    end
    t10 = mean(S.(measure1).(subjectList{iSub}).c0401) / mean(S.(measure2).(subjectList{iSub}).c0401);  
    t11 = mean(S.(measure1).(subjectList{iSub}).c0402) / mean(S.(measure2).(subjectList{iSub}).c0402);
    t12 = mean(S.(measure1).(subjectList{iSub}).c0403) / mean(S.(measure2).(subjectList{iSub}).c0403);
    
    a(iSub,1) = mean([t1 t2 t3],'omitnan');
    a(iSub,2) = mean([t4 t5 t6],'omitnan');
    a(iSub,3) = mean([t7 t8 t9],'omitnan');
    a(iSub,4) = mean([t10 t11 t12],'omitnan');
end

a = a * 100; % convert to %

t = array2table(a,'VariableNames',{'cond1','cond2','cond3','cond4'});

%% Create a table reflecting the within-subject factors

factorNames = {'posture','grip'};
factor1 = categorical([ 1; 1; 2; 2]); %posture
factor2 = categorical([ 1; 2; 1; 2]); %grip

within = table(factor1,factor2,'VariableNames',factorNames);

%% Create an interaction factor capturing all combinations of posture and grip
within.posture_grip = within.posture .* within.grip;

%% Fit the repeated measures ANOVA
rm = fitrm(t,'cond1-cond4~1','WithinDesign',within);

%% Run repeated measures analysis of variance
r = ranova(rm,'WithinModel','posture*grip')

%% Run multiple comparisons
c = multcompare(rm,'posture_grip', 'ComparisonType', 'tukey-kramer')

%% Calculate mean and standard deviation in each condition
statstbl = grpstats(rm,'posture_grip')

%% Calculate generalized eta squared for each effect
r.n2g = r.SumSq / sum(r.SumSq(3:2:end));

%% Test for sphericity
tbl = mauchly(rm)

%% Calculate effect size for multiple comparisons
G1 = findgroups(c.posture_grip_1);
G2 = findgroups(c.posture_grip_2);

for i = 1:size(c,1)
    A = t(:,G1(i));
    B = t(:,G2(i));
    c.SdDiff(i) = std(A{:,1} - B{:,1},'omitnan');
end

n = statstbl.GroupCount(1);
df = n-1;
c.d_av = c.Difference ./ c.SdDiff;
c.t = c.Difference ./ c.StdErr;
c.d_z = c.t/sqrt(n);
c.g_av = c.d_av * (1-(3/(4*df-1)));
c.CL = normcdf(c.d_z)

%% Figures
figure('color','w')
%plot seated
subplot(121)
plot(t{:,1:2}','*-')
hold on
plot(mean(t{:,1:2},'omitnan'),'k*-','LineWidth',2)
%edit appearance
xlim([0.5 2.5])
ylim(yLim)
% ylim auto
box off
ylabel(yLabel)
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Grip','No Grip'};
ax.XAxisLocation = 'origin';
title('Seated')

%plot non-seated
subplot(122)
plot(t{:,3:4}','*-')
hold on
plot(mean(t{:,3:4},'omitnan'),'k*-','LineWidth',2)
%edit appearance
xlim([0.5 2.5])
ylim(yLim)
% ylim auto
box off
ylabel('')
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Grip','No Grip'};
ax.XAxisLocation = 'origin';
title('Non-seated')

%% Export figure
cd(docDir)
export_fig(gcf,['fig_' figureName '_color'],'-eps','-cmyk','-r1200')
