%% Statistical analysis Grip vs No Grip - Pmax.i

%% Initialize
% clear;clc;close all

%% Set directories
expDir = '/Users/rosswilkinson/Google Drive/projects/grip-no-grip';
datDir = [expDir '/data'];
docDir = [expDir '/docs'];
codDir = [expDir '/matlab'];
resDir = [expDir '/results'];

%% Load data
[~,S] = sliceTablePmax(resDir,0);

%% Create a table storing the responses by column
subjectList = fieldnames(S.time);
nSubjects = numel(subjectList);
nConditions = 4;
a = NaN(nSubjects,nConditions);

measure = 'crankPowerT';

for iSub = 1:nSubjects
    t1 = max(S.(measure).(subjectList{iSub}).c0101);  
    t2 = max(S.(measure).(subjectList{iSub}).c0102);
    t3 = max(S.(measure).(subjectList{iSub}).c0103);
    t4 = max(S.(measure).(subjectList{iSub}).c0201);  
    t5 = max(S.(measure).(subjectList{iSub}).c0202);
    try
        t6 = max(S.(measure).(subjectList{iSub}).c0203);
    catch
        t6 = NaN;
    end
    t7 = max(S.(measure).(subjectList{iSub}).c0301);  
    t8 = max(S.(measure).(subjectList{iSub}).c0302);
    try
        t9 = max(S.(measure).(subjectList{iSub}).c0303);
    catch
        t9 = NaN;
    end
    t10 = max(S.(measure).(subjectList{iSub}).c0401);  
    t11 = max(S.(measure).(subjectList{iSub}).c0402);
    t12 = max(S.(measure).(subjectList{iSub}).c0403);
    
    a(iSub,1) = mean([t1 t2 t3],'omitnan');
    a(iSub,2) = mean([t4 t5 t6],'omitnan');
    a(iSub,3) = mean([t7 t8 t9],'omitnan');
    a(iSub,4) = mean([t10 t11 t12],'omitnan');
end

t = array2table(a,'VariableNames',{'SEG','SENG','STG','STNG'});

%% Create a table reflecting the within-subject factors

factorNames = {'posture','grip'};
factor1 = categorical([ 1; 1; 2; 2]); %posture
factor2 = categorical([ 1; 2; 1; 2]); %grip

within = table(factor1,factor2,'VariableNames',factorNames);

%% Create an interaction factor capturing all combinations of posture and grip
within.posture_grip = within.posture .* within.grip;

%% Fit the repeated measures ANOVA
rm = fitrm(t,'SEG-STNG~1','WithinDesign',within);

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
