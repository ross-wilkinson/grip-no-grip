%% Statistical analysis Grip vs No Grip - Pmax.c

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

%% Create a table storing the responses by column
subjectList = fieldnames(S.time);
nSubjects = numel(subjectList);
nConditions = 4;
a = NaN(nSubjects,nConditions);

measure = 'crankPowerT';

df = NaN(nSubjects*12,4);

sub = ones(12,1);
rows = 1:12;

n = 0;
for iSub = 1:nSubjects
    t1 = mean(S.(measure).(subjectList{iSub}).c0101);  
    df(rows(1)+n,:) = [iSub, 1, 1, t1];
    t2 = mean(S.(measure).(subjectList{iSub}).c0102);
    df(rows(2)+n,:) = [iSub, 1, 2, t2];
    t3 = mean(S.(measure).(subjectList{iSub}).c0103);
    df(rows(3)+n,:) = [iSub, 1, 3, t3];
    t4 = mean(S.(measure).(subjectList{iSub}).c0201);  
    df(rows(4)+n,:) = [iSub, 2, 1, t4];
    t5 = mean(S.(measure).(subjectList{iSub}).c0202);
    df(rows(5)+n,:) = [iSub, 2, 2, t5];
    try
        t6 = mean(S.(measure).(subjectList{iSub}).c0203);
        df(rows(6)+n,:) = [iSub, 2, 3, t6];
    catch
        t6 = NaN;
        df(rows(6)+n,:) = [iSub, 2, 3, t6];
    end
    t7 = mean(S.(measure).(subjectList{iSub}).c0301);  
    df(rows(7)+n,:) = [iSub, 3, 1, t7];
    t8 = mean(S.(measure).(subjectList{iSub}).c0302);
    df(rows(8)+n,:) = [iSub, 3, 2, t8];
    try
        t9 = mean(S.(measure).(subjectList{iSub}).c0303);
        df(rows(9)+n,:) = [iSub, 3, 3, t9];
    catch
        t9 = NaN;
        df(rows(9)+n,:) = [iSub, 3, 3, t9];
    end
    t10 = mean(S.(measure).(subjectList{iSub}).c0401);  
    df(rows(10)+n,:) = [iSub, 4, 1, t10];
    t11 = mean(S.(measure).(subjectList{iSub}).c0402);
    df(rows(11)+n,:) = [iSub, 4, 2, t11];
    t12 = mean(S.(measure).(subjectList{iSub}).c0403);
    df(rows(12)+n,:) = [iSub, 4, 3, t12];
    
    n = n + 12;
    
    a(iSub,1) = mean([t1 t2 t3],'omitnan');
    a(iSub,2) = mean([t4 t5 t6],'omitnan');
    a(iSub,3) = mean([t7 t8 t9],'omitnan');
    a(iSub,4) = mean([t10 t11 t12],'omitnan');
    
end

a0 = a;
a0(:,1) = a(:,1) - a(:,1);
a0(:,2) = a(:,2) - a(:,1);
a0(:,3) = a(:,3) - a(:,3);
a0(:,4) = a(:,4) - a(:,3);
aPerc = a0;
aPerc(:,2) = a0(:,2) ./ a(:,1) * 100;
aPerc(:,4) = a0(:,4) ./ a(:,3) * 100;

t = array2table(a,'VariableNames',{'SEGR','SENG','NSGR','NSNG'});

%% Create a table reflecting the within-subject factors

factorNames = {'posture','grip'};
factor1 = categorical([ 1; 1; 2; 2]); %posture
factor2 = categorical([ 1; 2; 1; 2]); %grip

within = table(factor1,factor2,'VariableNames',factorNames);

%% Create an interaction factor capturing all combinations of posture and grip
within.posture_grip = within.posture .* within.grip;

%% Fit the repeated measures ANOVA
rm = fitrm(t,'SEGR-NSNG~1','WithinDesign',within);

%% Run repeated measures analysis of variance
r = ranova(rm,'WithinModel','posture*grip');

%% Run multiple comparisons
c = multcompare(rm,'posture_grip', 'ComparisonType', 'tukey-kramer');

%% Calculate mean and standard deviation in each condition
statstbl = grpstats(rm,'posture_grip');

%% Calculate generalized eta squared for each effect
r.n2g = r.SumSq / sum(r.SumSq(3:2:end));

%% Test for sphericity
tbl = mauchly(rm);

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
