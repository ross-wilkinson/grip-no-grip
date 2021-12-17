function [meas_pdf, x_pdf, CI95, mu] = calculatePd(meas)
n = sum(~isnan(meas)); % sample size
df = n-1; % Degrees of freedom

[~,p,~,~] = jbtest(meas,[],0.001);

if p < 0.01 % non-paramteric: kernel dist., median, iqr
    pd = fitdist(meas,'kernel');
    mu = median(meas,'omitnan');
    Cd = [-1.57 1.57]*iqr(meas) / sqrt(n);
    CI95 = mu + Cd;
    x_pdf = linspace(mu-iqr(meas)*4,mu+iqr(meas)*4,1001);
    meas_pdf = pdf(pd,x_pdf);
else % parametric: normal dist., mean, sd
    pd = fitdist(meas,'kernel');
    mu = mean(meas,'omitnan');
    sigma1 = std(meas,'omitnan');
    ux = sigma1 / sqrt(n); % SEM
    Rux = ux / mu * 100; % Relative SEM
    tc = tinv([0.025 0.975],df); % t-score
    Cd = tc * ux; % confidence interval
    CI95 = mu + Cd; % 95% confidence interval
    x_pdf = linspace(mu-sigma1*4,mu+sigma1*4,1001);
    meas_pdf = pdf(pd,x_pdf);
end

end