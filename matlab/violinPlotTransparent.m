function violinPlotTransparent(meas, meas_pdf, x_pdf, CI95, mu, marker)

nSubjects = numel(meas);
x = 1;
y = meas;
xa = x - meas_pdf;
ya = x_pdf;

patch([xa fliplr(2-xa)], [ya fliplr(ya)],[1 1 1]*1,'facealpha',0.2,'linestyle','none')
hold on
set(gca,'ColorOrderIndex',1)
for i = 1:nSubjects
    scatter(x,y(i),30,'filled',marker)
end

kLow = find(ya<CI95(1),1,'last');
kHigh = find(ya>CI95(2),1,'first');
kMean = find(ya>mu,1,'first');

line([xa(kMean) 2-xa(kMean)],[mu mu],'color','w','linestyle','-')
line([xa(kLow) 2-xa(kLow)],[CI95(1) CI95(1)],'color','w','linestyle',':')
line([xa(kHigh) 2-xa(kHigh)],[CI95(2) CI95(2)],'color','w','linestyle',':')

box off

end