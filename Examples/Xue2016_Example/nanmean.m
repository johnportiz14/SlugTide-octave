function m=nanmean(x)
m=mean(x(isfinite(x)));
