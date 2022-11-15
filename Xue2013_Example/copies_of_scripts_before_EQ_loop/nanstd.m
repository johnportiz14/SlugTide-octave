function s=nanstd(x);
s=std(x(isfinite(x)));
