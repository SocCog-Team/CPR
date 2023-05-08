function [auroc0]= f_auroc(arr1,arr2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	 

arr1(isnan(arr1)) = [];
arr2(isnan(arr2)) = [];

nb = 500; %%number bins
nu1 = length(arr1);
nu2 = length(arr2);
maxr = max([max(arr1) max(arr2)]);
minr = min([min(arr1) min(arr2)]);
ib = minr-(maxr-minr)/nb:(maxr-minr)/nb:maxr+(maxr-minr)/nb;
n1 = histc(arr1,ib)/nu1;%%probability distr of 1
n2 = histc(arr2,ib)/nu2;%%probability distr of 2
mk = length(n1);
auroc0 =0;
if mk>0
    for jj = 1:mk
        auroc0 = auroc0 + n1(jj)*sum(n2(jj+1:mk))+(n1(jj)*n2(jj))/2; %%the term (n1(jj)*n2(jj))/2 is to account for when the bin for 1 and 2 are the same
    end
else
    auroc0=0.5;
end
