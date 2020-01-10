function [output]=fano_bootstrap(rep1,rep2,rep3)
%bootstrap to determine if the fano factor is likely to be close to one

all_cells=[rep1.Cells(:,[3 5]);rep2.Cells(:,[3 5]);rep3.Cells(:,[3 5])];

for i=1:16
    data=all_cells(all_cells(:,2)==i,1);
    [ci,bootstat]=bootci(10000,@fano_strp,data);
    fanoLower(i)=ci(1);
    fanoUpper(i)=ci(2);
    fanoMean(i)=mean(bootstat);
end
output(:,1)=fanoMean;
output(:,2)=fanoLower;
output(:,3)=fanoUpper;


function [fano]=fano_strp(x)
fano=var(x)./mean(x);
end

end
