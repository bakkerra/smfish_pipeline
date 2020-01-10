function [mean_space]=combine_discs(varargin)

cell_list=[0 0];
for i=1:nargin
    current=varargin{1};
    cell_list=[cell_list;current.Cells(:,2),current.Cells(:,3)];
end

for i=1:16
    locs=cell_list(:,1)==i;
    in_bin=cell_list(locs,2);
    se=std(in_bin)./sqrt(length(in_bin));
    mean_space(i,:)=[i,mean(in_bin),2.*se];
end