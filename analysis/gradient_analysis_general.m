function [mRNA_nums]=gradient_analysis_general(mRNA_nums,txn_nums,binsize)
%analyzes smfish data in bins on the x coordinate without any reference as
%zero point.  I use only on individual discs at a time

%make bin index based on x coordinate of image
F=ceil(mRNA_nums.Centroid(:,2)./binsize);
mRNA_nums.Bins=F;
txn_nums.Bins=F;

% F=mRNA_nums.Bins;
numbins=max(F);
x=1:numbins;%for graphing per bin data

%count cells per bin
for i=1:numbins
    locs=F==i;
    in_bin=mRNA_nums.Totals(locs);
    cellnumber(i)=length(in_bin);
end
figure;
scatter(x,cellnumber,100,[0 0 0],'filled')
set(gca,'FontSize',24,'LineWidth',3,'FontName','Arial');
xlabel('Bin Number');
ylabel('Number of Cells per Bin');

%make mRNA/cell histograms for each bin
figure;
for i=1:numbins
    locs=F==i;
    in_bin=mRNA_nums.Totals(locs);
    subplot(2,ceil(numbins./2),i)
    histogram(in_bin,'Normalization','Probability','binwidth',2,'facecolor','black')
    set(gca,'FontSize',12,'LineWidth',3,'FontName','Arial');
    xlim([0 100])
    ylim([0 .2])
end
figure;
% distributionPlot(mRNA_nums.Totals,'groups',mRNA_nums.Bins)
% set(gca,'FontSize',24,'LineWidth',3,'FontName','Arial');
% ylim([0 20])
%locate txn sites in each bin and graph their intensity by bin
figure;
sites_matrix=struct;sites_matrix.Bin=[];sites_matrix.Intensity=[];
for i=1:numbins
    locs=F==i;
    in_bin=txn_nums.totalnums(locs);
    site_list=[];
    k=[];
    for j=1:length(in_bin)
        sites=in_bin(j);
        siteses=sites{1,1};
        for jj=1:length(siteses)
            site_list=[site_list;siteses(jj)];
        end
    end
    k(1:length(site_list),1)=i;
    sites_matrix.Bin=[sites_matrix.Bin;k];
    sites_matrix.Intensity=[sites_matrix.Intensity;site_list];
    avg_txn_intensity(i)=mean(site_list);
    subplot(2,ceil(numbins./2),i)
    histogram(site_list,'binwidth',1)
    xlim([0 20])
    set(gca,'FontSize',12,'LineWidth',3,'FontName','Arial');
end
x=1:numbins;
figure;
scatter(x,avg_txn_intensity,100,[0 0 0],'filled')
set(gca,'FontSize',24,'LineWidth',3,'FontName','Arial');
xlabel('Bin Number');
ylabel('Mean Intensity of Transcription Site per Bin');

%frequency of cells w/txn site per bin
for i=1:numbins
    locs=F==i;
    in_bin=txn_nums.Totals(locs);
%     freq_txn(i)=nnz(in_bin)./length(in_bin);
    bootstat=bootstrp(10000,@proportion,in_bin);
    freq_txn(i)=mean(bootstat);
    err1(i)=2*std(bootstat);
end

figure;
scatter(x,freq_txn,100,[0 0 0],'filled');hold on
errorbar(x,freq_txn,err1,'Linestyle','none','CapSize',5,'Color','black');
set(gca,'FontSize',24,'LineWidth',3,'FontName','Arial');
xlabel('Bin Numbers');
ylabel('Frequency of Cell w/ Txn Site');

%frequency of cell w/txn site by mRNA/cell for each bin
% figure;
% for i=1:numbins
%     locs=F==i;
%     in_bin=struct;in_bin_txn=struct;
%     in_bin.Totals=mRNA_nums.Totals(locs);
%     in_bin_txn.Sum=txn_nums.Sum(locs);
%     [x1,y1,err,errx]=proportion_graph(in_bin,in_bin_txn);
%     subplot(2,ceil(numbins./2),i)
%     scatter(x1,y1,100,[0 0 0],'filled');hold on;
%     errorbar(x1,y1,err,'Linestyle','none','CapSize',5,'Color','black');
%     errorbar(x1,y1,errx,'Horizontal','Linestyle','none','CapSize',5,'Color','black');
%     set(gca,'FontSize',12,'LineWidth',3,'FontName','Arial');
%     xlabel('Mean mRNA/Cell in Bin');xlim([0 100]);
%     ylabel('Frequency Txn Cells');ylim([0 .7]);
% end

%graph each txn site vs mRNA/cell in each bin
% figure;
% for i=1:numbins
%     locs=F==i;
%     in_bin=struct;in_bin_txn=struct;
%     in_bin.Totals=mRNA_nums.Totals(locs);
%     in_bin_txn.Totals=txn_nums.Totals(locs);
%     in_bin_txn.totalnums=txn_nums.totalnums(locs);
%     [x2,y2]=numbers_graph(in_bin,in_bin_txn);
%     subplot(2,ceil(numbins./2),i)
%     scatter(x2,y2,100,[0 0 0],'filled');
%     set(gca,'FontSize',12,'LineWidth',3,'FontName','Arial');
%     xlabel('mRNA/Cell');
%     ylabel('Txn Site Intensity');
% end