function [data]=Kon_sweep(min1,max1,Ncells)

Kons=linspace(min1,max1,30);
%for log sweep use logspace function
Koff=2;
Kini=20; %rate of transcription initialization units: min^-1

Telong=3.05;
Kdeg=.04;
k=0;
for jj=Kons
    Kon=jj;
    k=k+1;
        for i=1:Ncells %simulate 2 alleles for each cell mature and nascent
            [temp1]=two_state_model_mature(Kini,Kdeg,Kon,Koff);
            [temp2]=two_state_model_mature(Kini,Kdeg,Kon,Koff);
            [temp3]=two_state_model_nascent(Kini,Kon,Koff,Telong);
            [temp4]=two_state_model_nascent(Kini,Kon,Koff,Telong);
            rna(i,k)=temp1.RNA(end)+temp2.RNA(end);
            nascent(i,k)=temp3.nascent(end)*.51+temp4.nascent(end)*.51;
            %average number of visible nascent rnas: fraction of txn
            %length between first and last probe* .5 + fraction of txn
            %length between last probe and txn stop
        end
    freq(k)=nnz(nascent(:,k)>=2)./Ncells;
    medianrna(k)=median(rna(:,k));
    kon(k)=Kon;
    kini(k)=Kini;
end

data.nascent=nascent;
data.rna=rna;
data.medianRNA(:,1)=medianrna;
data.freq(:,1)=freq;
data.kon(:,1)=kon;
data.kini(:,1)=kini;