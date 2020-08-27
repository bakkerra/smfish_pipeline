function [data]=Kini_sweep(min1,max1,Ncells)
%varies Kini between two values min1 and max1
%Ncells is number of cells in simulation



Kinis=linspace(min1,max1,30);%rate of transcription initialization units: min^-1
%for log sweep use logspace function 
Koff=2; %rate of state change from on>off units: min^-1
Kon=.04; %rate of state change from off>on units: min^-1
Telong=3.05;  %transcription time from 1st probe to stop site units: min
Kdeg=.04; %rate of rna degradation units min^-1
k=0; %initialize k

for jj=Kinis
    Kini=jj;
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
