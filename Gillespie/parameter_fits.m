function [data]=parameter_fits(min1,max1,min2,max2,Ncells)
%simulates combinations of burst size and frequency
%between 10^min and 10^max
Kinis=logspace(min1,max1,10); %txn start rate 
Kons=logspace(min2,max2,10); %state transition rate off>on

%all parameters are in minutes or minutes^-1
% Kon
Koff=1; %state transition rate on>off
Telong=3.05; %estimated txn time from 1st probe to txn stop
vis=.67;%avg fraction of total actual nascent rnas visible in fish data
%calculate length between probe 1 and txn stop
%vis= fraction of this length between probe 1 and last probe *.5 
% +fraction of length between last probe and txn stop
Kdeg=.04; %based on half life of 15 min for salm and sens
k=0;

for jj=Kinis
    Kini=jj;
    for j=Kons
        Kon=j;
        k=k+1;
        for i=1:Ncells %simulate 2 alleles for each cell mature and nascent
            [temp1]=two_state_model_mature(Kini,Kdeg,Kon,Koff);
            [temp2]=two_state_model_mature(Kini,Kdeg,Kon,Koff);
            [temp3]=two_state_model_nascent(Kini,Kon,Koff,Telong);
            [temp4]=two_state_model_nascent(Kini,Kon,Koff,Telong);
            rna(i,k)=temp1.RNA(end)+temp2.RNA(end);
            nascent(i,k)=temp3.nascent(end)*vis+temp4.nascent(end)*vis;
        end
        freq(k)=nnz(nascent(:,k)>=2)./Ncells;
        meanrna(k)=mean(rna(:,k));
        meannascent(k)=mean(find(nascent(:,k)>=2));
        kon(k)=Kon;
        kini(k)=Kini;
    end
end

data.nascent=nascent;
data.rna=rna;
data.meanRNA=meanrna;
data.meanNascent=meannascent;
data.freq=freq;
data.kon=kon;
data.kini=kini;


