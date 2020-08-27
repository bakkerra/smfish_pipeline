function [data]=two_state_model_nascent(Kini,Kon,Koff,Telong)
%simulates two state model of transcription for only the period of time that it takes to
%transcribe an rna

%simulation parameters
%sim_time=1000;
sim_time=Telong;%how long in minutes
data=struct;
%model parameters

%initialize at time=0 we have 1 DNA, and no nascent RNA, no nascent transcripts
d=1; %dna copies of gene
n=0; %nascent rnas
state=0;% 0 is promoter off 1 is promoter on
time=0;
counter=1;
timedata(counter)=time;
nascentdata(counter)=n;
                                                                                                                                                                                                                                                                                                                           
while time < sim_time
    if state==1
        sum(1)=0;
        sum(2)=0;
        sum(3)=0;
        a1=Kini*d;%likelihood of transcription initiation
        a2=Koff;%likelihood of state change
    
        a0=a1+a2;
        r1=rand(1,1);
        r2=rand(1,1);
        tau=-(1./a0)*log(r1);
        yr2=r2*a0;
        sum(2)=sum(1)+a1;
        sum(3)=sum(2)+a2;
        for ii=2:3
            if (sum(ii) >= yr2) && (sum(ii-1) < yr2) 
            mu=ii-1;
            end 
        end
        if (mu == 1)% && time>sim_time-deltat
            n=n+1;
        end
        if (mu == 2)
            state=0;
        end
      time=time + tau;
      counter=counter+1;
      timedata(counter)=time;
      nascentdata(counter)=n;
    end
    if state==0
        sum(1)=0;
        sum(2)=0;
        a1=Kon;%likelihood of state change, no initiation can occur         
        a0=a1;
        r1=rand(1,1);
        tau=-(1./a0)*log(r1);
        state=1;
        time=time + tau;
        counter=counter+1;
        timedata(counter)=time;
        nascentdata(counter)=n;
    end
end
data.time=timedata;
data.nascent=nascentdata;