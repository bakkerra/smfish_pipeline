function [data]=two_state_model_mature(Kini,Kdeg,Kon,Koff)
%simulates two state model of transcription 

%simulation parameters
sim_time=1000; %how long is simulation in minutes
data=struct;
%model parameters

%initialize at time=0 we have 1 DNA, and no nascent RNA, no transcripts
d=1; %dna copies of gene
n=0; %nascent rnas
r=0; %transcribed rnas
state=0;%0 is promoter off 1 is prmoter on
time=0;
counter=1;
timedata(counter)=time;
RNAdata(counter)=r;                                                                                                                                                                                                                                                                                                                            
while time < sim_time
    if state==1
        sum(1)=0;
        sum(2)=0;
        sum(3)=0;
        sum(4)=0;
        a1=Kini*d;%likelihood of transcription initiation
        a2=Kdeg*r;%likelihood of degradation
        a3=Koff;%likelihood of state change
    
        a0=a1+a2+a3; %liklihood that any event occurs
        r1=rand(1,1);
        r2=rand(1,1);
        tau=-(1./a0)*log(r1); %random time step based on a0, how long until an event occurs
        yr2=r2*a0;%random select which event has occured based on k values
        sum(2)=sum(1)+a1;
        sum(3)=sum(2)+a2;
        sum(4)=sum(3)+a3; 
        for ii=2:4
            if (sum(ii) >= yr2) && (sum(ii-1) < yr2) 
            mu=ii-1;
            end 
        end
        if (mu == 1) 
            r=r+1;
        end
        if (mu == 2)
            r=r-1;
        end  
        if (mu == 3) 
            state=0;
        end
      time=time + tau;
      counter=counter+1;
      timedata(counter)=time;
      RNAdata(counter)=r;
    end
    if state==0
        sum(1)=0;
        sum(2)=0;
        a1=Kon;  %likelihood of state change
        a2=Kdeg*r; %likelihood of degradation
        
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
        if (mu == 1) 
            state=1;
        end
        if (mu == 2)
            r=r-1;
        end  
        time=time + tau;
        counter=counter+1;
        timedata(counter)=time;
        RNAdata(counter)=r;
    end
end
data.time=timedata;
data.RNA=RNAdata;