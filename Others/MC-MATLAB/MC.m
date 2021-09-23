clear
dt=1/24.0;                 
S0=34.42;                     
r=0.015;                      
sigma=1.0;                    
expTerm=r*dt;                 
stddev=sigma*sqrt(dt);      
nDays1=90;                   
for nDays=1:nDays1          
nTrials=10000;               
for j=1:nTrials 
n = randn(1,nDays);          
S=S0; 
for i=1:nDays 
dS = S*(expTerm+stddev*n(i)); 
S=S+dS;                     
end
S1(nDays,j)=S;               
end
end
S2=mean(S1');                
plot(S2','-o')                   
figure(2)
hist(S1(24,:),0:0.5:35)      