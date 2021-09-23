clc
N=2; 
Iteration=2;
c1=2;
c2=2;
M=2; 
w_max=0.9;
w_min=0.4;
D=1;  
PL1=zeros(N,1);
PL2=zeros(N,1);
PL=[PL1,PL2];
f1=zeros(2*N,1);
f2=zeros(2*N,1);
f=[f1,f2];
 
Zh=zeros(N,24);
R=zeros(N,1);

 
    a=zeros(24,1);
    b=[(X(144:167)-Z(144:167)),a];
ZL(1:24,1)=2*(X(145:168)+max(b,[],2));
 
for i=1:59
tc(i)=i;
Ptc(i)=exp((-1)*i)+0.00339;
end
tc(60)=60;
Ptc(60)=1-sum(Ptc(1:59));
 

for i=1:N
 
    R(i,1)=floor((rand(1,1).*(12-2))/2)*2+2;
    
 
      if ZL(1,1)>0
     Zh(i,1)=ZL(1,1)/2;
      else
        Zh(i,1)=R(i,1); 
      end
     c=[R(i,1),Zh(i,1)];
     MIN(i)=min(c);  
 
end
 

 x_font=zeros(N,1);
 x_P=[R, x_font];
 v1=rand(N,1).*2.5;
 v=[v1];
 y=zeros(N,1);
 

 
           for i=1:N
               
for j=1:24
    customer(i,j)=0;
    for k=1:60
        for l=1:(floor(60/k))
    customer(i,j)=customer(i,j)+k*Ptc(k)*l*MIN(i)
        end
    end
    costaverage(i,j)=C(144+j)*Costnew(j,1)/((49*pa1+69*pa2+72*pa3)/35+customer(i,j)/60)    
end
 
PL1(i)=-(sum(sum(costaverage(i,:)))/24); 
 
 
PL2(i)=sum(sum(customer(i,:)))/24;
 
           end
               
        PL=[PL1,PL2];
        pb=x_P;
        n=find(PL1==min(PL1));
        pg=x_P(n,:);   
        
        f1=[PL1;PL1];
        f2=[PL2;PL2];
                   
          for t=1:Iteration
             x_Q=x_P;
            for i=1:N
                w=w_max-(t-1)*(w_max-w_min)/max((Iteration-1),1);
                v(i,:)=w*v(i,:)+c1*rand*(pb(i,1:D)-x_Q(i,1:D))+c2*rand*(pg(1:D)-x_Q(i,1:D));
              
                     if v(i,1)<=-2.5
                         v(i,1)=-2.5;
                     end
                     if v(i,1)>=2.5
                         v(i,1)=2.5;
                     end
                
                x_Q(i,1:D)=x_Q(i,1:D)+v(i,1:D);
       
             
                  if x_Q(i,1)<=2
                          x_Q(i,1)=2;
                  end
                  if x_Q(i,1)>=12
                          x_Q(i,1)=12;
                  end
             x_Q(i,1)=floor(x_Q(i,1)/2)*2;
 
             x=[x_Q;x_P];
             
             k=i+N;
             
f1(k)=f1(i); 
 
f2(k)=f2(i); 
 
K=i;
 
      if ZL(1,1)>0
     Zh(K,1)=ZL(1,1)/2;
      else
        Zh(K,1)=x(K,1); 
      end
     c=[x(K,1),Zh(K,1)];
     MIN(K)=min(c);  
 
for j=1:24
    customer(i,j)=0;
    for k=1:60
        for l=1:(floor(60/k))
    customer(i,j)=customer(i,j)+k*Ptc(k)*l*MIN(K)
        end
    end
    costaverage(i,j)=C(144+j)*Costnew(j,1)/((49*pa1+69*pa2+72*pa3)/35+customer(i,j)/60)    
end
 
f1(K)=-(sum(sum(costaverage(i,:)))/24);
 
 
f2(K)=sum(sum(customer(i,:)))/24; 
 
 
 
 
             f=[f1,f2];
         
      
             
            end
          
           front = 1;
           F(front).f = [];
           individual = [];
           x(:,D+1)=0;
           for i = 1 : 2*N
       
              individual(i).n = 0;
         
              individual(i).p = [];
              for j = 1 : 2*N
                    dom_less = 0;
                    dom_equal = 0;
                    dom_more = 0;
                  for k = 1 : M
                       if (f(i,k) < f(j,k))
                              dom_less = dom_less + 1;
                       elseif (f(i,k) == f(j,k))
                              dom_equal = dom_equal + 1;
                       else
                              dom_more = dom_more + 1;
                       end
                  end
                  if dom_less == 0 && dom_equal ~= M
                          individual(i).n = individual(i).n + 1;
                  elseif dom_more == 0 && dom_equal ~= M
                          individual(i).p = [individual(i).p j];
                  end
               end   
                 if individual(i).n == 0
                      x(i,D+1) = 1;
                      F(front).f = [F(front).f i];
                 end
            end

        while ~isempty(F(front).f)
              Q = [];
          for i = 1 : length(F(front).f)
                if ~isempty(individual(F(front).f(i)).p)
                    for j = 1 : length(individual(F(front).f(i)).p)
                        individual(individual(F(front).f(i)).p(j)).n = ...
                            individual(individual(F(front).f(i)).p(j)).n - 1;
                        if individual(individual(F(front).f(i)).p(j)).n == 0
                            x(individual(F(front).f(i)).p(j),D+1) = ...
                                    front + 1;
                            Q = [Q individual(F(front).f(i)).p(j)];
                         end
                    end
                end
          end
             front =  front + 1;
             F(front).f = Q;
        end
        [temp,index_of_fronts] = sort(x(:,D+1));
        for i = 1 : length(index_of_fronts)
           sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
        end
        
       
        current_index = 0;
       
        for i=1:N
             x_P(i,:)=sorted_based_on_front(i,:);
        end
       pg=sorted_based_on_front(1,:);
         for i=1:N
             
      if ZL(1,1)>0
     Zh(i,1)=ZL(1,1)/2;
      else
        Zh(i,1)=x_P(i,1); 
      end
     c=[x_P(i,1),Zh(i,1)];
     MIN(i)=min(c);  
 
for j=1:24
    customer(i,j)=0;
    for k=1:60
        for l=1:(floor(60/k))
    customer(i,j)=customer(i,j)+k*Ptc(k)*l*MIN(i)
        end
    end
    costaverage(i,j)=C(144+j)*Costnew(j,1)/((49*pa1+69*pa2+72*pa3)/35+customer(i,j)/60)    
end
 
f1(i)=-(sum(sum(costaverage(i,:)))/24); 
 
 
f2(i)=sum(sum(customer(i,:)))/24;
 
 
                  end
 
        f=[f1,f2];
        for i=1:N
          if  f1(i)<=PL1(i) | f2(i)<=PL2(i)  
                  if  f1(i)<=PL1(i)
                      PL1(i)=f1(i) ;   
                  end
                  if  f2(i)<=PL2(i)
                      PL2(i)=f2(i) ;   
                  end
                 pb(i,:)=x_P(i,:);
           end
        end
        
          end
 
                 n=find(PL1==min(PL1));
                 pg=x_P(n,:);
                 
      if ZL(1,1)>0
     Zh(n(1),1)=ZL(1,1)/2;
      else
        Zh(n(1),1)=pg(n(1),1); 
      end
     d=[pg(1,1),Zh(n(1),1)];
     MIN(1)=min(d);  
 
for j=1:24
    customer(n,j)=0;
    for k=1:60
        for l=1:(floor(60/k))
    customer(n,j)=customer(n,j)+k*Ptc(k)*l*MIN(1)
        end
    end
    costaverage(n(1),j)=C(144+j)*Costnew(j,1)/((49*pa1+69*pa2+72*pa3)/35+customer(n,j)/60)    
end
 
COSTaverage=abs(-(sum(sum(costaverage(n(1),:)))/24)); 
 
 
Customer=sum(sum(customer(n,:)))/24; 

                 
                 


                 
               