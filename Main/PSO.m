clc
%function [xm,fv]=PSO(fitness,N,c1,c2,w,M,D)
%待优化的目标函数：fitness
%粒子数目：N
%学习因子1：c1
%学习因子2：c2
%惯性权重：w
N=8; % 粒子数
Iteration=8;%迭代次数
c1=2;
c2=2;
M=1;%问题维数  
w_max=0.9;
w_min=0.4;
D=1;  %1个N
PL1=zeros(N,1);
PL=[PL1];
f1=zeros(2*N,1);
f=[f1];

%初始化优化量
for i=1:N

    NNN(i,1)=floor(rand(1,1).*59)+1;

end

%
 x_font=zeros(N,1);
 x_P=[NNN, x_font];
 v1=rand(N,1).*12;
 v=[v1];
 y=zeros(N,1);

 %pg=zeros(N,1);
 
           for i=1:N
for j=1:24               
LONG(j)=C(144+j)*(pa2*mxn2/(2.3+pa1*20)+pa3*mxn3/(2.4+pa1*20));
for o=1:60
ZZZ=o*NNN(i,1)*Ptc(o);
end
SHORT(j)=C(144+j)*pa1*mxn1/(1.633+17/35+ZZZ)
end
PL1(i)=abs(sum(LONG)-sum(SHORT))/24; %COSTaverage


           end
               
        PL=[PL1];
        pb=x_P;
        n=find(PL1==min(PL1));
        pg=x_P(n,:);   
        
        f1=[PL1];
        f2=[PL2];
                   
          for t=1:Iteration
             x_Q=x_P;
            for i=1:N
                w=w_max-(t-1)*(w_max-w_min)/max((Iteration-1),1); %线性递减权重法
                v(i,:)=w*v(i,:)+c1*rand*(pb(i,1:D)-x_Q(i,1:D))+c2*rand*(pg(1:D)-x_Q(i,1:D));
                %速度越界处理
                
                     if v(i,1)<=-12
                         v(i,1)=-12;
                     end
                     if v(i,1)>=12
                         v(i,1)=12;
                     end
                
                x_Q(i,1:D)=x_Q(i,1:D)+v(i,1:D);
             %位置越界处理
             
                  if x_Q(i,1)<=2
                          x_Q(i,1)=2;
                  end
                  if x_Q(i,1)>=12
                          x_Q(i,1)=12;
                  end
             x_Q(i,1)=floor(x_Q(i,1));

             x=[x_Q;x_P];
             
             k=i+N;
             
f1(k)=f1(i); %COST

f2(k)=f2(i); %voltage fluctuation+frequency fluctuation

K=i;

for j=1:24               
LONG(j)=C(144+j)*(pa2*mxn2/(2.3+pa1*20)+pa3*mxn3/(2.4+pa1*20));
for o=1:60
ZZZ=o*x(K,1)*Ptc(o);
end
SHORT(j)=C(144+j)*pa1*mxn1/(1.633+17/35+ZZZ)
end
f1(K)=abs(sum(LONG)-sum(SHORT))/24; %COSTaverage






             f=[f1];
         
      
             
            end
           %排序
           front = 1;
           F(front).f = [];
           individual = [];
           x(:,D+1)=0;
           for i = 1 : 2*N
          % Number of individuals that dominate this individual
              individual(i).n = 0;
          % Individuals which this individual dominate
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
% Find the subsequent fronts
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
        
        % 选择Gbest
        % pg=sorted_based_on_front(1,1:D);
        current_index = 0;
        %选择下一代粒子从经过排序后的种群R中按顺序选择N个粒子到种群P
        for i=1:N
             x_P(i,:)=sorted_based_on_front(i,:);
        end
       pg=sorted_based_on_front(1,:);
         for i=1:N
             
     for j=1:24               
LONG(j)=C(144+j)*(pa2*mxn2/(2.3+pa1*20)+pa3*mxn3/(2.4+pa1*20));
for o=1:60
ZZZ=o*x_P(i,1)*Ptc(o);
end
SHORT(j)=C(144+j)*pa1*mxn1/(1.633+17/35+ZZZ)
end
f1(i)=abs(sum(LONG)-sum(SHORT))/24; %COSTaverage



                  end

        f=[f1];
        for i=1:N
          if  f1(i)<=PL1(i)  
                  if  f1(i)<=PL1(i)
                      PL1(i)=f1(i) ;   
                  end

                 pb(i,:)=x_P(i,:);
           end
        end
        
          end

                 n=find(PL1==min(PL1));
                 pg=x_P(n,:);
for j=1:24                 
LONG(j)=C(144+j)*(pa2*mxn2/(2.3+pa1*20)+pa3*mxn3/(2.4+pa1*20));
for o=1:60
ZZZ=o*pg(1,1)*Ptc(o);
end
SHORT(j)=C(144+j)*pa1*mxn1/(1.633+17/35+ZZZ)
end
FFF=abs(sum(LONG)-sum(SHORT))/24; %COSTaverage                 
   

                 
                 


                 
               