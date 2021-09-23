close all   
clear   
echo on   
clc  
%读入数据定参数名
N=100;%训练次数
DATA=xlsread('DATAALLgai.xlsx',1,'B2:L169');             %读取所有H数据
J=floor(DATA(:,1)*0.1);  %进机场为原先0.1倍
C=floor(DATA(:,2)*0.135);%出机场为原先0.135倍
Z=zeros(168,1);
Z(1,1)=J(1)-C(1)
for i=1:167
Z((i+1),1)=J(i)-C(i)+Z(i,1);
end
X= (floor(DATA(:,6))*4);  %乘客为原来4倍
parameter=DATA(:,7:11);
P=[J,C,Z,X,parameter];% P 为输入矢量
T=[J,C,Z,X];% T 为目标矢量 
% 基于NARX时间序列，生成一个新的前向神经网络   
clc   
%  定义训练样本,8-13日      
Ptrain=P(1:143,:);   
Ptrain=Ptrain';
Ttrain=T(2:144,:);  
Ttrain=Ttrain';
clc   
for i=1:N
%  创建一个新的前向神经网络   
net=newff(Ptrain,Ttrain,[12],{'tansig','purelin'},'traingdm')     
%  设置训练参数    
net.trainParam.lr = 0.05;   
net.trainParam.mc = 0.9;   
net.trainParam.epochs = 1000;   
net.trainParam.goal = 1e-7;      
%  调用 TRAINGDM 算法训练 BP 网络   
[net,tr]=train(net,Ptrain,Ttrain);  
%  定义测试样本,14日 
Ptest=P(144:167,:);
Ptest=Ptest';
Ttest=T(145:168,:);
Ttest=Ttest';
%  对 BP 网络进行仿真   
Tpred = sim(net,Ptest) ;
for j=1:4
   for z=1:24
       if Tpred(j,z)<0
          Tpred(j,z)=0;
       end
   end
end
%  计算仿真误差   
 

Jpred(i,:)=floor(Tpred(1,:));
Cpred(i,:)=floor(Tpred(2,:));
Xpred(i,:)=floor(Tpred(4,:));
Zpred(i,:)=floor(Tpred(3,:));
E(((4*(i-1))+1):((4*(i-1))+4),:) = abs(Tpred - Ttest);  
clc   
echo off 
end
%数学模型计算结果
da1=17;
pa1=0.2069;
da2=37;
pa2=0.6294;
da3=40;
pa3=0.1637;
mxd1=49.0907;
mxd2=110.5573;
mxd3=120.7133;
mxn1=60.768;
mxn2=143.616;
mxn3=156.81;


%CASE1持续模型
for i=1:18
       xxx=Z((149+i))-X((149+i));
       n=0;
while xxx>0 
      n=n+1;
      xxx=xxx-(X(149+i)*n);
      if (i+n)>17
          Cost(i)=34.32;
         break
      end
end
     if ((Z(149+i))/(X(149+i)))>5  %超过5hour
         choose(i)=-100;
         Cost(i)=34.32;
     else
        choose(i)=pa1*(mxd1/(1.633+n))+pa2*(mxd2/(2.3+n))+pa3*(mxd3/(2.4+n))-34.32;
        Cost(i)=pa1*(mxd1/(1.633+n))+pa2*(mxd2/(2.3+n))+pa3*(mxd3/(2.4+n));
     end
end
%wanshang
for i=1:6
       xxx=Z((144+i))-X((144+i));
       n=0;
while xxx>0 
      n=n+1;
      xxx=xxx-(X(144+i)*n);
      if (i+n)>5
         Cost(i+18)=44.616; 
         break
      end
end
     if ((Z(144+i))/(X(144+i)))>5  %超过5小时
         choose(i+18)=-100;
         Cost(i+18)=44.616;
     else
        choose(i+18)=pa1*(mxn1/(1.633+n))+pa2*(mxn2/(2.3+n))+pa3*(mxn3/(2.4+n))-44.616;
        Cost(i+18)=pa1*(mxn1/(1.633+n))+pa2*(mxn2/(2.3+n))+pa3*(mxn3/(2.4+n));
     end
    
end
   Costaverage=sum(Cost)/24;
   
%CASE2本方案
for j=1:N
for i=1:18
       xxx=Z((149+i))-X((149+i));
       n=0;
while xxx>0 
      n=n+1;
      xxx=xxx-(Xpred(j,(4+i))*n);
      if (i+n)>17
          Costnew(i,j)=34.32;
         break
      end
end

    choosenew(i,j)=pa1*(mxd1/(1.633+n))+pa2*(mxd2/(2.3+n))+pa3*(mxd3/(2.4+n))-34.32;
    
    if choosenew(i,j)>0
        Costnew(i,j)=pa1*(mxd1/(1.633+n))+pa2*(mxd2/(2.3+n))+pa3*(mxd3/(2.4+n));
    else
        Costnew(i,j)=34.32;
    end
end
%wanshang
for i=1:6
       xxx=Z((144+i))-X((144+i));
       n=0;
while xxx>0 
      n=n+1;
      xxx=xxx-(Xpred(j,i)*n);
      if (i+n)>5
         Costnew(i+18,j)=44.616; 
         break
      end
end

    choosenew(i+18,j)=pa1*(mxn1/(1.633+n))+pa2*(mxn2/(2.3+n))+pa3*(mxn3/(2.4+n))-44.616;
   if choosenew(i+18,j)>0
        Costnew(i+18,j)=pa1*(mxn1/(1.633+n))+pa2*(mxn2/(2.3+n))+pa3*(mxn3/(2.4+n));
   else
       Costnew(i+18,j)=44.616; 
    end
    
end
   Costaveragenew(j)=sum(Costnew(:,j))/24;
end
Costaveragelast=sum(Costaveragenew)/N;
MSE=mse(E);  