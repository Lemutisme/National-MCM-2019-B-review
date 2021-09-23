close all   
clear   
echo on   
clc  
%�������ݶ�������
N=100;%ѵ������
DATA=xlsread('DATAALLgai.xlsx',1,'B2:L169');             %��ȡ����H����
J=floor(DATA(:,1)*0.1);  %������Ϊԭ��0.1��
C=floor(DATA(:,2)*0.135);%������Ϊԭ��0.135��
Z=zeros(168,1);
Z(1,1)=J(1)-C(1)
for i=1:167
Z((i+1),1)=J(i)-C(i)+Z(i,1);
end
X= (floor(DATA(:,6))*4);  %�˿�Ϊԭ��4��
parameter=DATA(:,7:11);
P=[J,C,Z,X,parameter];% P Ϊ����ʸ��
T=[J,C,Z,X];% T ΪĿ��ʸ�� 
% ����NARXʱ�����У�����һ���µ�ǰ��������   
clc   
%  ����ѵ������,8-13��      
Ptrain=P(1:143,:);   
Ptrain=Ptrain';
Ttrain=T(2:144,:);  
Ttrain=Ttrain';
clc   
for i=1:N
%  ����һ���µ�ǰ��������   
net=newff(Ptrain,Ttrain,[12],{'tansig','purelin'},'traingdm')     
%  ����ѵ������    
net.trainParam.lr = 0.05;   
net.trainParam.mc = 0.9;   
net.trainParam.epochs = 1000;   
net.trainParam.goal = 1e-7;      
%  ���� TRAINGDM �㷨ѵ�� BP ����   
[net,tr]=train(net,Ptrain,Ttrain);  
%  �����������,14�� 
Ptest=P(144:167,:);
Ptest=Ptest';
Ttest=T(145:168,:);
Ttest=Ttest';
%  �� BP ������з���   
Tpred = sim(net,Ptest) ;
for j=1:4
   for z=1:24
       if Tpred(j,z)<0
          Tpred(j,z)=0;
       end
   end
end
%  ����������   
 

Jpred(i,:)=floor(Tpred(1,:));
Cpred(i,:)=floor(Tpred(2,:));
Xpred(i,:)=floor(Tpred(4,:));
Zpred(i,:)=floor(Tpred(3,:));
E(((4*(i-1))+1):((4*(i-1))+4),:) = abs(Tpred - Ttest);  
clc   
echo off 
end
%��ѧģ�ͼ�����
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


%CASE1����ģ��
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
     if ((Z(149+i))/(X(149+i)))>5  %����5hour
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
     if ((Z(144+i))/(X(144+i)))>5  %����5Сʱ
         choose(i+18)=-100;
         Cost(i+18)=44.616;
     else
        choose(i+18)=pa1*(mxn1/(1.633+n))+pa2*(mxn2/(2.3+n))+pa3*(mxn3/(2.4+n))-44.616;
        Cost(i+18)=pa1*(mxn1/(1.633+n))+pa2*(mxn2/(2.3+n))+pa3*(mxn3/(2.4+n));
     end
    
end
   Costaverage=sum(Cost)/24;
   
%CASE2������
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