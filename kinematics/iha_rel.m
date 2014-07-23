function[N_IHA,V_IHA,S_IHA]=iha_rel(Rv1,Rv2,Pv1,Pv2)
% function[N_IHA,V_IHA,S_IHA]=iha_rel(Rv1,Rv2,Pv1,Pv2)
%
% estimation of IHA en FHA parameters (plus some extra's)
% based on the work by Woltring
% Input is Rv (m * 9) elements of the orientation-matrix
%          Pv (m * 9) elements of the position vector
% Output is theta (rotation)
%                 Euler (x,y,z) decomposed using rotxyz(Rv)
%                 = -rotzyx(Rv') 
%                 N_IHA position vector IHA
%                 V_IHA translation vector IHA
%                 S_IHA point on the IHA
% see Woltring (1990) of Spoor (1980)
%
% HEJV, based on F-C SU, 1994

I=[1 0 0;0 1 0;0 0 1];

[m,n]=size(Rv1);

N_IHA=zeros(m,3);
V_IHA=zeros(m,3);
S_IHA=zeros(m,3);
omega=zeros(m,3);

Rv1dot=afgnew(Rv1,10);
Pv1dot=afgnew(Pv1,10);

Rv2dot=afgnew(Rv2,10);
Pv2dot=afgnew(Pv2,10);

for i=1:m,
   Rm1=reshape(Rv1(i,:),3,3)'; 
   Vm1=[reshape(Rv1dot(i,:),3,3)'];
   w=0.5*(Vm1*Rm1'-Rm1*Vm1');
   Wu=[w(3,2);w(1,3);w(2,1)];
   omega1(i,:)=Wu';

   Rm2=reshape(Rv2(i,:),3,3)'; 
   Vm2=[reshape(Rv2dot(i,:),3,3)'];
   w=0.5*(Vm2*Rm2'-Rm2*Vm2');
   Wu=[w(3,2);w(1,3);w(2,1)];
   omega2(i,:)=Wu';

  Vu21=[Pv2dot(i,:)-Pv1dot(i,:)]';
  Pu21=[Pv2(i,:)-Pv1(i,:)]';
  Wu21=[omega2(i,:)-omega1(i,:)]';

   if norm(Wu)>.25, % extra limitation of a minimum rotation velocity  (.25 rad/sec)
      wu=sqrt(Wu21'*Wu21);
      Nu=Wu21/wu;
      N_IHA(i,:)=Nu';
      vu=Vu21'*Nu;
      Su=Pu21+vecprod(Wu21,Vu21/(wu*wu))+Pv1(i,:)';
      S_IHA(i,:)=Su';
      V_IHA(i)=vu;
   end
end
