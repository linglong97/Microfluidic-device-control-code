%function [V1,V2,V3,V4]=stag_point_control(xs,ys,X,K,V1,V2,V3,V4,P_MAX,P_min)  
global time Q2byQ1stor Q3byQ1stor xsstor ysstor center dim


% 
% X=[0.6,0.25];
% V1=1.5;V2=0.2;V3=1.5;V4=0.2;
% K=0.13/100;

x0=[X(1),X(2)]';
P_C=(V1+V2+V3+V4)/4;
Q1=(V2-P_C); Q2=(V4-P_C); Q3 = (V1-P_C);
[dudx,dudy]=graduvalues(X(1),X(2),[1 1+Q2/Q1 1+Q2/Q1+Q3/Q1]);
A=[dudx,dudy;dudy,-dudx];
maggradu=sqrt(dudx^2+dudy^2);
lmb1=-maggradu;
lmb2=maggradu;
mage1=sqrt((maggradu-dudx).^2+dudy^2);
mage2=sqrt((maggradu+dudx).^2+dudy^2);
e11=(maggradu-dudx)/mage1;
e12=-dudy/mage1;
e21=(maggradu+dudx)/mage2;
e22=dudy/mage2;

%t=toc(time);
t=1e-2;

B=[e11^2*exp(lmb1*t)+e21^2*exp(lmb2*t),e11*e12*exp(lmb1*t)+e21*e22*exp(lmb2*t);e11*e12*exp(lmb1*t)+e21*e22*exp(lmb2*t),e12^2*exp(lmb1*t)+e22^2*exp(lmb2*t)];

new_stag_point=(-A*B-K*B+K*eye(2))\((-A*B-K*B)*x0+K*[0.5;0.5]);

%%

Q2byQ1=interp2(xsstor,ysstor,Q2byQ1stor',new_stag_point(1,1),new_stag_point(2,1)); 
Q3byQ1=interp2(xsstor,ysstor,Q3byQ1stor',new_stag_point(1,1),new_stag_point(2,1));
   
if ((new_stag_point(2)>new_stag_point(1))&&((new_stag_point(2)+new_stag_point(1))<1))
    A=[0 -Q2byQ1 0 1 (Q2byQ1-1); 1 -Q3byQ1 0 0 (Q3byQ1-1); 1 1 1 1 -4; 0 0 1 0 0; 1 0 0 0 0];
    B=[0 0 0 P_min P_MAX]';
    [V]=A\B;
    V1=V(1);V2=V(2);V3=V(3);V4=V(4);
    
end

if ((new_stag_point(2)>new_stag_point(1))&&((new_stag_point(2)+new_stag_point(1))>1))
    A=[0 -Q2byQ1 0 1 (Q2byQ1-1); 1 -Q3byQ1 0 0 (Q3byQ1-1); 1 1 1 1 -4; 0 0 1 0 0; 0 1 0 0 0];
    B=[0 0 0 P_min P_MAX]';
    [V]=A\B;
    V1=V(1);V2=V(2);V3=V(3);V4=V(4);
    
end

if ((new_stag_point(2)<new_stag_point(1))&&((new_stag_point(2)+new_stag_point(1))<1))
    A=[0 -Q2byQ1 0 1 (Q2byQ1-1); 1 -Q3byQ1 0 0 (Q3byQ1-1); 1 1 1 1 -4; 0 0 0 1 0; 1 0 0 0 0];
    B=[0 0 0 P_min P_MAX]';
    [V]=A\B;
    V1=V(1);V2=V(2);V3=V(3);V4=V(4);
end

if ((new_stag_point(2)<new_stag_point(1))&&((new_stag_point(2)+new_stag_point(1))>1))
    A=[0 -Q2byQ1 0 1 (Q2byQ1-1); 1 -Q3byQ1 0 0 (Q3byQ1-1); 1 1 1 1 -4; 0 0 0 1 0; 0 1 0 0 0];
    B=[0 0 0 P_min P_MAX]';
    [V]=A\B;
    V1=V(1);V2=V(2);V3=V(3);V4=V(4);
%end