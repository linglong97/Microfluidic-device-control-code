function [P_new] = eqns_for_flowratepressure( P, Q )

 P0= (P(1)+P(2)+P(3)+P(4))/4;
 
 eqn1= P(1)-Q(1)-P0;
 eqn2= P(2)-Q(2)-P0;
 eqn3=P(3)-Q(3)-P0;

 eqn4=P(4)-Q(4)-P0;

%  P_new=[eqn1;eqn2;eqn3;eqn4;eqn5;eqn6];

P_new=sqrt(eqn1^2+eqn2^2+eqn3^2+eqn4^2);
 
end

