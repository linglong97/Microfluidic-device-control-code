%% 
function F=non_linear_eqns(stream_func,X)

N=100;


psi(1)=1; %fix Q1=1;
x=X(1);y=X(2); u=X(3); v=X(4);
psi(2)=stream_func(1); psi(3)=stream_func(2);

xm=ones(N,1)*x(:)';
ym=ones(N,1)*y(:)';
n=(1:N)'*ones(size(x));

c1=2*(psi(1)+(psi(2)-psi(1))*cos(n*pi))./sinh(n*pi);
c2=2*(psi(3)+(psi(2)-psi(3))*cos(n*pi))./sinh(n*pi);

eqn1= psi(2)+sum((c1.*sinh(n*pi.*(1-xm))+c2.*sinh(n*pi.*xm)).*cos(n*pi.*ym))-u;
eqn2=sum((c1.*cosh(n*pi.*(1-xm))-c2.*cosh(n*pi.*xm)).*sin(n*pi.*ym))-v;

% F(1)= psi(2)+sum((c1.*sinh(n*pi.*(1-xm))+c2.*sinh(n*pi.*xm)).*cos(n*pi.*ym))-u;
% F(2)=sum((c1.*cosh(n*pi.*(1-xm))-c2.*cosh(n*pi.*xm)).*sin(n*pi.*ym))-v;

 F=eqn1^2+eqn2^2;
end