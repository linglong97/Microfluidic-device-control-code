% function [U]=uvvalues_al(X,psi)
function [U]=uvvalues_al(X,psi)
% global N
N=100;
x=X(1);y=X(2);
xm=ones(N,1)*x(:)';
ym=ones(N,1)*y(:)';
n=(1:N)'*ones(size(x));
c1=2*(psi(1)+(psi(2)-psi(1))*cos(n*pi))./sinh(n*pi);
c2=2*(psi(3)+(psi(2)-psi(3))*cos(n*pi))./sinh(n*pi);

u=psi(2)+sum((c1.*sinh(n*pi.*(1-xm))+c2.*sinh(n*pi.*xm)).*cos(n*pi.*ym));
v=sum((c1.*cosh(n*pi.*(1-xm))-c2.*cosh(n*pi.*xm)).*sin(n*pi.*ym));

U= sqrt(u.^2+v.^2);
% pause
% U
end