function V = corrected_voltages(x,y,K,P_E,P_W,P_N,P_S)
V3=P_N;V1=P_S;
% if (x<=0)
%     V1=P_N;
%     V2=P_S+K*abs(x);
%     
% elseif (x>0)
%     V1=P_N+K*abs(x);
%     V2=P_S;
% end

if (y>=0)
    V2=P_W;
    V4=P_E+K*abs(y);
    
elseif (y<0)
    V2=P_W+K*abs(y);
    V4=P_E;
end

if (V2>=P_N)
    V2=0.75*P_N;
    
elseif (V4>=P_N)
    V4=0.75*P_N;
end
V=[V1 V2 V3 V4];