function [Ad,Bd] = getAdBd(xe,ue,dt,m,g,J)
syms o1 o2 o3 t1 t2 t3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4 real
    c1 = cos(t1);
    s1 = sin(t1);
    c2 = cos(t2);
    s2 = sin(t2);
    c3 = cos(t3);
    s3 = sin(t3);
    odot = [v1;v2;v3];
    N = [-s2 0 1;
        c2*s3 c3 0;
        c2*c3 -s3 0];
    tdot = N\[w1;w2;w3];
    Rz = [c1 -s1 0;
          s1 c1 0; 
          0 0 1];
    Ry = [c2 0 s2; 
          0 1 0;
         -s2 0 c2];
    Rx = [1 0 0; 
         0 c3 -s3; 
         0 s3 c3];
    R = Rz*Ry*Rx;
    vdot = [0;0;g] - R*[0;0;u4/m];
    wdot = J\([u1; u2; u3]-[0 -w3 w2; w3 0 -w1; -w2 w1 0]*J*[w1;w2;w3]);
    h_0 = [odot; tdot; vdot; wdot];
    
    x = [o1; o2; o3; t1; t2;t3;v1;v2;v3;w1;w2;w3];
    u = [u1;u2;u3;u4];
    vars = [o1, o2, o3, t1, t2, t3, v1,v2,v3,w1,w2,w3,u1,u2,u3,u4];
    vals = [xe(1),xe(2),xe(3),xe(4),xe(5),xe(6),xe(7),xe(8),xe(9),xe(10),xe(11),xe(12),ue(1),ue(2),ue(3),ue(4)];
    Ac = double(subs(jacobian(h_0,x), vars, vals));
% for reference, this is the matrix Ac:
%     Ac = [0 0 0 0 0 0 1 0 0 0 0 0;
%           0 0 0 0 0 0 0 1 0 0 0 0;
%           0 0 0 0 0 0 0 0 1 0 0 0;
%           0 0 0 0 0 0 0 0 0 0 0 1;
%           0 0 0 0 0 0 0 0 0 0 1 0;
%           0 0 0 0 0 0 0 0 0 1 0 0;
%           0 0 0 0 -g*cos(xe(4)) -g*sin(xe(4)) 0 0 0 0 0 0;
%           0 0 0 0 -g*sin(xe(4)) g*cos(xe(4)) 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0 0 0 0]
    Bc = double(subs(jacobian(h_0,u), vars, vals));
    [Ad,Bd] = ssdata(c2d(ss(Ac,Bc,[],[]),dt));