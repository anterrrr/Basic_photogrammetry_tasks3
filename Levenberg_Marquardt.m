function [numLoop, externalParameters] = Levenberg_Marquardt(pts3d, pts2d, externalParameter0, threshold, f)
%% 函数解释：高斯牛顿法最小化重投影误差
% pts3d:靶点空间坐标 pts2d:靶点像面坐标 externalParameter0：外方位参数初值  threshold：迭代停止条件

%% 参数初始化
phi0=externalParameter0(1); omega0=externalParameter0(2); k0=externalParameter0(3); tx0=externalParameter0(4); ty0=externalParameter0(5); tz0=externalParameter0(6);
I = eye(6); U=10; V=2;
%% 采用摄影测量中的先平移后旋转的坐标系运动关系定义共线方程
syms phi omega k tx ty tz Xw Yw Zw;%phi,omega,k分别对应y,x,z轴的旋转角度
syms fx(phi, omega, k, tx, ty,tz, Xw, Yw, Zw) fy(phi, omega, k, tx, ty,tz, Xw, Yw, Zw) R(phi, omega, k)
% 定义旋转矩阵
Rx = [1,0,0;0, cos(omega),sin(omega);0,-sin(omega),cos(omega)]; % x--omega
Ry = [cos(phi),0,-sin(phi);0,1,0;sin(phi),0, cos(phi)]; % y--phi
Rz = [cos(k),sin(k),0; -sin(k),cos(k),0;0,0,1];% z--k
R(phi, omega, k) = Rz*Rx*Ry;
temp = R(phi, omega, k);
r11 = temp(1,1);r12 = temp(1,2);r13 = temp(1,3);
r21 = temp(2,1);r22 = temp(2,2);r23 = temp(2,3);
r31 = temp(3,1);r32 = temp(3,2);r33 = temp(3,3);

fx(phi, omega, k, tx, ty, tz, Xw, Yw, Zw) = -f*(r11*(Xw-tx)+r12*(Yw-ty)+r13*(Zw-tz))/(r31*(Xw-tx)+r32*(Yw-ty)+r33*(Zw-tz));
fy(phi, omega, k, tx, ty, tz, Xw, Yw, Zw) = -f*(r21*(Xw-tx)+r22*(Yw-ty)+r23*(Zw-tz))/(r31*(Xw-tx)+r32*(Yw-ty)+r33*(Zw-tz));
% 对六个变量求导
fx_phi = diff(fx,phi);fx_omega = diff(fx,omega);fx_k = diff(fx,k);fx_tx = diff(fx,tx);fx_ty = diff(fx,ty);fx_tz = diff(fx,tz);
fy_phi = diff(fy,phi);fy_omega = diff(fy,omega);fy_k = diff(fy,k);fy_tx = diff(fy,tx);fy_ty = diff(fy,ty);fy_tz = diff(fy,tz);

%% 循环迭代
j=1;% while循环条件
numLoop=1;% 迭代次数
while j==1
    % 计算误差矩阵
    deltaF = [];
    for i=1:length(pts2d)
        temp1 = pts2d(i,1) - fx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3));
        temp2 = pts2d(i,2) - fy(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3));
        temp1 = double(temp1);
        temp2 = double(temp2);
        deltaF = [deltaF;temp1;temp2];
    end

    % 计算J雅可比矩阵
    J = [];
    for i=1:length(pts2d)
        temp1 = [double(fx_phi(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fx_omega(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fx_k(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fx_tx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fx_ty(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fx_tz(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3)))];
        temp2 = [double(fy_phi(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fy_omega(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fy_k(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fy_tx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fy_ty(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fy_tz(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3)))];
        J = [J;temp1;temp2];
    end
    
    % 计算外参数修正量
    deltaE = inv(J'*J+U*I)*J'*deltaF;
    deltaE = vpa(deltaE);
    
    % 初值计算像面误差
    f1 = [pts2d(1,1)-fx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(1,1), pts3d(1,2), pts3d(1,3));pts2d(1,2)-fy(phi0, omega0, k0, tx0, ty0, tz0, pts3d(1,1), pts3d(1,2), pts3d(1,3));
          pts2d(2,1)-fx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(2,1), pts3d(2,2), pts3d(2,3));pts2d(2,2)-fy(phi0, omega0, k0, tx0, ty0, tz0, pts3d(2,1), pts3d(2,2), pts3d(2,3))
          pts2d(3,1)-fx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(3,1), pts3d(3,2), pts3d(3,3));pts2d(3,2)-fy(phi0, omega0, k0, tx0, ty0, tz0, pts3d(3,1), pts3d(3,2), pts3d(3,3))
          pts2d(4,1)-fx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(4,1), pts3d(4,2), pts3d(4,3));pts2d(4,2)-fy(phi0, omega0, k0, tx0, ty0, tz0, pts3d(4,1), pts3d(4,2), pts3d(4,3))
          pts2d(5,1)-fx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(5,1), pts3d(5,2), pts3d(5,3));pts2d(5,2)-fy(phi0, omega0, k0, tx0, ty0, tz0, pts3d(5,1), pts3d(5,2), pts3d(5,3))
          pts2d(6,1)-fx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(6,1), pts3d(6,2), pts3d(6,3));pts2d(6,2)-fy(phi0, omega0, k0, tx0, ty0, tz0, pts3d(6,1), pts3d(6,2), pts3d(6,3))];
    f1 = vpa(f1);
    
    % 更新后的初值计算像面误差
    phi1 = phi0+deltaE(1);
    omega1 = omega0+deltaE(2);
    k1 = k0+deltaE(3);
    tx1 = tx0+deltaE(4);
    ty1 = ty0+deltaE(5);
    tz1 = tz0+deltaE(6);
    f2 = [pts2d(1,1)-fx(phi1, omega1, k1, tx1, ty1, tz1, pts3d(1,1), pts3d(1,2), pts3d(1,3));pts2d(1,2)-fy(phi1, omega1, k1, tx1, ty1, tz1, pts3d(1,1), pts3d(1,2), pts3d(1,3));
          pts2d(2,1)-fx(phi1, omega1, k1, tx1, ty1, tz1, pts3d(2,1), pts3d(2,2), pts3d(2,3));pts2d(2,2)-fy(phi1, omega1, k1, tx1, ty1, tz1, pts3d(2,1), pts3d(2,2), pts3d(2,3))
          pts2d(3,1)-fx(phi1, omega1, k1, tx1, ty1, tz1, pts3d(3,1), pts3d(3,2), pts3d(3,3));pts2d(3,2)-fy(phi1, omega1, k1, tx1, ty1, tz1, pts3d(3,1), pts3d(3,2), pts3d(3,3))
          pts2d(4,1)-fx(phi1, omega1, k1, tx1, ty1, tz1, pts3d(4,1), pts3d(4,2), pts3d(4,3));pts2d(4,2)-fy(phi1, omega1, k1, tx1, ty1, tz1, pts3d(4,1), pts3d(4,2), pts3d(4,3))
          pts2d(5,1)-fx(phi1, omega1, k1, tx1, ty1, tz1, pts3d(5,1), pts3d(5,2), pts3d(5,3));pts2d(5,2)-fy(phi1, omega1, k1, tx1, ty1, tz1, pts3d(5,1), pts3d(5,2), pts3d(5,3))
          pts2d(6,1)-fx(phi1, omega1, k1, tx1, ty1, tz1, pts3d(6,1), pts3d(6,2), pts3d(6,3));pts2d(6,2)-fy(phi1, omega1, k1, tx1, ty1, tz1, pts3d(6,1), pts3d(6,2), pts3d(6,3))];
    f2 = vpa(f2);
    
    % 更新U和V
    L = 0.5*deltaE'*(U*deltaE+J'*f1);
    L = vpa(L);
    rho = (f1'*f1-f2'*f2)/(2*L);
    rho = vpa(rho);
    if rho>0
        U = U*max(1/3, 1-(2*rho-1)^3);
        V = 2;
        phi0 = phi0+deltaE(1);
        omega0 = omega0+deltaE(2);
        k0 = k0+deltaE(3);
        tx0 = tx0+deltaE(4);
        ty0 = ty0+deltaE(5);
        tz0 = tz0+deltaE(6);
    else
        U = U*V;
        V = 2*V;
    end
    
    a = 1;   
    % 是否跳出循环
    if abs(deltaE(1))<threshold(1) && abs(deltaE(2))<threshold(2) && abs(deltaE(3))<threshold(3) && deltaE(4)<threshold(4) && deltaE(5)<threshold(5) && deltaE(6)<threshold(6)
        break;
    end

    % 更新迭代次数
    numLoop = numLoop+1;
end
externalParameters = [phi0;omega0;k0;tx0;ty0;tz0];
end