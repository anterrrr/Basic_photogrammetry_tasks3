clc;clear;close all;
%% 初始化数据
pts3d=load('./data/3d position.txt');
pts2d=load('./data/2d position.txt');
intrinsicParameter=load('./data/intrinsic parameter.txt');
xp=intrinsicParameter(1);  yp=intrinsicParameter(2); f=intrinsicParameter(3);% 主点偏移和焦距

% 对像面靶点像面坐标进行畸变矫正(2d)
[row,~]=size(pts2d);
pts2d=pts2d-ones([row,1])*[xp,yp]+FunDistortionCorrect(pts2d,intrinsicParameter);
K = [f,0,0;0,f,0;0,0,1];  % 畸变矫正之后的内参数矩阵

% 初始化外方位参数初值
k0=pi/2; omega0=0; phi0=pi/2; tx0=2000; ty0=0; tz0=0;

% 定义迭代停止条件
threshold = [1e-4,1e-4,1e-4,1e-4,1e-4,1e-4];

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
n=1;% 迭代次数
while j==1
    % 计算误差矩阵
    delta_f = [];
    for i=1:length(pts2d)
        temp1 = pts2d(i,1) - fx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3));
        temp2 = pts2d(i,2) - fy(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3));
        temp1 = double(temp1);
        temp2 = double(temp2);
        delta_f = [delta_f;temp1;temp2];
    end

    % 计算J雅可比矩阵
    J = [];
    for i=1:length(pts2d)
        temp1 = [double(fx_phi(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fx_omega(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fx_k(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fx_tx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fx_ty(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fx_tz(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3)))];
        temp2 = [double(fy_phi(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fy_omega(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fy_k(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fy_tx(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fy_ty(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3))), double(fy_tz(phi0, omega0, k0, tx0, ty0, tz0, pts3d(i,1), pts3d(i,2), pts3d(i,3)))];
        J = [J;temp1;temp2];
    end
    
    % 计算外参数修正量
    delta_e = inv(J.'*J)*J.'*delta_f;
    
    % 更新初值
    phi0 = phi0+delta_e(1);
    omega0 = omega0+delta_e(2);
    k0 = k0+delta_e(3);
    tx0 = tx0+delta_e(4);
    ty0 = ty0+delta_e(5);
    tz0 = tz0+delta_e(6);

    % 是否跳出循环
    if abs(delta_e(1))<threshold(1) && abs(delta_e(2))<threshold(2) && abs(delta_e(3))<threshold(3) && delta_e(4)<threshold(4) && delta_e(5)<threshold(5) && delta_e(6)<threshold(6)
        break;
    end

    % 更新迭代次数
    n = n+1;
end

sprintf('Gauss_Newton算法最终迭代%d次', n)
sprintf('相机朝向的三个欧拉角分别为%.2f, %.2f, %.2f', phi0,omega0,k0)
sprintf('相机世界坐标系中的三坐标分别为%.2f, %.2f, %.2f', tx0,ty0,tz0)



