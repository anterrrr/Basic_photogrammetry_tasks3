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
phi0=pi/2; omega0=0; k0=pi/2; tx0=2000; ty0=0; tz0=0;
externalParameter0 = [phi0;omega0;k0;tx0;ty0;tz0];
% 定义迭代停止条件
threshold = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];

% 高斯牛顿法进行位姿估计
[numLoop1, externalParameters1] = Gauss_Newton(pts3d, pts2d, externalParameter0, threshold, f);
disp(externalParameters1)

% 隆贝格马奎特法进行位姿估计
[numLoop2, externalParameters2] = Levenberg_Marquardt(pts3d, pts2d, externalParameter0, threshold, f);
disp(externalParameters2)


