clc
clear
close all
path = 'act_state_log.csv';

data = csvread(path);

x = data(:,1);
y = data(:,2);
z = data(:,3);

vx = data(:,4);
vy = data(:,5);
vz = data(:,6);

ax = data(:,7);
ay = data(:,8);
az = data(:,9);

qx = data(:,10);
qy = data(:,11);
qz = data(:,12);

wx = data(:,13);
wy = data(:,14);
wz = data(:,15);

% v = (x_dot .^ 2 + y_dot .^ 2) .^ 0.5;

plot(x,y)