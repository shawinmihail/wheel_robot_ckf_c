clc
clear
close all
path = 'ckf/ckf/wheel_robot_log.csv';

data = csvread(path);
x = data(:,1);
y = data(:,2);
th = data(:,3);
x_dot = data(:,4);
y_dot = data(:,5);
th_dot = data(:,6);

v = (x_dot .^ 2 + y_dot .^ 2) .^ 0.5;

plot(v)