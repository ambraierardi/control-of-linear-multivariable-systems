clc
clear
close all

%% ROV exercise
% underwater remotely operated vehicle

% define parameters of the system
alpha=
beta=
gamma=
m=
j=
theta=zeros(4,1);
J=zeros(4,1);

% define matrices
A=[-alpha/m 0 0; 0 -beta/m 0; 0 0 -gamma/j];
B=[cos(theta(1))/m cos(theta(2))/m cos(theta(3))/m cos(theta(4))/m; sin(theta(1))/m sin(theta(2))/m sin(theta(3))/m sin(theta(4))/m; J(1)/j J(2)/j J(3)/j] J(4)/j]
