clc; clear all; %clear console

%Example #1 using Newton's method on a single input, single valued function
syms x %create symbolic variable x
f = (x-1)^2 %we are going to find the zero of this symbolic function
x0 = 0 %initial value
eps = 1e-6 %tolerance
Newton_Raphson(x0, f, eps) %run Newton's method

%Example #2 using Newton's method on a multi input, multi valued function
syms x1 x2 %create symbolic variables [x1; x2]
f = [(x1+1)^2; (x2-1)^2] %we are going to find the zero of this 2D symbolic function
x0 = [0; 0] %initial value
eps = 1e-6 %tolerance
Newton_Raphson_multi(x0, f, eps) %run Newton's method