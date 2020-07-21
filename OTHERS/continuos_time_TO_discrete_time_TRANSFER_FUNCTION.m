%% 
% THIS PROGRAM CONVERTS A "S" TRANSFER FUNCTION TO A "Z" TRANSFER 
% FUNCTION WITH A SAMPLING TIME T
%%
clear all;
clc;

num=[1.3]
den=[13.9 1]

sysc=tf(num,den)
sysd=c2d(sysc,0.71)