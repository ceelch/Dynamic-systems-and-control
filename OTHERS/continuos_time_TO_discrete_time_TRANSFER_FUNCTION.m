%% ESTE PROGRAMA TRANSFORMA UNA FUNCION DE TRANSFERENCIA EN S A 
% UNA FUNCION DE TRANSFERENCIA EN Z CON UN TIEMPO DE MUESTRE T


clear all;
clc;

num=[1.3]
den=[13.9 1]

sysc=tf(num,den)
sysd=c2d(sysc,0.71)

