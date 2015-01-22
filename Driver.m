% Driver to call CSSA.m
% by Qian, 1/22/2015.

clear all;
close all;
clc;

k = 2; % the sparsity
b = rand(6, 5);

% Tree hierarchy
Tree = [0 0 0; 1 1 0; 1 1 0; 1 2 0; 2 2 3];


[x]=CSSA(b,k,Tree);