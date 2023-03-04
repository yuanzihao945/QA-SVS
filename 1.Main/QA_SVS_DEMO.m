% DEMO for Quantile-Adaptive Sufficient Variable Screening
%
% Author:  Zihao Yuan
% Update:  2022/12/30
% Contact: zihaoyuan@whut.edu.cn

clearvars
close all
clc

fun1 = @(x) 0.5 * (x(:, 1) + x(:, 2) + x(:, 101)) + randn([size(x, 1), 1]);
K = 4;
QP = (1 : (K - 1)) / K;
X = randn(500, 5000);
Y = fun(X);
Z = sum(Y >= quantile(Y, QP), 2) + 1;
[Fv0, Fvs0] = QA_SIS(X, Z);
[Fv1, Fvs1] = sort(QA_SIS(X, Z, 0.05), 'descend');
A_num = find(ismember(Fvs1, [1 2 101]), 1, 'last');

