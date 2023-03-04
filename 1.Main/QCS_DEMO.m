% DEMO for quantile correlation-based Variable Screening
%
% Author:  Zihao Yuan
% Update:  2021/10/16
% Contact: zihaoyuan@whut.edu.cn


clc;
clearvars;
close all;

%% Example Case 1 in paper
% Rand Samples
N = 300;                        % sample number
P = 3000;                       % feature number
R = 4;                          % class number
q = 10;                         % number of active feature of every class
prm = 1;                        % method of sampling probability
n2w = @(num) num / sum(num); dn = floor(N / log(N));


% Set parameters for generation probability
switch prm
    case 1      % Generating samples with equal difference probability
        pr = n2w(1 : R);
    case 2      % Generating samples with equal probability
        pr = n2w(ones(1, R));
    case 3      % Generating samples with random probability
        pr = n2w(rand(1, R));
end

% Intial correlation parameters
mu = zeros(R, P);
for i = 1 : R
    mu(i, q * (i - 1) + (1 : q)) = 1.5;
end

% Random class by multinomial random numbers
[Y, ~] = ind2sub([R N], find(mnrnd(1, pr, N)' == 1));

a0 = 1;
% From Y to random X: where Xi = mu(r, :) + N(0, 1), which Yi belong to yr
X1 = 0.2 * randn(N, P) + 0.8 * randn(N, P) + 0.6 * mu(Y, :);
% From Y to random X: where Xi = mu(r, :) + U(0, 1), which Yi belong to yr
X2 = mu(Y, :) + a0 * rand(N, P);
% From Y to random X: Generated from a LPA model
X3 = 0.9 * X1 + 0.1 * trnd(1, N, P);
% Choose X from X1~X3
X = X3;

% Feature screening
[Fv1, Forder1] = QVS(X, Y);
[Fv2, Forder2, ~, Fvc2, Fvcorder2] = QCS_FDR(X, Y, 0.1);
