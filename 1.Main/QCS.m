function [Vv, Vorder, Tau, Vvc, Vcorder, Tauc] = QCS(X, Y, dn, Xtau, Ytau)
% Expectile correlation-based Variable Screening for
%   Ultra-High Dimensional Heterogeneous Categorical Data
% 
% Input:
%   X      - Observation Samples(N * P matrix, N is numer of sample,
%           P is demension of feature)
%   Y      - Response variable(Y must be a N * 1 or 1 * N vector)
%
% Output:
%   Vv     - Variable selection result
%   Vorder - Sorting results of Variable(descend)
%   Q_j    - Estimator of Q_j
%
% Usage:
%   [Vv, Vorder, Q_j] = ECS(X, Y, dn)

% Author  : ZH.Yuan
% Update  : 2022/04/07 (First Version: 2022/04/07)
% Email   : zihaoyuan@whut.edu.cn (If any suggestions or questions)


[N, P] = size(X);               % size of sample matrix
if ~exist('dn','var') || isempty(dn)
    dn = floor(N / log(N));     % given the number screening last
end
if ~exist('Xtau','var') || isempty(Xtau)
    tn = 6;
    Xtau = (1:tn) / (tn+1);     % default setting for tau
end
if ~exist('Ytau','var') || isempty(Ytau)
    tn = 6;
    Ytau = (1:tn) / (tn+1);     % default setting for tau
end

% Discrete Y
Y = reshape(Y(:), N, 1);
if sum(abs(round(Y) - Y)) == 0
    YD = Y;
else
    Yqx = quantile(Y, Ytau);    % expectile of Y
    YD = sum(Y >= Yqx, 2) + 1;  % class of Y after sorted
end
YC = sort(unique(YD));          % class of Y after sorted
YR = length(YC);                % number of class
EjstY = zeros(1, YR);
for jYR = 1 : YR
    EjstY(jYR) = sum(YD == YC(jYR));
end

Tau = zeros(1, P);
if sum(abs(round(Y) - Y)) == 0
    Tauc = zeros(YR, P);
end
parfor j = 1 : P

    if sum(abs(round(X(:, j)) - X(:, j))) == 0
        XDj = X(:, j);
    else
        XDj = sum(X(:, j) >= quantile(X(:, j), Xtau), 2) + 1;
    end
    XCj = sort(unique(XDj));
    XRj = length(XCj);

    Tjst = zeros(XRj, YR);
    EjstX = zeros(XRj, 1);
    for jXRj = 1 : XRj
        EjstX(jXRj) = sum(XDj == XCj(jXRj));
        for jYR = 1 : YR
            Tjst(jXRj, jYR) = sum(XDj == XCj(jXRj) & YD == YC(jYR));
        end
    end

    Ejst = EjstX * EjstY / N;

    Tau(j) = sum((Tjst - Ejst).^2 ./ Ejst, "all");
    if sum(abs(round(Y) - Y)) == 0
        Tauc(:, j) = sum((Tjst - Ejst).^2 ./ Ejst, 1);
    end

end

[~, Vorder] = sort(Tau, 'descend');      % sort abs_tao_hat by descend
Vv = Vorder(1 : dn);

if sum(abs(round(Y) - Y)) == 0
    Vcorder = zeros(YR, P);
    Vvc = zeros(YR, dn);
    for r = 1 : YR
        [~, Vcorder(r, :)] = sort(Tauc(r, :), 'descend');      % sort abs_tao_hat by descend
        Vvc(r, :) = Vcorder(r, 1 : dn);
    end
end

end

