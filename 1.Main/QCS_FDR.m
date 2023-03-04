function [Vv, Vorder, Tau, Vvc, Vcorder, Tauc] = QCS_FDR(X, Y, alpha, Xtau, Ytau)
% Quantile Correlation-based Variable Selection
% 
% Input:
%   X      - Observation Samples(N * P matrix, N is numer of sample,
%           P is demension of feature)
%   Y      - Response variable(Y must be a N * 1 or 1 * N vector)
%   dn     - The number screening last or pre-specified FDR level
%   Xtau   - Grid points for X
%   Ytau   - Grid points for Y
%
% Output:
%   Vv     - Variable selection vector for response
%   Vorder - Sorting results of Variable(descend) for response
%   Tau    - Estimator of Tau in equation (2.2) in Tang et al.
%   Vvc    - Variable selection matrix for each grid
%   Vcorder- Sorting results of Variable(descend) for each grid
%   Tauc   - Estimator of  each grid of Tauc in equation (2.2) in Tang et al.
%
% Usage:
%   [Vv, Vorder, Tau, Vvc, Vcorder, Tauc] = QCS_FDR(X, Y, alpha, Xtau, Ytau)
% Reference:
%   Tang, W.; Xie, J.; Lin, Y.; Tang, N. Quantile Correlation-based
%   Variable Selection. Journal of Business & Economic Statistics 2021, pp.
%   641 1–13.

% Author  : Zihao Yuan
% Update  : 2022/04/07 (First Version: 2022/04/07)
% Contact : zihaoyuan@whut.edu.cn (If any suggestions or questions)


[N, P] = size(X);               % size of sample matrix
if ~exist('alpha','var') || isempty(alpha)
    alpha = max([0.05, 1/P]);   % given the defult alpha level
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
    XR = length(XCj);

    Tjst = zeros(XR, YR);
    EjstX = zeros(XR, 1);
    for jXRj = 1 : XR
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

[Tausort, Vorder] = sort(Tau, 'descend');   % sort abs_tao_hat by descend
FDR = P * (1 - chi2cdf(Tausort, (YR-1) * length(Xtau))) ./ (1 : P);
FDRoptindex = find(FDR <= alpha);
if isempty(FDRoptindex)
    Vv = [];
else
    Vv =  unique(Vorder(1 : FDRoptindex(end)));
end


if sum(abs(round(Y) - Y)) == 0
    Vcorder = zeros(YR, P);
    Vvc = cell(1, YR);
    for r = 1 : YR
        [Taucsort, Vcorder(r, :)] = sort(Tauc(r, :), 'descend');   % sort abs_tao_hat by descend
        FDRc = P * (1 - chi2cdf(Taucsort, length(Xtau))) ./ (1 : P);
        FDRcoptindex = find(FDRc <= alpha);
        if isempty(FDRcoptindex)
            Vvc{r} = [];
        else
            Vvc{r} = unique(Vcorder(r, 1 : FDRcoptindex(end)));
        end
    end
end

end



