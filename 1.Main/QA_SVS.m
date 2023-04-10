function [Fv, Fvs, upsilonk_hat, upsilon_hat] = QA_SVS(X, Y, dn)
% Quantile-Adaptive Sufficient Variable Screening by Controlling False Discovery
% 
% Input:
%   X      - Observation Samples(N * P matrix, N is numer of sample,
%               P is demension of feature)
%   Y      - Categorical responses(Y must be a N * 1 or 1 * N vector)
%   dn     - The number screening last or pre-specified FDR level
%
% Output:
%   Fv     - Feature selection result by every class
%   Fvs    - Sorting results of Feature(descend) by every class(K+1 * P)
%   upsilonk_hat
%          - Value of upsilonk_hat in equation (5)
%   upsilon_hat
%          - Value of upsilonk_hat in Corollary 1
%
% Usage:
%   [Fv, Fvs, upsilonk_hat, upsilon_hat] = QA_SVS(X, Y, dn)

% Author  : ZH.Yuan
% Update  : 2022/12/25 (First Version: 2021/10/14)
% Contact : zihaoyuan@whut.edu.cn (If any suggestions or questions)


[N, P] = size(X);           % size of sample matrix
if exist('dn', 'var') == 0 || isempty(dn)
    dn = min([round(N / log(N)) P]); % given the number screening last
end
if dn < 1
    alpha = dn;
    dn = 'FDR';
end

YC = sort(unique(Y));       % class of Y after sorted
K = length(YC);             % number of class

[~, Xp] = sort(X, 'descend');
[~, SX] = sort(Xp);

pk_hati = zeros(K, 1);
tauk_hat = zeros(K, P);
upsilonk_hat = zeros(K, P);

for k = 1 : K
    pk_hati(k) = sum(Y == YC(k)) / N;
    tauk_hat(k, :) = 1 / (N * (N+1) * pk_hati(k)) * sum(SX(Y == YC(k), :), 1) - 0.5;
    upsilonk_hat(k, :) = 12 * (N+1) * pk_hati(k) / (1 - pk_hati(k)) * tauk_hat(k, :).^2;
end

upsilon_hat = (K - 1) / K * sum(upsilonk_hat);
[upsilonk_hat_sort, Fvsk] = sort(upsilonk_hat, 2, 'descend');
[upsilon_hat_sort, Fvsall] = sort(upsilon_hat, 2, 'descend');
Fvs = [Fvsk; Fvsall];

switch dn

    case 'AFD'
        rhok = max([- norminv(1 / (2*P)), 3])^2;  % Filter operator
        rho = (K - 1) / K * rhok;

        % Screening the sufficient variables by adaptive false discovery
        Fv = cell(1, K+1);
        for k = 1 : K
            Fv{k} = find(upsilonk_hat(k, :) > rhok);
        end
        Fv{k+1} = find(upsilon_hat > rho);

    case 'FDR'
        % Screening the sufficient variables by false discovery rate
        Fv = cell(1, K+1);
        for k = 1 : K
            FDRk = P * (1 - chi2cdf(upsilonk_hat_sort(k, :), 1)) ./ (1 : P);
            FDRkoptindex = find(FDRk <= alpha);
            if ~isnan(FDRkoptindex)
                Fv{k} = Fvsk(k, 1 : FDRkoptindex(end));
            end
        end
        FDR = P * (1 - chi2cdf(upsilon_hat_sort, K - 1)) ./ (1 : P);
        FDRoptindex = find(FDR <= alpha);
        if ~isnan(FDRoptindex)
            Fv{k+1} = unique(Fvsall(1 : FDRoptindex(end)));
        end

    otherwise
        % Screening the sufficient variables by the dn top largest
        Fv = Fvs(:, 1 : dn);

end

