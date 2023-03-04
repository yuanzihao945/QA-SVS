clearvars
close all
clc

fun1 = @(x) 0.5 * (x(:, 1) + x(:, 2) + x(:, 101)) + randn([size(x, 1), 1]);
fun2 = @(x) 0.8 * x(:, 3) + 0.5 * (x(:, 4) + 1).^2 + tan(pi * (x(:, 102) + 1) / 4 )  + randn([size(x, 1), 1]);
fun3 = @(x) 0.5 * exp(x(:, 5)) + sin(pi * x(:, 6) / 2) ...
    + x(:, 103) .* (x(:, 103) > quantile(x(:, 103), 0.6)) + randn([size(x, 1), 1]);
fun4 = @(x) (1 + x(:, 7) + x(:, 8)).^(-3) .* exp(1 + 3 * sin(pi * x(:, 104) / 2)) + randn([size(x, 1), 1]);
fun5 = @(x) 0.5 * x(:, 9) + tan(pi * ((x(:, 10) - 1) .* (x(:, 105) + 1) / 4)) + randn([size(x, 1), 1]);
fun6 = @(x) 2 * (x(:, 11) + 2).^2 .* x(:, 12) .* (x(:, 12) > quantile(x(:, 12), 0.5)) ...
    .* x(:, 106)  .* (x(:, 106) < quantile(x(:, 106), 0.5)) + randn([size(x, 1), 1]);


for K = 4 : 10
    QP = (1 : (K - 1)) / K;
    p = 5000;
    tic
    for i = 1 : 100
        sigma = 0.5.^(abs(reshape(1:p,1,p)-reshape(1:p,1,p)'));
        X = mvnrnd(zeros(1, p), sigma, 500);
        Y1 = fun1(X);
        Y2 = fun2(X);
        Y3 = fun3(X);
        Y4 = fun4(X);
        Y5 = fun5(X);
        Y6 = fun6(X);

        Z1 = sum(Y1 >= quantile(Y1, QP), 2) + 1;
        Z2 = sum(Y2 >= quantile(Y2, QP), 2) + 1;
        Z3 = sum(Y3 >= quantile(Y3, QP), 2) + 1;
        Z4 = sum(Y4 >= quantile(Y4, QP), 2) + 1;
        Z5 = sum(Y5 >= quantile(Y5, QP), 2) + 1;
        Z6 = sum(Y6 >= quantile(Y6, QP), 2) + 1;

        [Fv1, Fvs1, upsilonk_hat1, upsilon_hat1] = QA_SVS(X, Z1);
        [Fv2, Fvs2, upsilonk_hat2, upsilon_hat2] = QA_SVS(X, Z2);
        [Fv3, Fvs3, upsilonk_hat3, upsilon_hat3] = QA_SVS(X, Z3);
        [Fv4, Fvs4, upsilonk_hat4, upsilon_hat4] = QA_SVS(X, Z4);
        [Fv5, Fvs5, upsilonk_hat5, upsilon_hat5] = QA_SVS(X, Z5);
        [Fv6, Fvs6, upsilonk_hat6, upsilon_hat6] = QA_SVS(X, Z6);

        A_num1(K - 3, i) = find(ismember(Fvs1(K + 1, :), [1 2 101]), 1, 'last');
        A_num2(K - 3, i) = find(ismember(Fvs2(K + 1, :), [3 4 102]), 1, 'last');
        A_num3(K - 3, i) = find(ismember(Fvs3(K + 1, :), [5 6 103]), 1, 'last');
        A_num4(K - 3, i) = find(ismember(Fvs4(K + 1, :), [7 8 104]), 1, 'last');
        A_num5(K - 3, i) = find(ismember(Fvs5(K + 1, :), [9 10 105]), 1, 'last');
        A_num6(K - 3, i) = find(ismember(Fvs6(K + 1, :), [11 12 106]), 1, 'last');

    end
    toc

end