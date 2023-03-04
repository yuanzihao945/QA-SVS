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


p = 5000;
tic
for i = 1 : 100
    
    X = randn(500, p);
    Y1 = fun1(X);
    Y2 = fun2(X);
    Y3 = fun3(X);
    Y4 = fun4(X);
    Y5 = fun5(X);
    Y6 = fun6(X);

    [~, Fvs1] = sort(abs(Y1' * X), 'descend');
    [~, Fvs2] = sort(abs(Y2' * X), 'descend');
    [~, Fvs3] = sort(abs(Y3' * X), 'descend');
    [~, Fvs4] = sort(abs(Y4' * X), 'descend');
    [~, Fvs5] = sort(abs(Y5' * X), 'descend');
    [~, Fvs6] = sort(abs(Y6' * X), 'descend');

    A_num1(i) = find(ismember(Fvs1, [1 2 101]), 1, 'last');
    A_num2(i) = find(ismember(Fvs2, [3 4 102]), 1, 'last');
    A_num3(i) = find(ismember(Fvs3, [5 6 103]), 1, 'last');
    A_num4(i) = find(ismember(Fvs4, [7 8 104]), 1, 'last');
    A_num5(i) = find(ismember(Fvs5, [9 10 105]), 1, 'last');
    A_num6(i) = find(ismember(Fvs6, [11 12 106]), 1, 'last');

end
toc
