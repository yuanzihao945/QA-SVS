clearvars
close all
clc

fun1 = @(x) sum(x(:, 1:10), 2) + randn([size(x, 1), 1]);
fun2 = @(x) sum(x(:, 1:50), 2) + randn([size(x, 1), 1]);
fun3 = @(x) exp(sum(x(:, 1:10), 2)) + randn([size(x, 1), 1]);
fun4 = @(x) exp(sum(x(:, 1:50), 2)) + randn([size(x, 1), 1]);
fun5 = @(x) sum((-1).^(0:9) .* x(:, 1:10), 2) + 0.1 * randn([size(x, 1), 1]);
fun6 = @(x) sum((-1).^(0:49) .* x(:, 1:50), 2) + 0.1 * randn([size(x, 1), 1]);
fun7 = @(x) sum(x(:, 1:10), 2) ./ (0.5 + (1.5 + sum((-1).^(0:2) .* x(:, 2:4), 2)).^2) + 0.1*randn([size(x, 1), 1]);
fun8 = @(x) sum(x(:, 1:50), 2) ./ (0.5 + (1.5 + sum((-1).^(0:19) .* x(:, 21:40), 2)).^2) + 0.1*randn([size(x, 1), 1]);


n = 1000;
p = 1000;
sigma1 = 0.5.^(abs(reshape(1:p,1,p) - reshape(1:p,1,p)'));
sigma2 = eye(p) + 0.2 * diag(ones(1,p-1),1) + 0.2 * diag(ones(1,p-1),-1);


for K = 2:6

    for i = 1 : 100

        X1 = mvnrnd(zeros(1, p), sigma1, n);
        X2 = mvnrnd(zeros(1, p), sigma1, n);
        X3 = mvnrnd(zeros(1, p), sigma1, n);
        X4 = mvnrnd(zeros(1, p), sigma1, n);
        X5 = mvnrnd(zeros(1, p), sigma2, n);
        X6 = mvnrnd(zeros(1, p), sigma2, n);
        X7 = mvnrnd(zeros(1, p), sigma2, n);
        X8 = mvnrnd(zeros(1, p), sigma2, n);

        Y1 = fun1(X1);
        Y2 = fun2(X2);
        Y3 = fun3(X3);
        Y4 = fun4(X4);
        Y5 = fun5(X5);
        Y6 = fun6(X6);
        Y7 = fun7(X7);
        Y8 = fun8(X8);

        QP = (1 : (K - 1)) / K;

        Z1 = sum(Y1 >= quantile(Y1, QP), 2) + 1;
        Z2 = sum(Y2 >= quantile(Y2, QP), 2) + 1;
        Z3 = sum(Y3 >= quantile(Y3, QP), 2) + 1;
        Z4 = sum(Y4 >= quantile(Y4, QP), 2) + 1;
        Z5 = sum(Y5 >= quantile(Y5, QP), 2) + 1;
        Z6 = sum(Y6 >= quantile(Y6, QP), 2) + 1;
        Z7 = sum(Y7 >= quantile(Y7, QP), 2) + 1;
        Z8 = sum(Y8 >= quantile(Y8, QP), 2) + 1;

        [Fv1AFD, Fvs1AFD] = QA_SVS(X1, Z1, 'AFD');
        [Fv2AFD, Fvs2AFD] = QA_SVS(X2, Z2, 'AFD');
        [Fv3AFD, Fvs3AFD] = QA_SVS(X3, Z3, 'AFD');
        [Fv4AFD, Fvs4AFD] = QA_SVS(X4, Z4, 'AFD');
        [Fv5AFD, Fvs5AFD] = QA_SVS(X5, Z5, 'AFD');
        [Fv6AFD, Fvs6AFD] = QA_SVS(X6, Z6, 'AFD');
        [Fv7AFD, Fvs7AFD] = QA_SVS(X7, Z7, 'AFD');
        [Fv8AFD, Fvs8AFD] = QA_SVS(X8, Z8, 'AFD');

        Fv1QAAFD = Fv1AFD{K+1};
        Fv2QAAFD = Fv2AFD{K+1};
        Fv3QAAFD = Fv3AFD{K+1};
        Fv4QAAFD = Fv4AFD{K+1};
        Fv5QAAFD = Fv5AFD{K+1};
        Fv6QAAFD = Fv6AFD{K+1};
        Fv7QAAFD = Fv7AFD{K+1};
        Fv8QAAFD = Fv8AFD{K+1};

        [Fv1FDR, Fvs1FDR] = QA_SVS(X1, Z1, 0.05);
        [Fv2FDR, Fvs2FDR] = QA_SVS(X2, Z2, 0.05);
        [Fv3FDR, Fvs3FDR] = QA_SVS(X3, Z3, 0.05);
        [Fv4FDR, Fvs4FDR] = QA_SVS(X4, Z4, 0.05);
        [Fv5FDR, Fvs5FDR] = QA_SVS(X5, Z5, 0.05);
        [Fv6FDR, Fvs6FDR] = QA_SVS(X6, Z6, 0.05);
        [Fv7FDR, Fvs7FDR] = QA_SVS(X7, Z7, 0.05);
        [Fv8FDR, Fvs8FDR] = QA_SVS(X8, Z8, 0.05);

        Fv1QAFDR = Fv1FDR{K+1};
        Fv2QAFDR = Fv2FDR{K+1};
        Fv3QAFDR = Fv3FDR{K+1};
        Fv4QAFDR = Fv4FDR{K+1};
        Fv5QAFDR = Fv5FDR{K+1};
        Fv6QAFDR = Fv6FDR{K+1};
        Fv7QAFDR = Fv7FDR{K+1};
        Fv8QAFDR = Fv8FDR{K+1};

        [Fv1QCFDR, Fvs1QCFDR] = QCS_FDR(X1, Z1, 0.05);
        [Fv2QCFDR, Fvs2QCFDR] = QCS_FDR(X2, Z2, 0.05);
        [Fv3QCFDR, Fvs3QCFDR] = QCS_FDR(X3, Z3, 0.05);
        [Fv4QCFDR, Fvs4QCFDR] = QCS_FDR(X4, Z4, 0.05);
        [Fv5QCFDR, Fvs5QCFDR] = QCS_FDR(X5, Z5, 0.05);
        [Fv6QCFDR, Fvs6QCFDR] = QCS_FDR(X6, Z6, 0.05);
        [Fv7QCFDR, Fvs7QCFDR] = QCS_FDR(X7, Z7, 0.05);
        [Fv8QCFDR, Fvs8QCFDR] = QCS_FDR(X8, Z8, 0.05);

        Snumber_QAAFD1(K-1, i) = length(Fv1QAAFD);
        Snumber_QAAFD2(K-1, i) = length(Fv2QAAFD);
        Snumber_QAAFD3(K-1, i) = length(Fv3QAAFD);
        Snumber_QAAFD4(K-1, i) = length(Fv4QAAFD);
        Snumber_QAAFD5(K-1, i) = length(Fv5QAAFD);
        Snumber_QAAFD6(K-1, i) = length(Fv6QAAFD);
        Snumber_QAAFD7(K-1, i) = length(Fv7QAAFD);
        Snumber_QAAFD8(K-1, i) = length(Fv8QAAFD);

        Snumber_QAFDR1(K-1, i) = length(Fv1QAFDR);
        Snumber_QAFDR2(K-1, i) = length(Fv2QAFDR);
        Snumber_QAFDR3(K-1, i) = length(Fv3QAFDR);
        Snumber_QAFDR4(K-1, i) = length(Fv4QAFDR);
        Snumber_QAFDR5(K-1, i) = length(Fv5QAFDR);
        Snumber_QAFDR6(K-1, i) = length(Fv6QAFDR);
        Snumber_QAFDR7(K-1, i) = length(Fv7QAFDR);
        Snumber_QAFDR8(K-1, i) = length(Fv8QAFDR);

        Snumber_QCFDR1(K-1, i) = length(Fv1QCFDR);
        Snumber_QCFDR2(K-1, i) = length(Fv2QCFDR);
        Snumber_QCFDR3(K-1, i) = length(Fv3QCFDR);
        Snumber_QCFDR4(K-1, i) = length(Fv4QCFDR);
        Snumber_QCFDR5(K-1, i) = length(Fv5QCFDR);
        Snumber_QCFDR6(K-1, i) = length(Fv6QCFDR);
        Snumber_QCFDR7(K-1, i) = length(Fv7QCFDR);
        Snumber_QCFDR8(K-1, i) = length(Fv8QCFDR);

        FDR_QAAFD1(K-1, i) = sum(Fv1QAAFD > 10) / length(Fv1QAAFD);
        FDR_QAAFD2(K-1, i) = sum(Fv2QAAFD > 50) / length(Fv2QAAFD);
        FDR_QAAFD3(K-1, i) = sum(Fv3QAAFD > 10) / length(Fv3QAAFD);
        FDR_QAAFD4(K-1, i) = sum(Fv4QAAFD > 50) / length(Fv4QAAFD);
        FDR_QAAFD5(K-1, i) = sum(Fv5QAAFD > 10) / length(Fv5QAAFD);
        FDR_QAAFD6(K-1, i) = sum(Fv6QAAFD > 50) / length(Fv6QAAFD);
        FDR_QAAFD7(K-1, i) = sum(Fv7QAAFD > 10) / length(Fv7QAAFD);
        FDR_QAAFD8(K-1, i) = sum(Fv8QAAFD > 50) / length(Fv8QAAFD);

        FDR_QAFDR1(K-1, i) = sum(Fv1QAFDR > 10) / length(Fv1QAFDR);
        FDR_QAFDR2(K-1, i) = sum(Fv2QAFDR > 50) / length(Fv2QAFDR);
        FDR_QAFDR3(K-1, i) = sum(Fv3QAFDR > 10) / length(Fv3QAFDR);
        FDR_QAFDR4(K-1, i) = sum(Fv4QAFDR > 50) / length(Fv4QAFDR);
        FDR_QAFDR5(K-1, i) = sum(Fv5QAFDR > 10) / length(Fv5QAFDR);
        FDR_QAFDR6(K-1, i) = sum(Fv6QAFDR > 50) / length(Fv6QAFDR);
        FDR_QAFDR7(K-1, i) = sum(Fv7QAFDR > 10) / length(Fv7QAFDR);
        FDR_QAFDR8(K-1, i) = sum(Fv8QAFDR > 50) / length(Fv8QAFDR);

        FDR_QCFDR1(K-1, i) = sum(Fv1QCFDR > 10) / length(Fv1QCFDR);
        FDR_QCFDR2(K-1, i) = sum(Fv2QCFDR > 50) / length(Fv2QCFDR);
        FDR_QCFDR3(K-1, i) = sum(Fv3QCFDR > 10) / length(Fv3QCFDR);
        FDR_QCFDR4(K-1, i) = sum(Fv4QCFDR > 50) / length(Fv4QCFDR);
        FDR_QCFDR5(K-1, i) = sum(Fv5QCFDR > 10) / length(Fv5QCFDR);
        FDR_QCFDR6(K-1, i) = sum(Fv6QCFDR > 50) / length(Fv6QCFDR);
        FDR_QCFDR7(K-1, i) = sum(Fv7QCFDR > 10) / length(Fv7QCFDR);
        FDR_QCFDR8(K-1, i) = sum(Fv8QCFDR > 50) / length(Fv8QCFDR);

        F1score_QAAFD1(K-1, i) = 2 * sum(Fv1QAAFD <= 10) / (length(Fv1QAAFD) + 10);
        F1score_QAAFD2(K-1, i) = 2 * sum(Fv2QAAFD <= 50) / (length(Fv2QAAFD) + 50);
        F1score_QAAFD3(K-1, i) = 2 * sum(Fv3QAAFD <= 10) / (length(Fv3QAAFD) + 10);
        F1score_QAAFD4(K-1, i) = 2 * sum(Fv4QAAFD <= 50) / (length(Fv4QAAFD) + 50);
        F1score_QAAFD5(K-1, i) = 2 * sum(Fv5QAAFD <= 10) / (length(Fv5QAAFD) + 10);
        F1score_QAAFD6(K-1, i) = 2 * sum(Fv6QAAFD <= 50) / (length(Fv6QAAFD) + 50);
        F1score_QAAFD7(K-1, i) = 2 * sum(Fv7QAAFD <= 10) / (length(Fv7QAAFD) + 10);
        F1score_QAAFD8(K-1, i) = 2 * sum(Fv8QAAFD <= 50) / (length(Fv8QAAFD) + 50);

        F1score_QAFDR1(K-1, i) = 2 * sum(Fv1QAFDR <= 10) / (length(Fv1QAFDR) + 10);
        F1score_QAFDR2(K-1, i) = 2 * sum(Fv2QAFDR <= 50) / (length(Fv2QAFDR) + 50);
        F1score_QAFDR3(K-1, i) = 2 * sum(Fv3QAFDR <= 10) / (length(Fv3QAFDR) + 10);
        F1score_QAFDR4(K-1, i) = 2 * sum(Fv4QAFDR <= 50) / (length(Fv4QAFDR) + 50);
        F1score_QAFDR5(K-1, i) = 2 * sum(Fv5QAFDR <= 10) / (length(Fv5QAFDR) + 10);
        F1score_QAFDR6(K-1, i) = 2 * sum(Fv6QAFDR <= 50) / (length(Fv6QAFDR) + 50);
        F1score_QAFDR7(K-1, i) = 2 * sum(Fv7QAFDR <= 10) / (length(Fv7QAFDR) + 10);
        F1score_QAFDR8(K-1, i) = 2 * sum(Fv8QAFDR <= 50) / (length(Fv8QAFDR) + 50);

        F1score_QCFDR1(K-1, i) = 2 * sum(Fv1QCFDR <= 10) / (length(Fv1QCFDR) + 10);
        F1score_QCFDR2(K-1, i) = 2 * sum(Fv2QCFDR <= 50) / (length(Fv2QCFDR) + 50);
        F1score_QCFDR3(K-1, i) = 2 * sum(Fv3QCFDR <= 10) / (length(Fv3QCFDR) + 10);
        F1score_QCFDR4(K-1, i) = 2 * sum(Fv4QCFDR <= 50) / (length(Fv4QCFDR) + 50);
        F1score_QCFDR5(K-1, i) = 2 * sum(Fv5QCFDR <= 10) / (length(Fv5QCFDR) + 10);
        F1score_QCFDR6(K-1, i) = 2 * sum(Fv6QCFDR <= 50) / (length(Fv6QCFDR) + 50);
        F1score_QCFDR7(K-1, i) = 2 * sum(Fv7QCFDR <= 10) / (length(Fv7QCFDR) + 10);
        F1score_QCFDR8(K-1, i) = 2 * sum(Fv8QCFDR <= 50) / (length(Fv8QCFDR) + 50);

    end
end

CriteriaAll = [mean(Snumber_QAAFD1, 2)' mean(Snumber_QAFDR1, 2)' mean(Snumber_QCFDR1, 2)';...
    mean(FDR_QAAFD1, 2)' mean(FDR_QAFDR1, 2)' mean(FDR_QCFDR1, 2)';...
    mean(F1score_QAAFD1, 2)' mean(F1score_QAFDR1, 2)' mean(F1score_QCFDR1, 2)';...
    mean(Snumber_QAAFD2, 2)' mean(Snumber_QAFDR2, 2)' mean(Snumber_QCFDR2, 2)';...
    mean(FDR_QAAFD2, 2)' mean(FDR_QAFDR2, 2)' mean(FDR_QCFDR2, 2)';...
    mean(F1score_QAAFD2, 2)' mean(F1score_QAFDR2, 2)' mean(F1score_QCFDR2, 2)';...
    mean(Snumber_QAAFD3, 2)' mean(Snumber_QAFDR3, 2)' mean(Snumber_QCFDR3, 2)';...
    mean(FDR_QAAFD3, 2)' mean(FDR_QAFDR3, 2)' mean(FDR_QCFDR3, 2)';...
    mean(F1score_QAAFD3, 2)' mean(F1score_QAFDR3, 2)' mean(F1score_QCFDR3, 2)';...
    mean(Snumber_QAAFD4, 2)' mean(Snumber_QAFDR4, 2)' mean(Snumber_QCFDR4, 2)';...
    mean(FDR_QAAFD4, 2)' mean(FDR_QAFDR4, 2)' mean(FDR_QCFDR4, 2)';...
    mean(F1score_QAAFD4, 2)' mean(F1score_QAFDR4, 2)' mean(F1score_QCFDR4, 2)';...
    mean(Snumber_QAAFD5, 2)' mean(Snumber_QAFDR5, 2)' mean(Snumber_QCFDR5, 2)';...
    mean(FDR_QAAFD5, 2)' mean(FDR_QAFDR5, 2)' mean(FDR_QCFDR5, 2)';...
    mean(F1score_QAAFD5, 2)' mean(F1score_QAFDR5, 2)' mean(F1score_QCFDR5, 2)';...
    mean(Snumber_QAAFD6, 2)' mean(Snumber_QAFDR6, 2)' mean(Snumber_QCFDR6, 2)';...
    mean(FDR_QAAFD6, 2)' mean(FDR_QAFDR6, 2)' mean(FDR_QCFDR6, 2)';...
    mean(F1score_QAAFD6, 2)' mean(F1score_QAFDR6, 2)' mean(F1score_QCFDR6, 2)';...
    mean(Snumber_QAAFD7, 2)' mean(Snumber_QAFDR7, 2)' mean(Snumber_QCFDR7, 2)';...
    mean(FDR_QAAFD7, 2)' mean(FDR_QAFDR7, 2)' mean(FDR_QCFDR7, 2)';...
    mean(F1score_QAAFD7, 2)' mean(F1score_QAFDR7, 2)' mean(F1score_QCFDR7, 2)';...
    mean(Snumber_QAAFD8, 2)' mean(Snumber_QAFDR8, 2)' mean(Snumber_QCFDR8, 2)';...
    mean(FDR_QAAFD8, 2)' mean(FDR_QAFDR8, 2)' mean(FDR_QCFDR8, 2)';...
    mean(F1score_QAAFD8, 2)' mean(F1score_QAFDR8, 2)' mean(F1score_QCFDR8, 2)';...
    ];

writematrix(CriteriaAll, 'Simulation2p1000.csv')