load lungCT.mat

for K = 2 : 10

    QP = (1 : (K - 1)) / K;

    Z95 = sum(Y95 >= quantile(Y95, QP), 2) + 1;
    Z05 = sum(Y05 >= quantile(Y05, QP), 2) + 1;

    [Fv95AFD, Fvs95AFD] = QA_SVS(lungCT, Z95, 'AFD');
    [Fv05AFD, Fvs05AFD] = QA_SVS(lungCT, Z05, 'AFD');

    Fv95QAAFD{K-1} = Fv95AFD{K+1};
    Fv05QAAFD{K-1} = Fv05AFD{K+1};

    [Fv95FDR, Fvs95FDR] = QA_SVS(lungCT, Z95, 0.05);
    [Fv05FDR, Fvs05FDR] = QA_SVS(lungCT, Z05, 0.05);

    Fv95QAFDR{K-1} = Fv95FDR{K+1};
    Fv05QAFDR{K-1} = Fv05FDR{K+1};

    [Fv95QCFDR{K-1}, Fvs95QCFDR] = QCS_FDR(lungCT, Z95, 0.05);
    [Fv05QCFDR{K-1}, Fvs05QCFDR] = QCS_FDR(lungCT, Z05, 0.05);

    disp(K-1)

end

lungCTmean = mean(lungCT);
for i = 1 : 5
    lungtem = zeros(512, 512);
    lungtem(Fv95QAFDR{i}) = 1;
    figure(i)
    imagesc(lungtem);
    set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    colormap('gray');
    saveas(figure(i), ['95PD_', num2str(i+1), '.pdf'], 'pdf')
end

for i = 1 : 5
    lungtem = zeros(512, 512);
    lungtem(Fv05QAFDR{i}) = 1;
    figure(i)
    imagesc(lungtem);
    set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    colormap('gray');
    saveas(figure(i), ['05PD_', num2str(i+1), '.pdf'], 'pdf')
end