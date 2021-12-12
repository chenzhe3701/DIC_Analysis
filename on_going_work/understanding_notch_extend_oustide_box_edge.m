% Notches can extend outside of the box edges.
% http://mohinish.blogspot.com/2010/10/matlab-boxplot-notch-error-resolved.html
a = [
    0
    0.2446
    0.2531
    0.2767
    0.2810
    0.2844
    0.3639
    0.4319
    0.4557
    0.4978
    0.5001
    0.5044
    0.5168
    0.5417
    0.5708
    0.5974
    0.7494
    0.9022
    ];

figure; boxplot(a,'notch','on');
figure; boxplot(a)