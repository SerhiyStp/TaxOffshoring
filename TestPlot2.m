clear;clc;
close all;
amcutoff = load('am_cutoff.txt');

mdat = load('mtmp.txt');
grida = mdat(:,1);
m = mdat(:,2);
c = mdat(:,3);
% h = mdat(:,4);
% V = mdat(:,5);
V = mdat(:,4);

figure
plot(mdat(:,1),mdat(:,2))

figure('Position', [20, 20, 1000,600])
subplot(1,2,1)
plot(m,c)
hold on
plot([amcutoff(1,2), amcutoff(1,2)], [min(c), max(c)], 'k--')
% xlim([30, 31])
title('Consumption')

subplot(1,2,2)
plot(m,V)
hold on
plot([amcutoff(1,2), amcutoff(1,2)], [min(V), max(V)], 'k--')
% xlim([30, 31])
title('Value fn')