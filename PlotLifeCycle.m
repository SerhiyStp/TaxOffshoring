clear; clc; close all;

data = load('LifeCycle.txt');

figure('Position',[60,60,900,600])
subplot(1,2,1)
plot(19+data(:,1), data(:,2))
xlabel('age')
title('Savings')
subplot(1,2,2)
plot(19+data(:,1), data(:,3))
xlabel('age')
title('Hours')