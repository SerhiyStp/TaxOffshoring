clear; clc; close all;

data = load('LifeCycle.txt');

figure('Position',[60,60,1200,600])
subplot(1,3,1)
plot(19+data(:,1), data(:,2))
xlabel('age')
title('Savings')
subplot(1,3,2)
plot(19+data(:,1), data(:,3))
xlabel('age')
title('Hours')
subplot(1,3,3)
plot(19+data(:,1), data(:,4))
xlabel('age')
title('Labor efficiency units')