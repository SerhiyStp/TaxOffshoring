clear; clc; %close all;

errs = load('ErrsOnGrid.txt');
GridA = load('GridX.txt');
GridH = load('GridY.txt');
na = size(GridA,1);
nh = size(GridH,1);
Errs = reshape(errs,[nh,na,2]);

TotErrSq = squeeze(Errs(:,:,1)).^2 + squeeze(Errs(:,:,2)).^2;
[M,I] = min(TotErrSq,[],"all");
[row,col] = ind2sub(size(TotErrSq),I);


figure('Position',[60,60,900,600])
subplot(1,2,1)
plot(GridH, TotErrSq(:,col))
xlabel('Hrs')
ylabel('focH^2 + bc^2')
title('Changing hours')

subplot(1,2,2)
plot(GridA, TotErrSq(row,:))
xlabel('A')
ylabel('focH^2 + bc^2')
title('Changing initial assets')

focH_a = squeeze(Errs(row,:,1));
focH_h = squeeze(Errs(:,col,1));
bc_a = squeeze(Errs(row,:,2));
bc_h = squeeze(Errs(:,col,2));

figure('Position',[60,60,900,600])
subplot(1,2,1)
plot(GridH, focH_h)
hold on
plot(GridH, bc_h, 'r--')
b = legend('foc(hrs)','bc')
set(b,'box','off','location','best')
xlabel('Hrs')
ylabel('Err')
title('Changing hours')

subplot(1,2,2)
plot(GridA, focH_a)
hold on
plot(GridA, bc_a, 'r--')
xlabel('A')
ylabel('Err')
title('Changing initial assets')
