clear;clc;
% close all;

data = load('test_plot_data.txt');
jc = data(1,1);

retired = false;
if retired
    jr = 1;
else
    jr = 100;
end

mtmp = load('mtmp.txt');
gridm = load('gridm.txt');
afun = load('afun.txt');
dA = afun(2:end)-afun(1:end-1); %./(gridm(2:end)-gridm(1:end-1));
cfun = load('cfun.txt');
Vfun = load('Vfun.txt');
dVfun =load('dVfun.txt');
ap_test = load('ap_test.txt');
offsh = load('offshoring.txt');
na = size(mtmp,1);
% ttt = amcutoff(1,2)*ones(na,1);

% figure('Position',[20,20,800,600])
% subplot(1,2,1)
% plot(mtmp(:,1),mtmp(:,2))
% xlabel('aprime')
% ylabel('m')
% subplot(1,2,2)
% plot(mtmp(:,1),mtmp(:,5))
% xlabel('aprime')
% ylabel('Offshoring')

% afun_test = zeros(501,1);
% afun_test(1:167,1) = afun(1:167,1);
% for i = 168 : 501
%     afun_test(i,1) = interp1(gridm(1:167,1),afun(1:167,1),gridm(i,1),'linear','extrap');
% end

if (jc < jr)
    hfun = load('hfun.txt');
end

alo = gridm(1,1); %30;
ahi = gridm(end,1); %60;

if (jc < jr)

    figure('Position',[20,20,1200,800])
    subplot(2,3,1)
    plot(gridm, hfun)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(hfun), max(hfun)], 'k--') 
    plot(gridm, offsh(:,1)*max(hfun), 'k--')
    plot(gridm, offsh(:,2)*max(hfun), 'g--')
    xlabel('GridA')
    xlim([alo, ahi])
    title('Hours')

    subplot(2,3,2)
    plot(gridm, dVfun)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(dVfun), max(dVfun)], 'k--')
    plot(gridm, offsh(:,1)*max(dVfun), 'k--')
    plot(gridm, offsh(:,2)*max(dVfun), 'g--')
    xlabel('GridA')
    xlim([alo, ahi])
    title('dV')

    subplot(2,3,3)
    plot(gridm, afun)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(afun), max(afun)], 'k--') 
    %plot(gridm, afun_test, 'k:')
    %plot(gridm, ttt)
    plot(gridm, offsh(:,1)*max(afun), 'k--')
    plot(gridm, offsh(:,2)*max(afun), 'g--')
    xlabel('GridA')
    xlim([alo, ahi])
    plot(gridm, gridm, 'k:')
    plot(gridm, ap_test, 'r-.')
    % xlim([20, 40])
    title('A(t+1)')

    subplot(2,3,4)
    %plot(gridm, cfun, '-bd', 'MarkerSize', 0.5)
    plot(gridm, cfun)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(cfun), max(cfun)], 'k--') 
    plot(gridm, offsh(:,1)*max(cfun), 'k--')
    plot(gridm, offsh(:,2)*max(cfun), 'g--')
    xlabel('GridA')
    xlim([alo, ahi])
    title('Consumption')

    subplot(2,3,5)
    plot(gridm, Vfun)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(Vfun), max(Vfun)], 'k--') 
    xlabel('GridM')
    xlim([alo, ahi])
    title('Value fn')

    subplot(2,3,6)
    plot(gridm(1:end-1), dA)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(Vfun), max(Vfun)], 'k--')
    plot(gridm, offsh(:,1)*max(dA), 'k--')
    plot(gridm, offsh(:,2)*max(dA), 'g--')
    xlabel('GridA')
    xlim([alo, ahi])
    title('dA(t+1)')    
    
    common_title = sprintf('Age = %d', jc);
    sgtitle(common_title)   

    figure
    plot(gridm, afun-ap_test)

else
    figure('Position',[20,20,1200,800])
    subplot(2,2,1)
    plot(gridm, dVfun)
    %hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(dVfun), max(dVfun)], 'k--') 
    xlabel('GridM')
    xlim([gridm(1,1), gridm(end,1)])
    title('dV')

    subplot(2,2,2)
    %plot(gridm, afun, '-bd', 'MarkerSize', 0.5)
    plot(gridm, afun)
    hold on
    plot(gridm, gridm, 'k:')
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(afun), max(afun)], 'k--')
    %plot(gridm, ttt, 'k--')
    xlabel('GridM')
    xlim([gridm(1,1), gridm(end,1)])
    % xlim([20, 40])
    title('A(t+1)')

    subplot(2,2,3)
    %plot(gridm, cfun, '-bd', 'MarkerSize', 0.5)
    plot(gridm, cfun)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(cfun), max(cfun)], 'k--') 
    %plot(mtmp(:,1),mtmp(:,5)*max(cfun),'k-.')
    xlabel('GridM')
    xlim([gridm(1,1), gridm(end,1)])
    title('Consumption')

    subplot(2,2,4)
    plot(gridm, Vfun)
    %hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(Vfun), max(Vfun)], 'k--') 
    xlabel('GridM')
    xlim([gridm(1,1), gridm(end,1)])
    title('Value fn')    
    
    common_title = sprintf('Age = %d', jc);
    sgtitle(common_title)
end
