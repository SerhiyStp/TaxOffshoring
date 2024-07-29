clear;clc;
close all;

data = load('test_plot_data.txt');
jc = data(1,1);

retired = false;

mtmp = load('mtmp.txt');
gridm = load('gridm.txt');
afun = load('afun.txt');
dA = afun(2:end)-afun(1:end-1); %./(gridm(2:end)-gridm(1:end-1));
cfun = load('cfun.txt');
Vfun = load('Vfun.txt');
dVfun =load('dVfun.txt');
ap_test = load('ap_test.txt');
offsh = load('offshoring.txt');
capinc = load('capinc.txt');
labinc = load('labinc.txt');
ctax = load('constax.txt');
na = size(mtmp,1);
gross_inc = load('pretax_inc.txt');
% ttt = amcutoff(1,2)*ones(na,1);

if (retired == false)
    hfun = load('hfun.txt');
end

alo = gridm(1,1); %30;
ahi = gridm(end,1); %60;

if (retired == false)

    figure
    plot(gridm, gross_inc)
    xlabel('GridA')
    title('Pre-tax income')    

    figure('Position',[20,20,1400,800])
    subplot(2,4,1)
    plot(gridm, hfun)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(hfun), max(hfun)], 'k--') 
    plot(gridm, offsh(:,1)*max(hfun), 'k--')
    plot(gridm, offsh(:,2)*max(hfun), 'g--')
    xlabel('GridA')
    xlim([alo, ahi])
    title('Hours')

    subplot(2,4,2)
    plot(gridm, dVfun)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(dVfun), max(dVfun)], 'k--')
    plot(gridm, offsh(:,1)*max(dVfun), 'k--')
    plot(gridm, offsh(:,2)*max(dVfun), 'g--')
    xlabel('GridA')
    xlim([alo, ahi])
    title('dV')

    subplot(2,4,3)
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

    subplot(2,4,4)
    %plot(gridm, cfun, '-bd', 'MarkerSize', 0.5)
    plot(gridm, cfun)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(cfun), max(cfun)], 'k--') 
    plot(gridm, offsh(:,1)*max(cfun), 'k--')
    plot(gridm, offsh(:,2)*max(cfun), 'g--')
    xlabel('GridA')
    xlim([alo, ahi])
    title('Consumption')

    subplot(2,4,5)
    plot(gridm, Vfun)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(Vfun), max(Vfun)], 'k--') 
    xlabel('GridA')
    xlim([alo, ahi])
    title('Value fn')

    subplot(2,4,6)
    plot(gridm(1:end-1), dA)
    hold on
    %plot([amcutoff(1,2), amcutoff(1,2)], [min(Vfun), max(Vfun)], 'k--')
    plot(gridm, offsh(:,1)*max(dA), 'k--')
    plot(gridm, offsh(:,2)*max(dA), 'g--')
    xlabel('GridA')
    xlim([alo, ahi])
    title('dA(t+1)')  

    subplot(2,4,7)
    plot(gridm, capinc)
    xlabel('GridA')
    xlim([alo, ahi])
    title('Capital income')   

    subplot(2,4,8)
    plot(gridm, labinc)
    xlabel('GridA')
    xlim([alo, ahi])
    title('Labor income')     
    
    common_title = sprintf('Age = %d', jc);
    sgtitle(common_title)   

    figure
    plot(gridm, afun-ap_test)

    figure
    plot(gridm, ctax)
    xlabel('GridA')
    xlim([alo, ahi])
    title('Consumption tax')    

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
