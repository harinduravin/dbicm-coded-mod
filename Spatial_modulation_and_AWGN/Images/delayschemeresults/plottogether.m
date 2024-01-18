fig1 = openfig('delay1.fig','invisible');
fig2 = openfig('delay2.fig','invisible');
fig3 = openfig('delay3.fig','invisible');
% fig4 = openfig('Hard_DBICM.fig','invisible');

figure
h(1) = subplot(1,1,1);
grid
legend('[0,0,1,0,1,1,0,1]','[0,0,1,1,0,1,1,0]','[0,0,1,1,0,0,0,0]')
xlabel('Eb/No (dB)')
ylabel('Total Channel Capacity')
% set(gca,'YScale','log')
% xline(2.226,'Label','Shannon limit = 2.23 dB (rate 3/4 | 8-APSK modulation)')
% xlim([3.8 5.3])
copyobj(allchild(get(fig1,'CurrentAxes')),h(1));
copyobj(allchild(get(fig2,'CurrentAxes')),h(1));
copyobj(allchild(get(fig3,'CurrentAxes')),h(1));
% copyobj(allchild(get(fig4,'CurrentAxes')),h(1));