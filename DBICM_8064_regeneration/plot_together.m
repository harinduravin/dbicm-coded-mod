

fig1 = openfig('BICM.fig','invisible');
fig2 = openfig('Soft_DBICM.fig','invisible');
fig3 = openfig('Lower_bound.fig','invisible');
% fig4 = openfig('Hard_DBICM.fig','invisible');

figure
h(1) = subplot(1,1,1);
grid
legend('a')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
set(gca,'YScale','log')
% xline(2.226,'Label','Shannon limit = 2.23 dB (rate 3/4 | 8-APSK modulation)')
xlim([3.8 5.3])
copyobj(allchild(get(fig1,'CurrentAxes')),h(1));
copyobj(allchild(get(fig2,'CurrentAxes')),h(1));
copyobj(allchild(get(fig3,'CurrentAxes')),h(1));
% copyobj(allchild(get(fig4,'CurrentAxes')),h(1));