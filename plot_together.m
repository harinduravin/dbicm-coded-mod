

fig1 = openfig('BICM_LDPC_DVB_S2_rate_2_3_v4.fig','invisible');
fig2 = openfig('BICM_LDPC_Wimax_rate_2_3_v3.fig','invisible');

figure
h(1) = subplot(1,1,1);
grid
legend('a')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
set(gca,'YScale','log')
xline(4.3517,'Label','Shannon limit = 4.35 dB (rate 2/3 | 32-ary modulation)')
xlim([0 7.8])
copyobj(allchild(get(fig1,'CurrentAxes')),h(1));
copyobj(allchild(get(fig2,'CurrentAxes')),h(1));