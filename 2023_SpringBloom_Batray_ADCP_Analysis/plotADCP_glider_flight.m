load('ADCP+glider_data.mat')

figure()
subplot(411)
plot(adat.time,adat.heading-adat.glider_heading,'.r'); hold on
ylabel('ADCP-glider heading')
datetick('x','dd.mmm','keeplimits')

subplot(412)
plot(adat.time,adat.pitch-adat.glider_pitch,'.m'); hold on
ylabel('ADCP-glider pitch')
datetick('x','dd.mmm','keeplimits')

subplot(413)
plot(adat.time,adat.roll-adat.glider_roll,'.b'); hold on
ylabel('ADCP-glider roll')
datetick('x','dd.mmm','keeplimits')

subplot(414)
plot(adat.time,adat.pres-adat.glider_pres,'.k'); hold on
ylabel('ADCP-glider pressure')
datetick('x','dd.mmm','keeplimits')

figure()
plot(adat.time,adat.ADCP_sos-adat.glider_CTD_sos,'.c'); hold on
ylabel('ADCP-glider Speed of Sound')
datetick('x','dd.mmm','keeplimits')