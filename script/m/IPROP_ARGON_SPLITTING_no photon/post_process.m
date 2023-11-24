addpath (canonicalize_file_name ("../"));
addpath (canonicalize_file_name ("../BolsigPlus"));

load results2.gz
load database.txt
x=linspace(0, 2e-3, 200);
%GRARFICO el
close all
figure()
semilogy(x, u(1:200, 150))
hold on
semilogy(x, u(1:200, 200))
semilogy(x, u(1:200, end))
legend('1','2','3')
title('Electron post-modifica jacobiano')
%GRARFICO AR2+
figure()
semilogy(x, u(801:1000, 150))
hold on
semilogy(x, u(801:1000, 200))
semilogy(x, u(801:1000, end))
legend('1','2','3')
title('Ar2+ post-modifica jacobiano')
%GRARFICO AR+
figure()
semilogy(x, u(401:600, 150))
hold on
semilogy(x, u(401:600, 200))
semilogy(x, u(401:600, end))
legend('1','2','3')
title('Ar+ post-modifica jacobiano')
%campo elettrico
E1=diff(u(1201:end, 150))./diff(x);
E2=diff(u(1201:end, 200))./diff(x);
figure
plot(x(2:end), E1)
hold on
plot(x(2:end), E2)
legend('1','2')
title('E post modifica jacobiano')
##%GRARFICO el

figure()
plot(x, u(1:200, 150))
hold on
plot(x, u(1:200, 200))
plot(x, u(1:200, end))
legend('1','2','3')
title('Electron post-modifica jacobiano')
##%GRARFICO AR2+
##figure()
##plot(x, u(801:1000, 150))
##hold on
##plot(x, u(801:1000, 200))
##plot(x, u(801:1000, end))
##legend('1','2','3')
##title('Ar2+ post-modifica jacobiano')
##%GRARFICO AR+
##figure()
##plot(x, u(401:600, 150))
##hold on
##plot(x, u(401:600, 200))
##plot(x, u(401:600, end))
##legend('1','2','3')
##title('Ar+ post-modifica jacobiano')






load results1.gz
%GRARFICO el
figure()
semilogy(x, u(1:200, 150))
hold on
semilogy(x, u(1:200, 200))
legend('1','2','3')
title('Electron pre-modifica jacobiano')
%GRARFICO AR2+
figure()
semilogy(x, u(801:1000, 150))
hold on
semilogy(x, u(801:1000, 200))

legend('1','2')
title('Ar2+ pre-modifica jacobiano')
%GRARFICO AR+
figure()
semilogy(x, u(401:600, 150))
hold on
semilogy(x, u(401:600, 200))
legend('1','2')
title('Ar+ pre-modifica jacobiano')
%campo elettrico
E1=diff(u(1201:end, 150))./diff(x);
E2=diff(u(1201:end, 200))./diff(x);
figure
ploy(x(2:end), E1)
hold on
plot(x(2:end), E2)
legend('1','2')
title('E pre modifica jacobiano')
##%GRARFICO el
##figure()
##plot(x, u(1:200, 150))
##hold on
##plot(x, u(1:200, 200))
##legend('1','2','3')
##title('Electron pre-modifica jacobiano')
##%GRARFICO AR2+
##figure()
##plot(x, u(801:1000, 150))
##hold on
##plot(x, u(801:1000, 200))
##
##legend('1','2')
##title('Ar2+ pre-modifica jacobiano')
##%GRARFICO AR+
##figure()
##plot(x, u(401:600, 150))
##hold on
##plot(x, u(401:600, 200))
##legend('1','2')
##title('Ar+ pre-modifica jacobiano')
