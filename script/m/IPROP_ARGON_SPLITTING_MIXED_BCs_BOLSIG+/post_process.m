load results.gz
figure
plot(x, u(1:200, 1), 'LineWidth', 1.5)
hold on
plot(x, u(1:200, 15), 'LineWidth', 1.5)
plot(x, u(1:200, 30), 'LineWidth', 1.5)
plot(x, u(1:200, 60), 'LineWidth', 1.5)
plot(x, u(1:200, end), 'LineWidth', 1.5)
legend('Location', 'north')
title('electron density')
print ("-dpng", sprintf ("electron_density.png"))

figure
plot(x, u(401:600, 1), 'LineWidth', 1.5)
hold on
plot(x, u(401:600, 15), 'LineWidth', 1.5)
plot(x, u(401:600, 30), 'LineWidth', 1.5)
plot(x, u(401:600, 60), 'LineWidth', 1.5)
plot(x, u(401:600, end), 'LineWidth', 1.5)
legend('Location', 'north')
title('Ar+')
print ("-dpng", sprintf ("Arplus.png"))
figure
plot(x, u(801:1000, 1), 'LineWidth', 1.5)
hold on
plot(x, u(801:1000, 15), 'LineWidth', 1.5)
plot(x, u(801:1000, 30), 'LineWidth', 1.5)
plot(x, u(801:1000, 60), 'LineWidth', 1.5)
plot(x, u(801:1000, end), 'LineWidth', 1.5)
legend('Location', 'north')
title('Ar2+')
print ("-dpng", sprintf ("Ar2plus.png"))

E1=-gradient(u(1201:1400,1),x); %[V/m]
E2=-gradient(u(1201:1400,15),x);%[V/m]
E3=-gradient(u(1201:1400,30),x);%[V/m]
E4=-gradient(u(1201:1400,60),x);%[V/m]
E5=-gradient(u(1201:1400,end),x);%[V/m]

E1Td=E1./u(201:400, 1).*1e21; %[Td]
E2Td=E2./u(201:400, 15).*1e21;%[Td]
E3Td=E3./u(201:400, 30).*1e21;%[Td]
E4Td=E4./u(201:400, 60).*1e21;%[Td]
E5Td=E5./u(201:400, end).*1e21;%[Td]


figure
plot(x, E1, 'LineWidth', 1.5)
hold on
plot(x, E2, 'LineWidth', 1.5)
plot(x, E3, 'LineWidth', 1.5)
plot(x, E4, 'LineWidth', 1.5)
plot(x, E5, 'LineWidth', 1.5)
legend('Location', 'north')
title('E=-gradient(phi,x)')
print ("-dpng", sprintf ("E_SI.png"))

figure
plot(x, E1Td, 'LineWidth', 1.5)
hold on
plot(x, E2Td, 'LineWidth', 1.5)
plot(x, E3Td, 'LineWidth', 1.5)
plot(x, E4Td, 'LineWidth', 1.5)
plot(x, E5Td, 'LineWidth', 1.5)
legend('Location', 'north')
title('E(Td)')
print ("-dpng", sprintf ("E_Td.png"))
figure
plot(x, Te)
title('Electron temperature at the last iteration')







print ("-dpng", sprintf ("Te.png"))
