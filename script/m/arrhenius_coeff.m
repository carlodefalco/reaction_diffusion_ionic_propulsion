
k_b= 1.380649e-23;%[J K^-1]
Mean_Energy=5; %[eV]
T_e=2/3*Mean_Energy;%[eV]
T_eKelvin=T_e*1.60218e-19/k_b;


Activation_Excitation=11.5; %[eV]
Excitation=4.9e-15*T_e^0.5*exp(-Activation_Excitation/T_e); %[m^3/s]

Activation_Ionization=15.8;%[eV]
Ionization=1.27e-14*T_e^0.5*exp(-Activation_Excitation/T_e); %[m^3/s]

Activation_IonizationArex=4.2;%[eV]
Ionization_Arex=1.37e-13*T_e^0.5*exp(-Activation_Excitation/T_e); %[m^3/s]

Excitation_Ar2plus=7e-13*(300/T_eKelvin)^0.5;




