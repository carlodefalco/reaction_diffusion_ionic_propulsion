function reactions=compute_coeff( y, x,reactions, index,N_species, database)
 N=numel(x);
 L=x(end);
 k_b     = 1.380649e-23; %Boltzmann constant [J K^-1]


 numx = numel(unique(database(:,1)));
 numy = numel(unique(database(:,3)));
 numz = numel(unique(database(:,4)));

 EF_Reduced            =reshape(database(:,1),[numy, numx, numz]); %Electric_field/N(IN) (Td)
 Plasma_density        =reshape(database(:,3),[numy, numx, numz]);% Plasma Density(IN) [1/m^3]
 Ar_Ex_Fraction        =reshape(database(:,4),[numy, numx, numz]);%Excited Argon Fraction(IN) [1]
 %Ionization_Degree     =reshape(database(:,2),[numy, numx, numz]);% Ionization_Degree(IN) [1]
 Mean_Energy           =reshape(database(:,6),[numy, numx, numz]);% Electron Mean Energy [eV]
 Coeff_Effective       =reshape(database(:,7),[numy, numx, numz]);% [m^3/s]
 Coeff_Excitation      =reshape(database(:,8),[numy, numx, numz]);% [m^3/s]
 Coeff_IonizationAr    =reshape(database(:,9),[numy, numx, numz]);% [m^3/s]
 Coeff_IonizationArex  =reshape(database(:,10),[numy, numx, numz]);% [m^3/s]

 ne=y(1+N*(index.("e")-1):index.("e")*N);
 nAr=2.5e25;
 nArex=y(1+N*(index.("Ar*")-1):index.("Ar*")*N);
 phi=y(1+N*N_species:end);
 E=abs(-gradient(phi,x));%[V/m]
 %E=[0; abs(-diff(phi)./diff(x))];
 reduced_E=E./nAr.*1e21;%[Td] 1Td=1e-21 V*m^2
## if any(reduced_E>200)
##   reduced_E(find(reduced_E>200))=200;
## endif
 ArExFrac=ne./(nAr+nArex);
 coeff_eff            =interp3(EF_Reduced, Plasma_density, Ar_Ex_Fraction, Coeff_Effective, reduced_E, ne , ArExFrac, "linear",1e-16);
 coeff_excitation     =interp3(EF_Reduced, Plasma_density, Ar_Ex_Fraction, Coeff_Excitation, reduced_E, ne , ArExFrac,  "linear",1e-16);
 coeff_ionizationAr   =interp3(EF_Reduced, Plasma_density, Ar_Ex_Fraction, Coeff_IonizationAr, reduced_E, ne , ArExFrac,  "linear",1e-16);
 coeff_ionizationArEx =interp3(EF_Reduced, Plasma_density, Ar_Ex_Fraction, Coeff_IonizationArex, reduced_E, ne , ArExFrac,  "linear",1e-16);
 mean_energy          =interp3(EF_Reduced, Plasma_density, Ar_Ex_Fraction, Mean_Energy, reduced_E, ne , ArExFrac,  "linear",1e-16);

 Te=(2/3)*mean_energy*1.60218e-19/k_b;
 reactions(1).rate_coeffs(:,1)=coeff_ionizationAr;
 reactions(2).rate_coeffs(:,1)=coeff_excitation;
 reactions(3).rate_coeffs(:,1)=coeff_ionizationArEx;
 reactions(5).rate_coeffs(:,1)=7*1e-13*(300./Te).^0.5;
 reactions(8).rate_coeffs(:,1)=coeff_eff;

endfunction
