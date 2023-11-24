function Te=compute_temperature( y, x, index,N_species, database)
 N=numel(x);

 L=x(end);
 k_b     = 1.380649e-23; %Boltzmann constant [J K^-1]


 numx = numel(unique(database(:,1)));
 numy = numel(unique(database(:,3)));
 numz = numel(unique(database(:,4)));

 EF_Reduced            =reshape(database(:,1),[numy, numx, numz]); %Electric_field/N(IN) (Td)
 Plasma_density        =reshape(database(:,3),[numy, numx, numz]); %Plasma Density(IN) [1/m^3]
 Ar_Ex_Fraction        =reshape(database(:,4),[numy, numx, numz]); %Excited Argon Fraction(IN) [1]
 Mean_Energy           =reshape(database(:,6),[numy, numx, numz]); %Electron Mean Energy [eV]

 ne=y(1+N*(index.("e")-1):index.("e")*N);
 nAr=2.5e25;
 nArex=y(1+N*(index.("Ar*")-1):index.("Ar*")*N);
 phi=y(1+end-N:end);

 E=abs(-gradient(phi,x));%[V/m]
 %E=[0; abs(-diff(phi)./diff(x))];

 reduced_E=E./nAr.*1e21;%[Td] 1Td=1e-21 V*m^2
##  if any(reduced_E>200)
##   reduced_E(find(reduced_E>200))=200;
##  endif
 ArExFrac=ne./(nAr+nArex);
 mean_energy          =interp3(EF_Reduced, Plasma_density, Ar_Ex_Fraction, Mean_Energy, reduced_E, ne , ArExFrac,  "linear",15);

 Te=(2/3)*mean_energy*1.60218e-19/k_b;
endfunction
