function mobility=compute_electron_mobility( y, x,index,N_species, database)
 N=numel(x);
 N_species-=1;
 L=x(end);

 numx = numel(unique(database(:,1)));
 numy = numel(unique(database(:,3)));
 numz = numel(unique(database(:,4)));

 EF_Reduced            =reshape(database(:,1),[numy, numx, numz]); %Electric_field/N(IN) (Td)
 Plasma_density        =reshape(database(:,3),[numy, numx, numz]);% Plasma Density(IN) [1/m^3]
 Ar_Ex_Fraction        =reshape(database(:,4),[numy, numx, numz]);%Excited Argon Fraction(IN) [1]
 %Ionization_Degree     =reshape(database(:,2),[numy, numx, numz]);% Ionization_Degree(IN) [1]
 Mobility_Reduced      =reshape(database(:,5),[numy, numx, numz]);%Mobility*N(OUT) [1/m/V/s]

 ne=y(1+N*(index.("e")-1):index.("e")*N);
 nAr=2.5e25;
 nArex=y(1+N*(index.("Ar*")-1):index.("Ar*")*N);
 phi=y(1+N*N_species:end);
 E=abs(-gradient(phi,x));%[V/m]
 %E=[0; abs(-diff(phi)./diff(x))];
 reduced_E=E./nAr.*1e21;%[Td] 1Td=1e-21 V*m^2
##  if any(reduced_E>200)
##   reduced_E(find(reduced_E>200))=200;
## endif
 ArExFrac=ne./(nAr+nArex);
 mobility=interp3(EF_Reduced, Plasma_density, Ar_Ex_Fraction, Mobility_Reduced, reduced_E, ne , ArExFrac,  "linear",1e21)./nAr;

 mobility = 2 ./ (1./mobility(1:end-1) + 1./mobility(2:end));

 endfunction
