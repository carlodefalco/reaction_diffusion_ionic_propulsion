clc
clear
close all
n_ar=2.5e25; %reference ARGON density

E=linspace(1, 1000, 10); % Reduced electric field array(Td)
plasma_n=[0 logspace(1, 25, 25)]; %electron number density Array span from 0 to 0.001*n_ar
%plasma_n=linspace(0, 2.5e25,10);
Ar_ex_frac=[0 logspace(-18,0,19)]; %fraction of ecxtided ARGON
%n_argon is assumed constant, so that its fraction is one. at most. the TOTAL GAS COMPOSITION=2, considering
% the variation of exctided argon
[X, Y, Z]=meshgrid(E, plasma_n, Ar_ex_frac); %genereate the possible combination of the three parameters

%adjust X,Y,Z in column array
l=numel(E)*numel(plasma_n)*numel(Ar_ex_frac);
x=reshape(X,l,1);
y=reshape(Y,l,1);
z=reshape(Z,l,1);
%compute the ionization_degree
Ar_ex=z.*n_ar;
n_tot=Ar_ex+ones(numel(Ar_ex),1)*n_ar;
ionization_degree= y./n_tot;

%BOLSIG SIMULATION
maxRUN=1000; %max number of successive simulation that BOLSIG+ can perform
nRUN=ceil(numel(x)/maxRUN);
##for ii=1:nRUN
##  %create a .txt file with 1000 variables value at most
##  bolsigVAR=sprintf("bolsigVAR%d.txt", ii)
##  fid = fopen( bolsigVAR , 'wt' );
##  %if the number of residual values of VAR is less then 1000 the file txt is built accounting for it
##  if (numel(x)-ii*maxRUN>=0)
##    end_run=ii*maxRUN;
##  else
##    end_run=numel(x);
##  endif
##    fprintf(fid, '%7.0f. %14.0e %24.0e %28.0e \n',...
##            [x(1 +(ii-1)*maxRUN:end_run),ionization_degree(1 +(ii-1)*maxRUN:end_run),...
##             y(1 +(ii-1)*maxRUN:end_run),z(1 +(ii-1)*maxRUN:end_run)]');
##
##  %put the values of VAR in the .dat input file for BOLSIG
##  input_empty=fileread("inputARGON_empty.dat"); %read the precompiled input in which some parts needs to be fullfilled
##  run_cond=fileread(bolsigVAR);
##  outputname=sprintf("output%d.dat",ii);
##  input_file_new=strrep(input_empty, "/putVARhere", run_cond);
##  input_file_new=strrep(input_file_new, "namefile.dat", outputname);
##  input_full=sprintf("inputARGON%d.dat", ii)
##  fileid=fopen(input_full, "w+");
##  fprintf(fileid, input_file_new);
##  fclose(fileid);
##  cmd=sprintf("echo inputARGON%d.dat|bolsigminus.exe ",ii);
##  system(cmd)
##endfor

%write the output in a single file
Electric_field        =cell(numel(x),1);
Ionization_Degree     =cell(numel(x),1);
Plasma_Density        =cell(numel(x),1);
Ar_Ex_Fraction        =cell(numel(x),1);
Mobility              =cell(numel(x),1);
Mean_Energy           =cell(numel(x),1);
Coeff_Effective       =cell(numel(x),1);
Coeff_Excitation      =cell(numel(x),1);
Coeff_IonizationAr    =cell(numel(x),1);
Coeff_IonizationArex  =cell(numel(x),1);

for ii=1:nRUN
  str=sprintf("output%d.dat",ii);
  input=fileread(str);
  %CONDITIONS
  if (numel(x)-ii*maxRUN>=0)
    end_run=ii*maxRUN;
  else
    end_run=numel(x);
  endif
  [~,~,~,~, EF]    = regexp (input, 'Electric field / N \(Td\)[ ]*([0-9]+\.*[0-9]*)');
  [~,~,~,~,ID]     = regexp (input, 'Ionization degree[ ]*([+-]?([0-9]+\.*[0-9]*[eE]?[+-]?[0-9]+)*)'); %matches 0.000 0.4000E-24, 100, 100.32, 0.1000E+05
  [~,~,~,~,PD]     = regexp (input, 'Plasma density \(1/m3\)[ ]*([+-]?([0-9]+\.*[0-9]*[eE]?[+-]?[0-9]+)*)');
  [~,~,~,~,ArExF]  = regexp (input, 'Mole fraction Ar\*[ ]*([+-]?([0-9]+\.*[0-9]*[eE]?[+-]?[0-9]+)*)');

  %OUTPUTS
  [~,~,~,~,M]      = regexp (input, 'Mobility \*N \(1/m/V/s\)[ ]*([+-]?([0-9]+\.*[0-9]*[eE]?[+-]?[0-9]+)*)');
  [~,~,~,~,ME]     = regexp (input, 'Mean energy \(eV\)[ ]*([0-9]+\.*[0-9]*)');
  [~,~,~,~,CEf]    = regexp (input, 'C1[ ]*Ar[ ]*Effective[ ]*\(momentum\)[ ]*([+-]?([0-9]+\.*[0-9]*[eE]?[+-]?[0-9]+)*)');
  [~,~,~,~,CEx]    = regexp (input, 'C2[ ]*Ar[ ]*Excitation[ ]*11.50[ ]*eV[ ]*([+-]?([0-9]+\.*[0-9]*[eE]?[+-]?[0-9]+)*)');
  [~,~,~,~,CIAr]   = regexp (input, 'C3[ ]*Ar[ ]*Ionization[ ]*15.80[ ]*eV[ ]*([+-]?([0-9]+\.*[0-9]*[eE]?[+-]?[0-9]+)*)');
  [~,~,~,~,CIArex] = regexp (input, 'C4[ ]*Ar\*[ ]*Ionization[ ]*4.20[ ]*eV[ ]*([+-]?([0-9]+\.*[0-9]*[eE]?[+-]?[0-9]+)*)');

  Electric_field(1+(ii-1)*maxRUN:end_run)       =EF;
  Ionization_Degree(1+(ii-1)*maxRUN:end_run)    =ID;
  Plasma_Density (1+(ii-1)*maxRUN:end_run)      =PD;
  Ar_Ex_Fraction(1+(ii-1)*maxRUN:end_run)       =ArExF;
  Mobility(1+(ii-1)*maxRUN:end_run)             =M;
  Mean_Energy(1+(ii-1)*maxRUN:end_run)          =ME;
  Coeff_Effective (1+(ii-1)*maxRUN:end_run)     =CEf;
  Coeff_Excitation(1+(ii-1)*maxRUN:end_run)     =CEx;
  Coeff_IonizationAr(1+(ii-1)*maxRUN:end_run)   =CIAr;
  Coeff_IonizationArex(1+(ii-1)*maxRUN:end_run) =CIArex;
endfor

Electric_field        =(cell2mat(Electric_field));
Ionization_Degree     =(cell2mat(Ionization_Degree));
Plasma_Density        =cell2mat(Plasma_Density);
Ar_Ex_Fraction        =cell2mat(Ar_Ex_Fraction);
Mobility              =cell2mat(Mobility);
Mean_Energy           =cell2mat(Mean_Energy);
Coeff_Effective       =cell2mat(Coeff_Effective);
Coeff_Excitation      =cell2mat(Coeff_Excitation);
Coeff_IonizationAr    =cell2mat(Coeff_IonizationAr);
Coeff_IonizationArex  =cell2mat(Coeff_IonizationArex);
a=[Electric_field, Ionization_Degree, Plasma_Density, Ar_Ex_Fraction,...
  Mobility, Mean_Energy, Coeff_Effective, Coeff_Excitation,...
  Coeff_IonizationAr, Coeff_IonizationArex];
final_output=fopen("database.txt", "w+");
fprintf(final_output, "%%database proprieta' calcolate con BOLSIG+ \n");
fprintf(final_output, "%%Electric_field/N(IN) [Td] Ionization_Degree(IN) [1]...
       Plasma Density(IN) [1/m^3] Excited Argon Fraction(IN) [1] Mobility*N(OUT) [1/m/V/s]...
       Mean_Energy(OUT) [1/m/V/s], Coeff_Effective(OUT) [m^3/s] Coeff_Excitation(OUT) [m^3/s]...
       Coeff_IonizationAr(OUT) [m^3/s] Coeff_IonizationArex(OUT) [m^3/s] \n");
##fprintf(final_output,'%1s %26s %53s %70s %87s %104s %121s %138s %155s %172s \n',Electric_field{:},...
##        Ionization_Degree{:},Plasma_Density{:}, Ar_Ex_Fraction{:},...
##        Mobility{:}, Mean_Energy{:}, Coeff_Effective{:},Coeff_Excitation{:},...
##        Coeff_IonizationAr{:}, Coeff_IonizationArex{:});
for ii=1:numel(x)
  fprintf(final_output,'%1s %26s %30s %30s %30s %30s %30s %35s %35s %35s  \n',a{ii, :});
endfor

fclose(final_output);
