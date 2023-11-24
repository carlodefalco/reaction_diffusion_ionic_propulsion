function current=compute_flux( y,q,Vth,mobility,N,L)


phi=y(end-N+1:end);
n=y(1:N);
[bp, bn]=bimu_bernoulli((phi(2:end)-phi(1:end-1))/Vth);
current=(mobility*Vth/L)*(n(2:end).*bp-n(1:end-1).*bn);
endfunction
