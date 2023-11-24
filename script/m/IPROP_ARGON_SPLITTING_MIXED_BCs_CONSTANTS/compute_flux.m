function current=compute_flux( n,phi,q,Vth,mobility,x)

  dx=diff(x)(end);

 [bp, bn]=bimu_bernoulli((phi(2:end)-phi(1:end-1))/Vth);%);
  current=(q*mobility*Vth./dx).*(n(2:end).*bp-n(1:end-1).*bn);

endfunction
