eval.EF = function(k, t, m){
  
  W = pnorm(m, 30, 10)
  
  out = W*(sqrt(2)/sqrt(m))*sin(2*k*pi*t/m) + (1-W)*(sqrt(2)/sqrt(m))*cos(2*pi*t/m)
  
  
  out
}