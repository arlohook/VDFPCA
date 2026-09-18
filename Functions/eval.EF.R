eval.EF = function(k, t, m, case = 1){
  
  if(case == 1){
    W = pnorm(m, 30, 10)
    
    out = W*(sqrt(2)/sqrt(m))*sin(2*k*pi*t/m) + (1-W)*(sqrt(2)/sqrt(m))*cos(2*pi*t/m)
  }
  
  if(case == 2){
    
    
    Wm = pnorm(m, mean = 40, sd = 30)
    alpha = pi * sin(8*Wm^3)
    
    scale = sqrt(2) / sqrt(m)
    
    out = scale * sin(2 * k * pi * t / m + alpha) 
    
  }
  
  out
}