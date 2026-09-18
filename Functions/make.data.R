make.data = function(N, e.sig, case){
  
  Mu = function(m){
    t = seq(0,m, length = 101)
    0.0001*(t-120)^2 + 3*sin(pi*t/60)
  }
  
  if(case ==1){
    Ef = function(m, k){
      
      t = seq(0,m, length = 101)
      W = pnorm(m, mean = 30, sd = 10)
      W*(sqrt(2)/sqrt(m))*sin(2*k*pi*t/m) + (1-W)*(sqrt(2)/sqrt(m))*cos(2*k*pi*t/m)
      
    }
    
    M = c(2, 80)
    
    mi = runif(n = N, min = M[1], M[2])
    
    sigma = diag(0.5^(0:9))
    
  }
  
  if(case == 2){
    Ef = function(m, k){
      
      t = seq(0,m, length = 101)
      Wm = pnorm(m, mean = 40, sd = 30)
      alpha = pi * sin(8*Wm^3)
      
      scale = sqrt(2) / sqrt(m)
      
      out = scale * sin(2 * k * pi * t / m + alpha)
      
    }
    
    
    
    mi = rgamma(n = N, shape = 50, rate = 1)
    
    sigma = diag(0.8^(0:9))
    
  }
  
  
  Scores = rmvnorm(N, sigma = sigma)
  
  
  X = sapply(1:N, function(i){
    
    xi = Mu(mi[i])+ Scores[i,1]*Ef(mi[i], k = 1)+
      Scores[i,1]*Ef(mi[i], k = 1) +
      Scores[i,2]*Ef(mi[i], k = 2) + 
      Scores[i,3]*Ef(mi[i], k = 3) +
      Scores[i,4]*Ef(mi[i], k = 4) +
      Scores[i,5]*Ef(mi[i], k = 5) +
      Scores[i,6]*Ef(mi[i], k = 6) + 
      Scores[i,7]*Ef(mi[i], k = 7) +
      Scores[i,8]*Ef(mi[i], k = 8) +
      Scores[i,9]*Ef(mi[i], k = 9) +
      Scores[i,10]*Ef(mi[i], k = 10)
  })
  
  X = X + matrix(rnorm(N*101, 0, e.sig), 101, N)
  
  out = list("X" = X, "mi" = mi)
}
