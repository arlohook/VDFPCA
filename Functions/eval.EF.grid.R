eval.EF.grid = function(k, m_grid, case, Tn = 101){
  
  do.call(rbind, lapply(1:length(m_grid), function(i){
    
    data.frame("pc" = k,
               "m" = unique(m_grid)[i], 
               "t" = seq(0,unique(m_grid)[i], length = Tn),
               "value" = eval.EF(k = k, 
                                 t = seq(0,unique(m_grid)[i], length = Tn), 
                                 m = unique(m_grid)[i],
                                 case = case))
    
  }))
  
}