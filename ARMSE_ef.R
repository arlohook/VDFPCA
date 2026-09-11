ARMSE_ef = function(eig_df, nK = 2){
  
  unlist(sapply(1:nK, function(k){
    
    
    c(left_join(eval.EF.grid(k = k, m_grid = unique(eig_df$m)), filter(eig_df, pc == k), by = c("pc", "m", "t")) %>%
        group_by(m) %>%
        mutate(flip = ifelse(sum(value.x*value.y) < 0, -1, 1)) %>%
        ungroup() %>%
        mutate(value.y = value.y*flip) %>% 
        select(-c(flip)) %>% summarise(ARMSE = sum((value.x-value.y)^2)/101/100))
  }))
  
}