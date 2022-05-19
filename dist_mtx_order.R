# Function to order distance matrices
# Code source https://github.com/leylabmpi/animal_gut_16S-uni

dist_mtx_order = function(d, x){
  # Ordering distance matrixes
  # d = distance matrix (dist class)
  # x = vector to order dist. obj. by
  m = d %>% as.matrix
  d = as.dist(m[x,x])
  return(d)
}
