
for(i in c(3,1,2)){
  f <- list.files(path = "/data/4C/high_res_methods/new_analysis/reciprocal/peakC/", recursive = T, pattern=paste0("pos",i), full.names = T)
  p <- readMultiple(f[1], vp.pos = vp[i], window = 1.5e6)
  r <- single.analysis(p,  vp.pos=vp[i], qWd = 2.5)
  plot_C(r, xlim=c(64e6,67e6))
}
