select_markers_wpp <- function(norm.dat, cl, n.markers=20,de.genes=NULL, ...)                           
{
  if(is.null(de.genes)){
    print("de_stats_all_pairs_wpp")
    de.genes=de_stats_all_pairs_wpp(norm.dat, cl, ...)
  }
  pairs = names(de.genes)
  pairs.df = gsub("cl","", do.call("rbind",strsplit(pairs, "_")))
  row.names(pairs.df)=pairs
  select.pairs = pairs[pairs.df[,1] %in% cl & pairs.df[,2]%in% cl]
  de.markers = sapply(select.pairs, function(s){
    tmp = de.genes[[s]]
    c(head(tmp$up.genes,n.markers), head(tmp$down.genes,n.markers))
  },simplify=F)
  markers = intersect(unlist(de.markers),row.names(norm.dat))
  return(list(markers=markers, de.genes=de.genes[select.pairs]))
}