library(tidyverse)
setwd("/home/luomeng/DATA/BRAIN/SnRNA/")

cellAlignment.files <-
  list.files(
    "/home/luomeng/DATA/REMOTE_Broken/ST_results_lb/cell_deconvolution/sc2st/",
    pattern = ".csv"
  )

library(doParallel)
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)


cell_spots <-
  foreach(filename = cellAlignment.files, .packages = c("tidyverse")) %dopar% {
    read_csv(
      file.path(
        "/home/luomeng/DATA/REMOTE_Broken/ST_results_lb/cell_deconvolution/sc2st",
        filename
      )
    )
  }
parallel::stopCluster(cl)

cl <- parallel::makeCluster(50)
doParallel::registerDoParallel(cl)

for (i in 1:length(cell_spots)) {
  df <- cell_spots[[i]]
  n <- 5000
  nr <- nrow(df)
  df %<>% column_to_rownames("...1")
  df.list <-
    split(df, rep(1:ceiling(nr / n), each = n, length.out = nr))
  timestamp()
  sample = str_remove_all(cellAlignment.files[[i]], "_sc2st.csv")
  
  cells_to_spots <-
    foreach(
      data = df.list,
      .packages = c("tidyverse"),
      .combine = "rbind"
    ) %dopar% {
      cells_list = rownames(data)
      spot_list = colnames(data)
      belonging_list <-
        purrr::map_dfr(1:length(cells_list), function(i) {
          max_prob = max(data[i,])
          data.frame(
            cellID = cells_list[i],
            max_prob = max_prob,
            spotID =
              spot_list[data[i,] >= max_prob],
            sample = sample
          )
        })
    }
  
  saveRDS(cells_to_spots,str_c(sample,"_cellAlignment.rds"))
  print(timestamp())
}
  
  parallel::stopCluster(cl)
  