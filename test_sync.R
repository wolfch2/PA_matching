library(dplyr)
library(stringr)
library(aws.s3) 
library(purrr)

# aws.s3::s3sync(path = "data_input",
#                bucket = "fbedecarrats",
#                prefix = "Replication_wolf", # diffusion to be able to share
#                create = FALSE,
#                region = "")

my_files <- get_bucket_df(bucket = "fbedecarrats",
                          prefix = "Replication_wolf",
                          region = "")

my_tifs <- my_files %>%
  filter(str_ends(Key, ".tif")) %>%
  pluck("Key")

my_tifs_dest <- str_replace(my_tifs, "Replication_wolf", "data_input")

save_from_s3 <- function(x, y) {
  aws.s3::save_object(
    object = x,
    bucket = "fbedecarrats",
    file = y,
    overwrite = FALSE,
    region = "")
}

map2(my_tifs, my_tifs_dest, save_from_s3)


save_object(
  file = "data_input/other.zip",
  object = "Replication_wolf/other2.zip",
  bucket = "fbedecarrats",
  region = "")
