# sudo apt-get update
# sudo apt-get install libcurl4-openssl-dev
# sudo apt-get install libssl-dev
install.packages(c("dplyr", "stringr", "aws.s3", "purrr"))

library(dplyr)
library(stringr)
library(aws.s3)
library(purrr)

# aws.s3::s3sync(path = "data_input",
#                bucket = "fbedecarrats",
#                prefix = "Replication_wolf", # diffusion to be able to share
#                create = FALSE,
#                region = "")

save_object(
  file = "other.zip",
  object = "Replication_wolf/other2.zip",
  bucket = "fbedecarrats",
  region = "")

# unzip other.zip

file.remove("other.zip")

my_files <- get_bucket_df(bucket = "fbedecarrats",
                          prefix = "Replication_wolf",
                          region = "")

my_tifs <- my_files %>%
  filter(str_ends(Key, ".tif")) %>%
  pluck("Key")

my_tifs_dest <- str_replace(my_tifs, "Replication_wolf", "fil_input")

save_from_s3 <- function(x, y) {
  aws.s3::save_object(
    object = x,
    bucket = "fbedecarrats",
    file = y,
    overwrite = FALSE,
    region = "")
}

map2(my_tifs, my_tifs_dest, save_from_s3)

file.info("data_processed/threatened.tif")

put_object(
  file = "data_processed/threatened.tif",
  object = "Replication_wolf/threatened.tif",
  bucket = "fbedecarrats",
  region = "")