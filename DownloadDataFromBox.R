library("boxr")
#library("devtools")
#devtools::install_github("brendan-r/boxr")
## follow this: https://ucsf.app.box.com/developers/console
box_auth()
box_ls()
box_setwd(54040178619)
files = as.data.frame(box_ls())$id
for(file in files) {
  print(file)
  box_dl(file_id = file)
}

