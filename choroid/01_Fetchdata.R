# Fetch study data from GEO database for Age-related macular degeneration 

# Load Library ------------------------------------------------------------


suppressWarnings({
  lapply(
    c(
      "dplyr",
      "Seurat",
      "HGNChelper",
      'GEOquery',
      'Matrix',
      'gridExtra',
      'tidyverse',
      'ggplot2'
    ),
    library,
    character.only = T
  )
})

# Fetch data  -------------------------------------------------------------

# set working directory 
workdir= "~/choroid" 
setwd(workdir)
accession_no = "GSE188280"
# download raw file to the directory
getGEOSuppFiles(accession_no)
# To open the tar file, access the folder with
paths=paste0(workdir,"/",accession_no)
setwd(paths)
# Open the tar file 
untar(paste0(accession_no,"_RAW.tar"))


# Create a vector list of all the samples ---------------------------------



#retina_sample = c("GSM5676874","GSM5676876","GSM5676878","GSM5676880","GSM5676882","GSM5676884")
choriod_sample = c("GSM5676873","GSM5676875","GSM5676877","GSM5676879","GSM5676881","GSM5676883")

samples = choriod_sample
# Create sub folder  for each sample 


# File manipulation to prepare for seurat object creation -----------------



for (j in 1:length(samples)){
  folder<-dir.create(paste0(paths,'/',samples[j]))
}
# Move each sample file to the folder 

library(filesstrings)
dir.list= list.dirs(paths,recursive = FALSE)


for (i in 1:length(samples)){
  file_list<-list.files(paths, pattern  = paste0(samples[i],'_'))
  file.move(file_list, dir.list[i], overwrite=TRUE)
  
}

# Rename files in directory to fit into 10X format
for (i in dir.list){
  x=list.files(i, pattern = '.gz')
  for (a in x){
    b=gsub(".*_", "", a)
    setwd(i)
    file.rename(a,b)
  }
  
}

setwd(paths)

# Merge all together to prepare seurat object 

# get data location
dirs <- list.dirs(path = paths, recursive = F, full.names = F)

for(x in dirs){
  name <- x
  
  cts <- Read10X(data.dir = x)
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts, min.cells = 3, min.genes =200,
                                  project = 'AMD_choriod'))
}




# Merge dataset -----------------------------------------------------------



merge_max = merge(GSM5676873, y=c(GSM5676875,GSM5676877,GSM5676879,GSM5676881,
                                  GSM5676883),
                    add.cell.ids=c("GSM5676873", "GSM5676875","GSM5676877",
                                   "GSM5676879","GSM5676881","GSM5676883"),
                   project = 'AMD_choroid')


# # Save seurat object in file
#
save(merge_max,
      file=paste0(workdir,'/','choroid_cell.RData'))




