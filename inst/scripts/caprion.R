HOME <- Sys.getenv("HOME")
load(paste(HOME,"Caprion","caprion.rda",sep="/"))
caprion <- protein_list
save(caprion,file='caprion.rda',compress='xz')
