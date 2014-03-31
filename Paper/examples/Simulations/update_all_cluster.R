for (i in 1:length(cl)){
        
    system(paste('ssh',cl[i],'FASTIN-R/git pull'))

}
