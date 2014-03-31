for (i in 1:length(cl)){
    
    system(paste('scp FASTIN_cluster_install.sh philipp@',cl[i],':/home/philipp/',sep=''))
    
    system(paste('ssh',cl[i],'chmod 700 FASTIN_cluster_install.sh'))

    system(paste('ssh',cl[i],'./FASTIN_cluster_install.sh'))

}
