##### Load packages #####

pkgTest <- function(x)
{
  if (!require(x, character.only = TRUE))
  {
    install.packages(x, dep = TRUE)
    if (!require(x, character.only = TRUE))
      stop("Package not found")
  }
}

# pkgTest("shinyLP")
pkgTest("shinyBS")
pkgTest("seqinr")
pkgTest("plotly")
pkgTest("DT")
pkgTest("shiny")
pkgTest("gtools")
pkgTest("plyr")

map <- function(v, mu) {
  mapply(
    x = gsub("(^\\d+)(.*$)", "\\1", v),
    y = strsplit(gsub("\\d+", "", v), "/"),
    FUN = function(x, y) {
      xy <- paste(x, y, sep = "")
      m <- sapply(xy, FUN = grepl, mu)
      apply(m, MARGIN = 2, any)
    }
  )
}

library(plyr)

# get all the zip files
zipF <- list.files(path = paste0(dirname(rstudioapi::getSourceEditorContext()$path)), pattern = "*.zip", full.names = TRUE)

# unzip all your files
ldply(.data = zipF, .fun = unzip, exdir = paste0(dirname(rstudioapi::getSourceEditorContext()$path)))

options(shiny.trace=TRUE) #show what happens in shiny

server <- function(input, output){
  
  #run Segminator II
  observeEvent(input$do, {
    setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/segminator_II"))
    system("java -jar segminator_II_0.1.1.jar")
  })
  
  datasetInput <- eventReactive(input$Load, {
    
    if(input$Load == 0){return()}
    inFile <- input$file1
    if (is.null(inFile)){return(NULL)}
    
    my_data <- read.csv(inFile$datapath, header = T,sep = "\t", quote = input$quote, stringsAsFactors = FALSE)
    my_data <- data.frame(my_data)
    t <- my_data[,1:8] #remove unnecessary columns
    pos <- 1807:5096 #pol gene positions
    tab2 <- cbind(pos, t) #aggregating positions
    
    
    if (input$region == "Protease"){
      
      #PR 
      pr1 <- tab2[pos >= 2253 & pos <= 2549, ] 
      pr1[,4:7]/pr1[,3]
      pr3 <- cbind(pr1, pr1[,4:7]/pr1[,3])
      pr4 <- cbind(pr3[,1], pr3[,9:13], pr3[,10:13]) #position, reference, frequency
      colnames(pr4) <- c("pos", "ref", "freqa", "freqt", "freqg", "freqc", "a", "t", "g", "c")
      
      pr4$a[which(pr4$a!=0)] <- "A"
      pr4$t[which(pr4$t!=0)] <- "T"
      pr4$g[which(pr4$g!=0)] <- "G"
      pr4$c[which(pr4$c!=0)] <- "C"
      
      datafim <- data.frame() #dataframe to store the results
      
      i <- 1
      while (i <= (dim(pr4)[1])-2){ #length protease
        d <- pr4[i:(i+2),] #dataframe with one codon 
        aap <- (d$pos[3]-2252)/3 #amino acids positions
        aar <- translate(s2c(paste(d$ref, collapse = ""))) #amino acids reference
        
        #assembling nucleotides
        dn <- data.frame(expand.grid(as.character(d[1,7:10]), as.character(d[2,7:10]), as.character(d[3,7:10]))) #combinations of codons
        df <- data.frame(expand.grid(as.numeric(d[1,3:6]), as.numeric(d[2,3:6]), as.numeric(d[3,3:6]))) #combinations of frequencies
        dd <- data.frame(cbind(dn, df))
        ddind <- apply(dd, 1, function(row) all(row !='0' )) #remove rows with zero
        dd2 <- dd[ddind,] #new dataframe with codons and frequencies
        colnames(dd2) <- c("pos1", "pos2", "pos3", "freq1", "freq2", "freq3") 
        
        if (dim(dd2)[1]!=0){
          for (ff in 1:dim(dd2)[1]){
            dd2$pos[ff] <- aap
            dd2$codonref[ff] <- paste(d$ref, collapse = "")
            dd2$aaref[ff] <- aar
            dd2$codon[ff] <- paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = "")
            dd2$aa[ff] <- translate(s2c(paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = ""))) #codification in amino acids
            dd2$freqf[ff] <- dd2$freq1[ff] * dd2$freq2[ff] * dd2$freq3[ff]
          }
        } else {
          ddd2 <- data.frame()
        }
        datafim <- rbind(datafim, dd2)
        i <- i+3
      }
      
      #Analyzing amino acids with frequency greater than or equal to 0.01%
      
      data <- datafim[,7:12] #limiting the number of columns
      data$freqt <- data$freqf*100 #frequency to percentage
      data2 <- data[data$freqt >=0.01,] #removing amino acids frequency <= 0.01%
      data3 <- data2
      data3$codon <- NULL; data3$codonref <- NULL #dataframe with only important columns 
      data3 <- subset(data3, !is.na(freqf)) #removing rows with NA
      
      tt <- data.frame()
      
      for (r in 1:as.numeric(data3$pos[length(data3$pos)])){
        temp <- data3[which(data3$pos==r),]
        repaa <- names(which((table(temp$aa)>=2)==TRUE)) #find repeated amino acids
        ll <- length(repaa)
        for (x in 1:ll) {
          ff <- sum(temp[temp$aa==repaa[x],]$freqf) #sum of equal amino acids frequencies
          temp$freqf[which(temp$aa==repaa[x])[1]] <- ff
          temp2 <- temp[-which(temp$aa==repaa[x])[-1],] 
          temp <- temp2
        }
        tt <- rbind(tt, temp)
      }
      
      tabpos <- unique(tt$pos)
      datapos <- unique(data3$pos)
      falt <- datapos[!datapos%in%tabpos] #positions that do not appear in table because they do not have repeated amino acids
      complem <- data3[data3$pos%in%falt,]
      tabela <- rbind(tt, complem)
      tabela <- tabela[order(tabela$pos, -tabela$freqf), ]
      
      #coverage
      covmat <- as.data.frame(matrix(pr1$cover, ncol = 3, byrow = T)) #coverage
      #taking the lowest coverage of each codon
      min_covmat <- apply(covmat,1,min)
      #repeating each coverage the number of times that each position appears in the database datafim
      vec_cov <- NULL
      for(iii in 1: length(tabela$pos)){
        tamanho <- length(which(tabela$pos == iii))
        repete <- rep(min_covmat[iii], tamanho)
        vec_cov <- c(vec_cov, repete)
      }
      
      #table with the cover (we used the minimum coverage between the 3 nucleotides of each codon)
      tabela <- as.data.frame(cbind(tabela[,1:4], vec_cov))
      posmutpr <- c(10,11,13,16,20,23,24,30,32,33,34,35,36,43,46,47,48,50,53,54,58,60,62,63,64,69,71,73,74,76,77,82,83,84,85,88,89,90,93)
      dposmut <- tabela[tabela$pos%in%posmutpr,]
      
      aarr <- c('L','V','I','G','K','L','L','D','V','L','E','E','M','K','M','I','G','I','F','I','Q','D','I','L','I','H','A','G','T','L','V','V','N','I','I','N','L','L','I')
      
      #amino acid reference subtype B
      for(m in 1:length(aarr)){
        for(l in 1:dim(dposmut)[1]){
          if (dposmut[l,1] == posmutpr[m]){
            dposmut[l,2] <- aarr[m]
          }
        }
      }
      
      #names(dposmut)[5] <- "pr"
      dposmut <- dposmut[dposmut$freqf >= 0.01,]
      dposmut <- as.data.frame(cbind(dposmut[,5], dposmut[,2], dposmut[,1], dposmut[,3], dposmut[,4]))
      colnames(dposmut) <- c("Coverage", "Amino Acid in Protease", "Position", "Amino Acid", "Frequency")
      dposmut[,5] <- as.numeric(as.character(dposmut[,5]))
    } 
    
    if (input$region == "Reverse Transcriptase"){
      
      #RT
      
      rt1 <- tab2[pos >= 2550 & pos <= 3869, ]
      rt1[,4:7]/rt1[,3]
      rt3 <- cbind(rt1, rt1[,4:7]/rt1[,3])
      rt4 <- cbind(rt3[,1], rt3[,9:13], rt3[,10:13]) #position, reference, frequency
      colnames(rt4) <- c("pos", "ref", "freqa", "freqt", "freqg", "freqc", "a", "t", "g", "c")
      
      rt4$a[which(rt4$a!=0)] <- "A"
      rt4$t[which(rt4$t!=0)] <- "T"
      rt4$g[which(rt4$g!=0)] <- "G"
      rt4$c[which(rt4$c!=0)] <- "C"
      
      datafim <- data.frame() #dataframe to store the results
      
      i <- 1
      while (i <= (dim(rt4)[1])-2){ #length transcriptase reverse
        d <- rt4[i:(i+2),] #dataframe with one codon 
        aap <- (d$pos[3]-2549)/3 #amino acids positions 
        aar <- translate(s2c(paste(d$ref, collapse = ""))) #amino acids reference
        
        #assembling nucleotides
        dn <- data.frame(expand.grid(as.character(d[1,7:10]), as.character(d[2,7:10]), as.character(d[3,7:10]))) #combinations of codons
        df <- data.frame(expand.grid(as.numeric(d[1,3:6]), as.numeric(d[2,3:6]), as.numeric(d[3,3:6]))) #combinations of frequencies
        dd <- data.frame(cbind(dn, df))
        ddind <- apply(dd, 1, function(row) all(row !='0' )) #remove rows with zero
        dd2 <- dd[ddind,] #new dataframe with codons and frequencies
        colnames(dd2) <- c("pos1", "pos2", "pos3", "freq1", "freq2", "freq3") 
        
        if (dim(dd2)[1]!=0){
          for (ff in 1:dim(dd2)[1]){
            dd2$pos[ff] <- aap
            dd2$codonref[ff] <- paste(d$ref, collapse = "")
            dd2$aaref[ff] <- aar
            dd2$codon[ff] <- paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = "")
            dd2$aa[ff] <- translate(s2c(paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = ""))) #codification in amino acids
            dd2$freqf[ff] <- dd2$freq1[ff] * dd2$freq2[ff] * dd2$freq3[ff]
          }
        } else {
          ddd2 <- data.frame()
        }
        datafim <- rbind(datafim, dd2)
        i <- i+3
      }
      
      #Analyzing amino acids with frequency greater than or equal to 0.01%
      
      data <- datafim[,7:12] #limiting the number of columns
      data$freqt <- data$freqf*100 #frequency to percentage
      data2 <- data[data$freqt >=0.01,] #removing amino acids frequency <= 0.01%
      data3 <- data2
      data3$codon <- NULL; data3$codonref <- NULL #dataframe with only important columns 
      data3 <- subset(data3, !is.na(freqf)) #removing rows with NA
      
      tt <- data.frame()
      
      for (r in 1:as.numeric(data3$pos[length(data3$pos)])){
        temp <- data3[which(data3$pos==r),]
        repaa <- names(which((table(temp$aa)>=2)==TRUE)) #find repeated amino acids
        ll <- length(repaa)
        for (x in 1:ll) {
          ff <- sum(temp[temp$aa==repaa[x],]$freqf) #sum of equal amino acids frequencies
          temp$freqf[which(temp$aa==repaa[x])[1]] <- ff
          temp2 <- temp[-which(temp$aa==repaa[x])[-1],] 
          temp <- temp2
        }
        tt <- rbind(tt, temp)
      }
      
      tabpos <- unique(tt$pos)
      datapos <- unique(data3$pos)
      falt <- datapos[!datapos%in%tabpos] #positions that do not appear in table because they do not have repeated amino acids
      complem <- data3[data3$pos%in%falt,]
      tabela <- rbind(tt, complem)
      tabela <- tabela[order(tabela$pos, -tabela$freqf), ]
      
      #coverage
      covmat <- as.data.frame(matrix(rt1$cover, ncol = 3, byrow = T)) #coverage
      #taking the lowest coverage of each codon
      min_covmat <- apply(covmat,1,min)
      #repeating each coverage the number of times that each position appears in the database datafim
      vec_cov <- NULL
      for(iii in 1: length(tabela$pos)){
        tamanho <- length(which(tabela$pos == iii))
        repete <- rep(min_covmat[iii], tamanho)
        vec_cov <- c(vec_cov, repete)
      }
      
      #table with the coverage (we used the minimum coverage between the 3 nucleotides of each codon)
      tabela <- as.data.frame(cbind(tabela[,1:4], vec_cov))
      posmutrt <- c(40,41,44,62,65,67,68,69,70,74,75,77,90,98,100,101,103,106,108,115,116,118,138,151,179,181,184,188,190,210,215,219,221,225,227,230,234,236,238,318,348)
      dposmut <- tabela[tabela$pos%in%posmutrt,]
      
      aarr_rt <- c('E','M','E','A','K','D','S','T','K','L','V','F','V','A','L','K','K','V','V','Y','F','V','E','Q','V','Y','M','Y','G','L','T','K','H','P','F','M','L','P','K','Y','N')
      
      #amino acids reference subtype B
      for(m in 1:length(aarr_rt)){
        for(l in 1:dim(dposmut)[1]){
          if (dposmut[l,1] == posmutrt[m]){
            dposmut[l,2] <- aarr_rt[m]
          }
        }
      }
      
      # names(dposmut)[5] <- "rt"
      dposmut <- dposmut[dposmut$freqf >= 0.01,]
      dposmut <- as.data.frame(cbind(dposmut[,5], dposmut[,2], dposmut[,1], dposmut[,3], dposmut[,4]))
      colnames(dposmut) <- c("Coverage", "Amino Acid in Reverse Transcriptase", "Position", "Amino Acid", "Frequency (%)")
      dposmut[,5] <- as.numeric(as.character(dposmut[,5]))
      
    } 
    
    if (input$region == "Integrase"){
      
      #IN
      
      int1 <- tab2[pos >= 4230 & pos <= 5096, ] 
      int1[,4:7]/int1[,3]
      int3 <- cbind(int1, int1[,4:7]/int1[,3])
      int4 <- cbind(int3[,1], int3[,9:13], int3[,10:13]) #position, reference, frequency
      colnames(int4) <- c("pos", "ref", "freqa", "freqt", "freqg", "freqc", "a", "t", "g", "c")
      
      int4$a[which(int4$a!=0)] <- "A"
      int4$t[which(int4$t!=0)] <- "T"
      int4$g[which(int4$g!=0)] <- "G"
      int4$c[which(int4$c!=0)] <- "C"
      
      datafim <- data.frame() #dataframe to store the results
      
      i <- 1
      while (i <= (dim(int4)[1])-2){ #length integrase
        d <- int4[i:(i+2),] #dataframe with one codon 
        aap <- (d$pos[3]-4229)/3 #amino acids positions
        aar <- translate(s2c(paste(d$ref, collapse = ""))) #amino acids reference
        
        #assembling nucleotides
        dn <- data.frame(expand.grid(as.character(d[1,7:10]), as.character(d[2,7:10]), as.character(d[3,7:10]))) #combinations of codons
        df <- data.frame(expand.grid(as.numeric(d[1,3:6]), as.numeric(d[2,3:6]), as.numeric(d[3,3:6]))) #combinations of frequencies
        dd <- data.frame(cbind(dn, df))
        ddind <- apply(dd, 1, function(row) all(row !='0' )) #remove rows with zero
        dd2 <- dd[ddind,] #new dataframe with codons and frequencies
        colnames(dd2) <- c("pos1", "pos2", "pos3", "freq1", "freq2", "freq3") 
        
        if (dim(dd2)[1]!=0){
          for (ff in 1:dim(dd2)[1]){
            dd2$pos[ff] <- aap
            dd2$codonref[ff] <- paste(d$ref, collapse = "")
            dd2$aaref[ff] <- aar
            dd2$codon[ff] <- paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = "")
            dd2$aa[ff] <- translate(s2c(paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = ""))) #codification in amino acids
            dd2$freqf[ff] <- dd2$freq1[ff] * dd2$freq2[ff] * dd2$freq3[ff]
          }
        } else {
          ddd2 <- data.frame()
        }
        datafim <- rbind(datafim, dd2)
        i <- i+3
      }
      
      #Analyzing amino acids with frequency greater than or equal to 0.01%
      
      data <- datafim[,7:12] #limiting the number of columns
      data$freqt <- data$freqf*100 #frequency to percentage
      data2 <- data[data$freqt >=0.01,] #removing amino acids frequency <= 0.01%
      data3 <- data2
      data3$codon <- NULL; data3$codonref <- NULL #dataframe with only important columns
      data3 <- subset(data3, !is.na(freqf)) #removing rows with NA
      
      tt <- data.frame()
      
      for (r in 1:as.numeric(data3$pos[length(data3$pos)])){
        temp <- data3[which(data3$pos==r),]
        repaa <- names(which((table(temp$aa)>=2)==TRUE)) #find repeated amino acids
        ll <- length(repaa)
        for (x in 1:ll) {
          ff <- sum(temp[temp$aa==repaa[x],]$freqf) #sum of equal amino acids frequencies
          temp$freqf[which(temp$aa==repaa[x])[1]] <- ff
          temp2 <- temp[-which(temp$aa==repaa[x])[-1],] 
          temp <- temp2
        }
        tt <- rbind(tt, temp)
      }
      
      tabpos <- unique(tt$pos)
      datapos <- unique(data3$pos)
      falt <- datapos[!datapos%in%tabpos] #positions that do not appear in table because they do not have repeated amino acids
      complem <- data3[data3$pos%in%falt,]
      tabela <- rbind(tt, complem)
      tabela <- tabela[order(tabela$pos, -tabela$freqf), ]
      
      #coverage
      covmat <- as.data.frame(matrix(int1$cover, ncol = 3, byrow = T)) #coverage
      #taking the lowest coverage of each codon
      min_covmat <- apply(covmat,1,min)
      #repeating each coverage the number of times that each position appears in the database datafim
      vec_cov <- NULL
      for(iii in 1: length(tabela$pos)){
        tamanho <- length(which(tabela$pos == iii))
        repete <- rep(min_covmat[iii], tamanho)
        vec_cov <- c(vec_cov, repete)
      }
      
      #table with the coverage (we used the minimum coverage between the 3 nucleotides of each codon)
      tabela <- as.data.frame(cbind(tabela[,1:4], vec_cov))
      posmutint <- c(51,66,74,92,95,97,101,114,118,121,128,138,140,143,145,146,147,148,151,153,155,157,163,193,230,263)
      dposmut <- tabela[tabela$pos%in%posmutint,]
      
      aarr_int <- c('H','T','L','E','Q','T','L','H','G','F','A','E','G','Y','P','Q','S','Q','V','S','N','E','G','G','S','R')
      
      #amino acids reference subtype B
      for(m in 1:length(aarr_int)){
        for(l in 1:dim(dposmut)[1]){
          if (dposmut[l,1] == posmutint[m]){
            dposmut[l,2] <- aarr_int[m]
          }
        }
      }
      
      # names(dposmut)[5] <- "int"
      dposmut <- dposmut[dposmut$freqf >= 0.01,]
      dposmut <- as.data.frame(cbind(dposmut[,5], dposmut[,2], dposmut[,1], dposmut[,3], dposmut[,4]))
      colnames(dposmut) <- c("Coverage", "Amino Acid in Integrase", "Position", "Amino Acid", "Frequency (%)")
      dposmut[,5] <- as.numeric(as.character(dposmut[,5]))
    }
    
    dposmut
    
  })
  
  output$table1 <- DT::renderDataTable({server= FALSE
  
  if(input$Load == 0){return()} #just to not show warning messages

  dposmut <- datasetInput()
  
  if (colnames(dposmut)[2] == 'Amino Acid in Protease' | colnames(dposmut)[2] == 'Amino Acid in Reverse Transcriptase' | colnames(dposmut)[2] == 'Amino Acid in Integrase'){
    
    aa_mont <- paste(dposmut[,2], dposmut[,3], dposmut[,4], sep = "")
    arre <- round(as.numeric(dposmut[,5])*100, digits=2)
    pas <- paste(arre, "%", sep = "")
    mmm <- cbind(aa_mont, pas, as.numeric(as.character(dposmut[,1])))
    colnames(mmm) <- c("Amino acids", "Frequency", "Coverage")
    mmm
    
  }
  
  table1 <- datatable(mmm, extensions = 'Buttons', options = list(searching = T, pageLength = nrow(mmm),
                                                                  dom = 'lfBrtip',
                                                                  buttons =
                                                                    list(list(
                                                                      extend = 'collection',
                                                                      buttons = list(list(extend='csv',
                                                                                          filename = paste(input$text, "-", input$region), title = "Drug Resistance Positions"),
                                                                                     list(extend='excel',
                                                                                          filename = paste(input$text, "-", input$region), title = "Drug Resistance Positions"),
                                                                                     list(extend='pdf',
                                                                                          filename= paste(input$text, "-", input$region), title = "Drug Resistance Positions", messageTop = "This table shows the resistance positions, the frequency found for each amino acid, and the coverage.", 
                                                                                          messageBottom = "All sequence interpretations are for research use only."),
                                                                                     list(extend='print',
                                                                                          filename= paste(input$text, "-", input$region), title = "Drug Resistance Positions", messageTop = "This table shows the resistance positions, the frequency found for each amino acid, and the coverage.", 
                                                                                          messageBottom = "All sequence interpretations are for research use only.")),
                                                                      text = 'Download'
                                                                      # customize the length menu
                                                                      
                                                                    )),
                                                                  scrollX = TRUE,
                                                                  initComplete = JS("function(settings, json) {",
                                                                                    "$(this.api().table().header()).css({'background-color': '#2c3e50', 'color': '	#ffffff'});",
                                                                                    "}")))
  
  
  })
  
  datasetInput2 <- eventReactive(input$Load, {
    if(input$Load == 0){return()}
    inFile <- input$file1
    if (is.null(inFile)){return(NULL)}
    
    my_data <- read.csv(inFile$datapath, header = T,sep = "\t", quote = input$quote,stringsAsFactors =FALSE)
    my_data <- data.frame(my_data)
    t <- my_data[,1:8] #remove unnecessary columns
    pos <- 1807:5096 #pol gene positions
    tab2 <- cbind(pos, t) #aggregating positions
    
    if (input$region == "Protease"){
      #PR 
      pr1 <- tab2[pos >= 2253 & pos <= 2549, ] 
      pr1[,4:7]/pr1[,3]
      pr3 <- cbind(pr1, pr1[,4:7]/pr1[,3])
      pr4 <- cbind(pr3[,1], pr3[,9:13], pr3[,10:13]) #position, reference, frequency
      colnames(pr4) <- c("pos", "ref", "freqa", "freqt", "freqg", "freqc", "a", "t", "g", "c")
      
      pr4$a[which(pr4$a!=0)] <- "A"
      pr4$t[which(pr4$t!=0)] <- "T"
      pr4$g[which(pr4$g!=0)] <- "G"
      pr4$c[which(pr4$c!=0)] <- "C"
      
      datafim <- data.frame() #dataframe to store the results
      
      i <- 1
      while (i <= (dim(pr4)[1])-2){ #length protease
        d <- pr4[i:(i+2),] #dataframe with one codon 
        aap <- (d$pos[3]-2252)/3 #amino acids positions
        aar <- translate(s2c(paste(d$ref, collapse = ""))) #amino acids reference
        
        #assembling nucleotides
        dn <- data.frame(expand.grid(as.character(d[1,7:10]), as.character(d[2,7:10]), as.character(d[3,7:10]))) #combinations of codons
        df <- data.frame(expand.grid(as.numeric(d[1,3:6]), as.numeric(d[2,3:6]), as.numeric(d[3,3:6]))) #combinations of frequencies
        dd <- data.frame(cbind(dn, df))
        ddind <- apply(dd, 1, function(row) all(row !='0' )) #remove rows with zero
        dd2 <- dd[ddind,] #new dataframe with codons and frequencies
        colnames(dd2) <- c("pos1", "pos2", "pos3", "freq1", "freq2", "freq3") 
        
        if (dim(dd2)[1]!=0){
          for (ff in 1:dim(dd2)[1]){
            dd2$pos[ff] <- aap
            dd2$codonref[ff] <- paste(d$ref, collapse = "")
            dd2$aaref[ff] <- aar
            dd2$codon[ff] <- paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = "")
            dd2$aa[ff] <- translate(s2c(paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = ""))) #codification in amino acids
            dd2$freqf[ff] <- dd2$freq1[ff] * dd2$freq2[ff] * dd2$freq3[ff]
          }
        } else {
          ddd2 <- data.frame()
        }
        datafim <- rbind(datafim, dd2)
        i <- i+3
      }
      
      #Analyzing amino acids with frequency greater than or equal to 0.01%
      
      data <- datafim[,7:12] #limiting the number of columns
      data$freqt <- data$freqf*100 #frequency to percentage
      data2 <- data[data$freqt >=0.01,] #removing amino acids frequency <= 0.01%
      data3 <- data2
      data3$codon <- NULL; data3$codonref <- NULL #dataframe with only important columns 
      data3 <- subset(data3, !is.na(freqf)) #removing rows with NA
      
      tt <- data.frame()
      
      for (r in 1:as.numeric(data3$pos[length(data3$pos)])){
        temp <- data3[which(data3$pos==r),]
        repaa <- names(which((table(temp$aa)>=2)==TRUE)) #find repeated amino acids
        ll <- length(repaa)
        for (x in 1:ll) {
          ff <- sum(temp[temp$aa==repaa[x],]$freqf) #sum of equal amino acids frequencies
          temp$freqf[which(temp$aa==repaa[x])[1]] <- ff
          temp2 <- temp[-which(temp$aa==repaa[x])[-1],] 
          temp <- temp2
        }
        tt <- rbind(tt, temp)
      }
      
      tabpos <- unique(tt$pos)
      datapos <- unique(data3$pos)
      falt <- datapos[!datapos%in%tabpos] #positions that do not appear in table because they do not have repeated amino acids
      complem <- data3[data3$pos%in%falt,]
      tabela <- rbind(tt, complem)
      tabela <- tabela[order(tabela$pos, -tabela$freqf), ]
      
      #coverage
      covmat <- as.data.frame(matrix(pr1$cover, ncol = 3, byrow = T)) #coverage
      #taking the lowest coverage of each codon
      min_covmat <- apply(covmat,1,min)
      data_graph <- cbind.data.frame(seq(1:length(min_covmat)),min_covmat) 
      
    } 
    
    if (input$region == "Reverse Transcriptase"){
      rt1 <- tab2[pos >= 2550 & pos <= 3869, ]
      rt1[,4:7]/rt1[,3]
      rt3 <- cbind(rt1, rt1[,4:7]/rt1[,3])
      rt4 <- cbind(rt3[,1], rt3[,9:13], rt3[,10:13]) #position, reference, frequency
      colnames(rt4) <- c("pos", "ref", "freqa", "freqt", "freqg", "freqc", "a", "t", "g", "c")
      
      rt4$a[which(rt4$a!=0)] <- "A"
      rt4$t[which(rt4$t!=0)] <- "T"
      rt4$g[which(rt4$g!=0)] <- "G"
      rt4$c[which(rt4$c!=0)] <- "C"
      
      datafim <- data.frame() #dataframe to store the results
      
      i <- 1
      while (i <= (dim(rt4)[1])-2){ #length transcriptase reverse
        d <- rt4[i:(i+2),] #dataframe with one codon 
        aap <- (d$pos[3]-2549)/3 #amino acids positions 
        aar <- translate(s2c(paste(d$ref, collapse = ""))) #amino acids reference
        
        #assembling nucleotides
        dn <- data.frame(expand.grid(as.character(d[1,7:10]), as.character(d[2,7:10]), as.character(d[3,7:10]))) #combinations of codons
        df <- data.frame(expand.grid(as.numeric(d[1,3:6]), as.numeric(d[2,3:6]), as.numeric(d[3,3:6]))) #combinations of frequencies
        dd <- data.frame(cbind(dn, df))
        ddind <- apply(dd, 1, function(row) all(row !='0' )) #remove rows with zero
        dd2 <- dd[ddind,] #new dataframe with codons and frequencies
        colnames(dd2) <- c("pos1", "pos2", "pos3", "freq1", "freq2", "freq3") 
        
        if (dim(dd2)[1]!=0){
          for (ff in 1:dim(dd2)[1]){
            dd2$pos[ff] <- aap
            dd2$codonref[ff] <- paste(d$ref, collapse = "")
            dd2$aaref[ff] <- aar
            dd2$codon[ff] <- paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = "")
            dd2$aa[ff] <- translate(s2c(paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = ""))) #codification in amino acids
            dd2$freqf[ff] <- dd2$freq1[ff] * dd2$freq2[ff] * dd2$freq3[ff]
          }
        } else {
          ddd2 <- data.frame()
        }
        datafim <- rbind(datafim, dd2)
        i <- i+3
      }
      
      #Analyzing amino acids with frequency greater than or equal to 0.01%
      
      data <- datafim[,7:12] #limiting the number of columns
      data$freqt <- data$freqf*100 #frequency to percentage
      data2 <- data[data$freqt >=0.01,] #removing amino acids frequency <= 0.01%
      data3 <- data2
      data3$codon <- NULL; data3$codonref <- NULL #dataframe with only important columns 
      data3 <- subset(data3, !is.na(freqf)) #removing rows with NA
      
      tt <- data.frame()
      
      for (r in 1:as.numeric(data3$pos[length(data3$pos)])){
        temp <- data3[which(data3$pos==r),]
        repaa <- names(which((table(temp$aa)>=2)==TRUE)) #find repeated amino acids
        ll <- length(repaa)
        for (x in 1:ll) {
          ff <- sum(temp[temp$aa==repaa[x],]$freqf) #sum of equal amino acids frequencies
          temp$freqf[which(temp$aa==repaa[x])[1]] <- ff
          temp2 <- temp[-which(temp$aa==repaa[x])[-1],] 
          temp <- temp2
        }
        tt <- rbind(tt, temp)
      }
      
      tabpos <- unique(tt$pos)
      datapos <- unique(data3$pos)
      falt <- datapos[!datapos%in%tabpos] #positions that do not appear in table because they do not have repeated amino acids
      complem <- data3[data3$pos%in%falt,]
      tabela <- rbind(tt, complem)
      tabela <- tabela[order(tabela$pos, -tabela$freqf), ]
      
      #coverage
      covmat <- as.data.frame(matrix(rt1$cover, ncol = 3, byrow = T)) #coverage
      #taking the lowest coverage of each codon
      min_covmat <- apply(covmat,1,min)
      data_graph <- cbind.data.frame(seq(1:length(min_covmat)), min_covmat)
      
    } 
    
    if (input$region == "Integrase"){
      int1 <- tab2[pos >= 4230 & pos <= 5096, ] 
      int1[,4:7]/int1[,3]
      int3 <- cbind(int1, int1[,4:7]/int1[,3])
      int4 <- cbind(int3[,1], int3[,9:13], int3[,10:13]) #position, reference, frequency
      colnames(int4) <- c("pos", "ref", "freqa", "freqt", "freqg", "freqc", "a", "t", "g", "c")
      
      int4$a[which(int4$a!=0)] <- "A"
      int4$t[which(int4$t!=0)] <- "T"
      int4$g[which(int4$g!=0)] <- "G"
      int4$c[which(int4$c!=0)] <- "C"
      
      datafim <- data.frame() #dataframe to store the results
      
      i <- 1
      while (i <= (dim(int4)[1])-2){ #length integrase
        d <- int4[i:(i+2),] #dataframe with one codon 
        aap <- (d$pos[3]-4229)/3 #amino acids positions
        aar <- translate(s2c(paste(d$ref, collapse = ""))) #amino acids reference
        
        #assembling nucleotides
        dn <- data.frame(expand.grid(as.character(d[1,7:10]), as.character(d[2,7:10]), as.character(d[3,7:10]))) #combinations of codons
        df <- data.frame(expand.grid(as.numeric(d[1,3:6]), as.numeric(d[2,3:6]), as.numeric(d[3,3:6]))) #combinations of frequencies
        dd <- data.frame(cbind(dn, df))
        ddind <- apply(dd, 1, function(row) all(row !='0' )) #remove rows with zero
        dd2 <- dd[ddind,] #new dataframe with codons and frequencies
        colnames(dd2) <- c("pos1", "pos2", "pos3", "freq1", "freq2", "freq3") 
        
        if (dim(dd2)[1]!=0){
          for (ff in 1:dim(dd2)[1]){
            dd2$pos[ff] <- aap
            dd2$codonref[ff] <- paste(d$ref, collapse = "")
            dd2$aaref[ff] <- aar
            dd2$codon[ff] <- paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = "")
            dd2$aa[ff] <- translate(s2c(paste(dd2[ff,1], dd2[ff,2], dd2[ff,3], collapse = "", sep = ""))) #codification in amino acids
            dd2$freqf[ff] <- dd2$freq1[ff] * dd2$freq2[ff] * dd2$freq3[ff]
          }
        } else {
          ddd2 <- data.frame()
        }
        datafim <- rbind(datafim, dd2)
        i <- i+3
      }
      
      #Analyzing amino acids with frequency greater than or equal to 0.01%
      
      data <- datafim[,7:12] #limiting the number of columns
      data$freqt <- data$freqf*100 #frequency to percentage
      data2 <- data[data$freqt >=0.01,] #removing amino acids frequency <= 0.01%
      data3 <- data2
      data3$codon <- NULL; data3$codonref <- NULL #dataframe with only important columns
      data3 <- subset(data3, !is.na(freqf)) #removing rows with NA
      
      tt <- data.frame()
      
      for (r in 1:as.numeric(data3$pos[length(data3$pos)])){
        temp <- data3[which(data3$pos==r),]
        repaa <- names(which((table(temp$aa)>=2)==TRUE)) #find repeated amino acids
        ll <- length(repaa)
        for (x in 1:ll) {
          ff <- sum(temp[temp$aa==repaa[x],]$freqf) #sum of equal amino acids frequencies
          temp$freqf[which(temp$aa==repaa[x])[1]] <- ff
          temp2 <- temp[-which(temp$aa==repaa[x])[-1],] 
          temp <- temp2
        }
        tt <- rbind(tt, temp)
      }
      
      tabpos <- unique(tt$pos)
      datapos <- unique(data3$pos)
      falt <- datapos[!datapos%in%tabpos] #positions that do not appear in table because they do not have repeated amino acids
      complem <- data3[data3$pos%in%falt,]
      tabela <- rbind(tt, complem)
      tabela <- tabela[order(tabela$pos, -tabela$freqf), ]
      
      #coverage
      covmat <- as.data.frame(matrix(int1$cover, ncol = 3, byrow = T)) #coverage
      #taking the lowest coverage of each codon
      min_covmat <- apply(covmat,1,min)
      data_graph <- cbind.data.frame(seq(1:length(min_covmat)), min_covmat) 
      
    }
    
    data_graph
    
  })
  
  output$plot <- renderPlotly({
    if(input$Loadd == 0){return()} #just to not show warning messages
    
    d_graph <- datasetInput2()
    colnames(d_graph) <- c("Position","Coverage")
    plot_ly(d_graph, x = ~Position, y = ~Coverage, type = "bar") 
  })
  
  datahivdb <- eventReactive(input$load, {
    
    if(input$Load == 0){return()} #just to not show warning messages
    
    dposmut <- datasetInput()
    
    colnames(dposmut) <- c("cover", "ref", "pos", "aa", "freqf")
    
    if(input$region == "Protease"){
      
      #reading score files to PR
      
      Scores_PI <- read.table(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/HIVdb/Scores_PI_2018.txt"), header = T)
      
      drugName <- c("DRV.r","DRV.r","DRV.r","DRV.r","DRV.r","DRV.r","DRV.r","DRV.r","DRV.r","DRV.r","DRV.r",   
                    "LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r","LPV.r",   
                    "ATV.r","ATV.r","ATV.r","ATV.r","ATV.r","ATV.r","ATV.r","ATV.r","ATV.r","ATV.r","ATV.r",     
                    "SQV.r","SQV.r","SQV.r","SQV.r","SQV.r","SQV.r","SQV.r","SQV.r",   
                    "IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r","IDV.r",   
                    "FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r","FPV.r",    
                    "TPV.r",   
                    "NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV","NFV")
      
      p1 <- c(11,11,32,32,32,32,32,47,47,54,54,  
              11,11,32,32,32,32,32,46,46,46,47,47,54,54,54,54,82,   
              32,32,46,46,46,47,53,54,54,73,82,   
              46,46,46,53,54,54,73,82,
              11,11,32,32,32,32,32,46,46,46,46,47,47,53,54,54,54,54,73,82,   
              11,11,32,32,32,32,32,46,46,46,46,47,47,53,54,54,54,54,73,82,
              46,  
              11,11,32,32,32,32,32,46,46,46,46,47,47,53,54,54,54,54,73,82)
      
      a1 <- c("IL","IL","I","I","I","I","I","VA","VA","LM","LM",  
              "IL","IL","I","I","I","I","I","IL","ILV","ILV","AV","AV","ALMSTV","ALMSTV","LM","LM","ACFLMST",  
              "I","I","IL","ILV","ILV","AV","L","ALMSTV","ALMSTV","ACSTV","ACFLMST",   
              "IL","ILV","ILV","L","ALMSTV","ALMSTV","ACSTV","ACFLMST",    
              "IL","IL","I","I","I","I","I","IL","ILV","ILV","ILV","AV","AV","L","ALMSTV","ALMSTV","LM","LM","ACSTV","ACFLMST",  
              "IL","IL","I","I","I","I","I","IL","ILV","ILV","ILV","AV","AV","L","ALMSTV","ALMSTV","LM","LM","ACSTV","ACFLMST",   
              "IL",  
              "IL","IL","I","I","I","I","I","IL","ILV","ILV","ILV","AV","AV","L","ALMSTV","ALMSTV","LM","LM","ACSTV","ACFLMST")
      
      p2 <- c(32,54,47,54,76,84,89,54,84,84,89,  
              32,54,47,54,76,84,89,84,76,82,54,84,82,90,84,89,90,  
              47,54,84,82,90,54,90,82,90,90,90,    
              84,82,90,90,82,90,90,90,   
              32,54,47,54,76,84,89,84,76,82,90,54,84,90,82,90,84,89,90,90,   
              32,54,47,54,76,84,89,84,76,82,90,54,84,90,82,90,84,89,90,90,   
              84,  
              32,54,47,54,76,84,89,84,76,82,90,54,84,90,82,90,84,89,90,90)
      
      a2 <- c("I","ML","VA","ML","V","V","V","ML","V","V","V",   
              "I","ML","VA","ML","V","V","V","V","V","ACFLMST","LM","V","ACFLMST","M","V","V","M",  
              "AV","LM","V","ACFLMST","M","LM","M","ACFLMST","M","M","M",   
              "V","ACFLMST","M","M","ACFLMST","M","M","M",  
              "I","LM","AV","LM","V","V","V","V","V","ACFLMST","M","LM","V","M","ACFLMST","M","V","V","M","M",
              "I","LM","AV","LM","V","V","V","V","V","ACFLMST","M","LM","V","M","ACFLMST","M","V","V","M","M",    
              "V",   
              "I","LM","AV","LM","V","V","V","V","V","ACFLMST","M","LM","V","M","ACFLMST","M","V","V","M","M")
      
      p3 <- c(0,0,0,0,0,0,0,0,0,0,0,   
              0,0,0,0,0,0,0,90,0,0,0,0,0,0,0,0,0,   
              0,0,90,0,0,0,0,0,0,0,0,    
              90,0,0,0,0,0,0,0,  
              0,0,0,0,0,0,0,90,0,0,0,0,0,0,0,0,0,0,0,0,    
              0,0,0,0,0,0,0,90,0,0,0,0,0,0,0,0,0,0,0,0,   
              90,   
              0,0,0,0,0,0,0,90,0,0,0,0,0,0,0,0,0,0,0,0)
      
      a3 <- c("O","O","O","O","O","O","O","O","O","O","O",   
              "O","O","O","O","O","O","O","M","O","O","O","O","O","O","O","O","O",   
              "O","O","M","O","O","O","O","O","O","O","O",    
              "M","O","O","O","O","O","O","O",   
              "O","O","O","O","O","O","O","M","O","O","O","O","O","O","O","O","O","O","O","O",
              "O","O","O","O","O","O","O","M","O","O","O","O","O","O","O","O","O","O","O","O",    
              "M",   
              "O","O","O","O","O","O","O","M","O","O","O","O","O","O","O","O","O","O","O","O")
      
      score <- c(5,5,5,5,5,5,5,5,5,5,5,   
                 5,5,5,5,5,5,5,5,10,10,5,5,10,5,5,5,5,  
                 5,5,5,10,10,5,10,10,10,10,10,  
                 5,10,5,10,10,10,10,10,   
                 5,5,5,5,5,5,5,5,10,10,10,5,5,10,10,10,5,5,10,10,    
                 5,5,5,5,5,5,5,5,10,10,10,5,5,10,10,10,5,5,10,10,   
                 5,  
                 5,5,5,5,5,5,5,5,10,10,10,5,5,10,10,10,5,5,10,10)
      
      comb_PI <- cbind.data.frame(drugName, p1, a1, p2, a2, p3, a3, score)
      
      if (input$hivdb == ">=20%") {
        cutoff <- subset(dposmut, freqf >= 0.2)
      }
      
      if (input$hivdb == ">=1%") {
        cutoff <- subset(dposmut, freqf >= 0.01)
      }
      
      ##### HIVdb #####
      
      pa <- with(Scores_PI, paste0(Position, AA))
      pro <- cbind(pa, Scores_PI)
      aa <- cutoff$aa
      posi <- cutoff$pos
      mut <- paste0(posi, aa)
      new <- pro[pro$pa %in% mut,]
      res <- colSums(new[,4:dim(new)[2]]) #punctuation simple resistance
      z <- data.frame(names(res), res) #reorganizing the simple answer
      row.names(z) <- NULL
      names(z) <- c("drugName", "score")
      cc <- data.frame()
      tam <- length(new$pa) #items number found in table
      
      #punctuation combined resistance
      
      library(gtools)
      
      if (tam >= 3){
        ind <- permutations(tam, 2, seq(1:tam))
        
        for(k in 1:dim(ind)[1]){
          ind <- permutations(tam, 2, seq(1:tam))[k,] #to generate the combinations between the positions
          #combination
          a <- comb_PI[grep(new$Position[ind[1]], comb_PI[,2]),] #checks column by column in the data set with the combined stanford mutations; grep finds a certain pattern
          b <- a[grep(new$AA[ind[1]], a[,3]),] 
          c <- b[grep(new$Position[ind[2]], b[,4]),] 
          d <- c[grep(new$AA[ind[2]], c[,5]),]
          
          indNA <- apply(d, 1, function(x) all(is.na(x)))
          d <- d[!indNA, ]
          
          cc <- rbind(cc, d)
          cc <- unique(cc) #remove repeated lines
        }
        
        if (dim(cc)[1] == 0){ #in case there is no combination, the data set will be empty
          cfim <- z
          
        } else {
          cc1 <- cc[with(cc, order(drugName)), ] #orders the drugs to stay close to each other
          cc2 <- cc1[,c("drugName", "score")] 
          cc3 <- aggregate(. ~ drugName, data=cc2, FUN=sum) #sum the values of the same drugs
          
          cfim <- aggregate(. ~ drugName, data=rbind(cc3, z), FUN=sum)
        }
      }
      
      
      if (tam == 2){
        
        ind <- permutations(tam, 2, seq(1:tam))
        
        for(k in 1:dim(ind)[1]){
          ind <- permutations(tam, 2, seq(1:tam))[k,]
          a <- comb_PI[grep(new$Position[1], comb_PI[,2]),] #checks column by column
          b <- a[grep(new$AA[1], a[,3]),] 
          c <- b[grep(new$Position[2], b[,4]),] 
          d <- c[grep(new$AA[2], c[,5]),]
          
          indNA <- apply(d, 1, function(x) all(is.na(x)))
          d <- d[!indNA, ]
          
          cc <- d
          cc <- unique(cc)
          
        }
        
        if (dim(cc)[1] == 0){ #in case there is no combination, the data set will be empty
          cfim <- z
          
        } else {
          cc1 <- cc[with(cc, order(drugName)), ] #orders the drugs to stay close to each other
          cc2 <- cc1[,c("drugName", "score")] 
          cc3 <- aggregate(. ~ drugName, data=cc2, FUN=sum) #sum the values of the same drugs
          
          cfim <- aggregate(. ~ drugName, data=rbind(cc3, z), FUN=sum)
        }
      }
      
      if (tam <= 1){
        cfim <- z
      }
      
      cfim$class <- rep(NA, nrow(cfim)) #added a column in dataframe "cfim"
      
      #Stanford (HIVdb) classifications
      
      for (f in 1:dim(cfim)[1]){
        cfim$score[f]
        if (cfim$score[f] <= 9){cfim$class[f] <- "Susceptible"}
        if (cfim$score[f] >= 10 & cfim$score[f] <= 14){cfim$class[f] <- "Potential low-level resistance"}
        if (cfim$score[f] >= 15 & cfim$score[f] <= 29){cfim$class[f] <- "Low-level resistance"}
        if (cfim$score[f] >= 30 & cfim$score[f] <= 59){cfim$class[f] <- "Intermediate resistance"}
        if (cfim$score[f] >= 60){cfim$class[f] <- "High-level resistance"}
      }
      
      cfim$drugName <- gsub("ATV.r", "atazanavir/r (ATV/r)", cfim$drugName)
      cfim$drugName <- gsub("DRV.r", "darunavir/r (DRV/r)", cfim$drugName)
      cfim$drugName <- gsub("FPV.r", "fosamprenavir/r (FPV/r)", cfim$drugName)
      cfim$drugName <- gsub("IDV.r", "indinavir/r (IDV/r)", cfim$drugName)
      cfim$drugName <- gsub("LPV.r", "lopinavir/r (LPV/r)", cfim$drugName)
      cfim$drugName <- gsub("NFV", "nelfinavir (NFV)", cfim$drugName)
      cfim$drugName <- gsub("SQV.r", "saquinavir/r (SQV/r)", cfim$drugName)
      cfim$drugName <- gsub("TPV.r", "tipranavir/r (TPV/r)", cfim$drugName)
      
      names(cfim) <- c("Protease Inhibitors", "Mutation Scoring", "HIVdb")
      
      res_hivdb <- cfim
      
      ##### ANRS #####
      
      mut <- paste0(posi, aa)
      
      #ATV 
      
      atvclass <- NA
      
      atv1 <- c("50L","88S")
      atv1_len <- sum(atv1%in%mut)
      
      atv2 <- c("10F/I/V", "16E", "33F/I/V", "46I/L", "60E", "71V/T", "84V", "85V", "90M")
      atv2_out <- map(atv2, mut)
      atv2_len <- sum(sapply(atv2_out, any))
      
      if(atv1_len >= 1 | atv2_len >= 3){
        atvclass <- "Resistant"
      } else if (atv2_len == 2){
        atvclass <- "Possible Resistant"
      } else {
        atvclass <- "Susceptible"
      }
      
      #DRV
      
      drvclass <- NA
      
      drv1 <- c("11I", "32I", "33F", "47V", "50V", "54L/M", "74P", "76V", "84V", "89V")
      drv1_out <- map(drv1, mut)
      drv1_len <- sum(sapply(drv1_out, any))
      
      if(drv1_len >= 4 ){
        drvclass <- "Resistant"
      } else if (drv1_len == 3){
        drvclass <- "Possible Resistant"
      } else {
        drvclass <- "Susceptible"
      }
      
      #FPV
      
      fpvclass <- NA
      
      fpv1 <- c("50V")
      fpv1_len <- sum(fpv1%in%mut)
      
      fpv2 <- c("32I")
      fpv2_len <- sum(fpv2%in%mut)
      
      fpv3 <- c("47A/V")
      fpv3_out <- map(fpv3, mut)
      fpv3_len <- ifelse(sum(sapply(fpv3_out, any))>=1, 1,0) #when it has only one position but more than one possibility of aa
      
      fpv4 <- c("10F/I/V", "33F", "36I", "54A/L/M/S/T/V", "62V", 
                "82A/C/F/G", "84V", "90M")
      fpv4_out <- map(fpv4, mut)  
      fpv4_len <- sum(sapply(fpv4_out, any))
      
      if(fpv1_len == 1 | (fpv2_len == 1 & fpv3_len == 1) | fpv4_len >= 4){
        fpvclass <- "Resistant"
      } else {
        fpvclass <- "Susceptible"
      }
      
      #IDV
      
      idvclass <- NA
      
      idv1 <- c("46I/L", "82A/F/M/S/T", "84A", "84V")
      idv1_out <- map(idv1, mut)
      idv1_len <- sum(sapply(idv1_out, any))
      
      idv2 <- c("90M")
      idv2_len <- sum(idv2%in%mut)
      
      idv3 <- c("20M/R", "24I", "32I", "36I","54V/L/M/T", "71V/T", "73S/A", "77I")
      idv3_out <- map(idv3, mut)
      idv3_len <- sum(sapply(idv3_out, any))
      
      if(idv1_len >= 1 | (idv2_len == 1 & idv3_len >= 2)){
        idvclass <- "Resistant"
      } else if ((idv2_len == 1)){
        idvclass <- "Possible Resistant"
      } else {
        idvclass <- "Susceptible"
      }
      
      #LPV
      
      lpvclass <- NA
      
      lpv1 <- c("47A", "76V")
      lpv1_len <- sum(lpv1%in%mut)
      
      lpv2 <- c("10F/I/R/V", "20M/R", "24I", "33F", "46I/L", "50V",
                "53L","54M/L/T/V","63P","71I/L/T/V", "82A/F/S/T", "84V","90M")
      lpv2_out <- map(lpv2, mut)
      lpv2_len <- sum(sapply(lpv2_out, any))
      
      if(lpv1_len >= 1 | lpv2_len >= 4){
        lpvclass <- "Resistant"
      } else if (lpv2_len == 3){
        lpvclass <- "Possible Resistant"
      } else {
        lpvclass <- "Susceptible"
      }
      
      #SQV
      
      sqvclass <- NA
      
      sqv1 <- c("10F/I/M/R/V", "15A/V", "20I/M/R/T", "24I","62V",
                "73S/T", "82A/F/T/S", "84V", "90M")
      sqv1_out <- map(sqv1, mut)
      sqv1_len <- sum(sapply(sqv1_out, any))
      
      sqv2 <- c("48V")
      sqv2_len <- sum(sqv2%in%mut)
      
      if(sqv1_len >= 3 | sqv2_len == 1){
        sqvclass <- "Resistant"
      } else if (sqv1_len == 2){
        sqvclass <- "Possible Resistant"
      } else {
        sqvclass <- "Susceptible"
      }
      
      #TPV
      
      tpvclass <- NA
      
      tpv1 <- c("36I/L/V")
      tpv1_out <- map(tpv1, mut)
      tpv1_len <- ifelse(sum(sapply(tpv1_out, any))>=1,1,0)
      
      tpv2 <- c("53L/W/Y")
      tpv2_out <- map(tpv2, mut)  
      tpv2_len <- ifelse(sum(sapply(tpv2_out, any))>=1,1,0)
      
      tpv3 <- c("58E")
      tpv3_len <- sum(tpv3%in%mut)
      
      tpv4 <- c("69I/K/N/Q/R/Y")
      tpv4_out <- map(tpv4, mut)
      tpv4_len <- ifelse(sum(sapply(tpv4_out, any))>=1,1,0)
      
      tpv5 <- c("89I/M/R/T/V")
      tpv5_out <- map(tpv5, mut)
      tpv5_len <- ifelse(sum(sapply(tpv5_out, any))>=1,1,0)
      
      scoretpv <- tpv1_len - tpv2_len + tpv3_len + tpv4_len + tpv5_len
      
      if(scoretpv >= 3){
        tpvclass <- "Resistant"
      } else if (scoretpv == 2){
        tpvclass <- "Possible Resistant"
      } else {
        tpvclass <- "Susceptible"
      }
      
      #NFV
      
      nfvclass <- NA
      
      nfv1 <- c("30N", "84A/V", "88S/D", "90M")
      nfv1_out <- map(nfv1, mut)
      nfv1_len <- sum(sapply(nfv1_out, any))
      
      nfv2 <- c("82A/F/S/T")
      nfv2_out <- map(nfv2, mut)
      nfv2_len <- ifelse(sum(sapply(nfv2_out, any))>=1,1,0)
      
      nfv3 <- c("10I", "36I", "46I/L", "54V/L/M/T", "71V/T/I")
      nfv3_out <- map(nfv3, mut)
      nfv3_len <- sum(sapply(nfv3_out, any))
      
      if(nfv1_len >= 1){
        nfvclass <- "Resistant"
      } else if (nfv2_len == 1 & nfv3_len >= 2){
        nfvclass <- "Possible Resistant"
      } else {
        nfvclass <- "Susceptible"
      }
      
      res_anrs <- rbind.data.frame(c("atazanavir/r (ATV/r)", atvclass),
                                   c("darunavir/r (DRV/r)", drvclass),
                                   c("fosamprenavir/r (FPV/r)", fpvclass),
                                   c("indinavir/r (IDV/r)", idvclass),
                                   c("lopinavir/r (LPV/r)", lpvclass),
                                   c("nelfinavir (NFV)", nfvclass),
                                   c("saquinavir/r (SQV/r)", sqvclass),
                                   c("tipranavir/r (TPV/r)", tpvclass)
      )
      
      colnames(res_anrs) <- c("drug", "ANRS")
      
      ##### Rega #####
      
      mut <- paste0(posi, aa)
      
      #ATV
      
      atvclass <- NA
      
      atv1 <- c("48M/V","50L")
      atv1_out <- map(atv1, mut)
      atv1_len <- sum(sapply(atv1_out, any))*2
      
      atv2 <- c("54A/S/T", "88S")
      atv2_out <- map(atv2, mut)
      atv2_len <- sum(sapply(atv2_out, any))*1.5
      
      atv3 <- c("10Y","20T","47V","54V","74P","82L/T","84V","88D","90M")
      atv3_out <- map(atv3, mut)
      atv3_len <- sum(sapply(atv3_out, any))*1
      
      atv4 <- c("10F","24I","32I","46I/L","54L/M","58E","71I/L","73A/C/F/S/T","82A/M/S/F","85V")
      atv4_out <- map(atv4, mut)
      atv4_len <- sum(sapply(atv4_out, any))*0.5
      
      atv5 <- c("10I/V","20I/M/R/V","33F/I","71T/V","74A/S")
      atv5_out <- map(atv5, mut)
      atv5_len <- sum(sapply(atv5_out, any))*0.25
      
      atv6 <- c("76V")
      atv6_len <- sum(atv6%in%mut)
      
      scoreatv <- atv1_len + atv2_len + atv3_len + atv4_len + atv5_len + atv6_len
      
      if(scoreatv >= 3.5){
        atvclass <- "Resistant"
      } else if (scoreatv >= 2 & scoreatv < 3.5){
        atvclass <- "Intermediate Resistant"
      } else {
        atvclass <- "Susceptible"
      }
      
      #DRV
      
      drvclass <- NA
      
      drv1 <- c("50V/M","76V","84A/C/V")
      drv1_out <- map(drv1, mut)
      drv1_len <- sum(sapply(drv1_out, any))*1.5
      
      drv2 <- c("11I","32I","33F","41T","47A/V","54L","70E","74P","89V")
      drv2_out <- map(drv2, mut)
      drv2_len <- sum(sapply(drv2_out, any))*1
      
      drv3 <- c("32L","33M/V","34V","35G","41I","46I/L","48M","54S/T/V","73A/C/F/S/T/V","74E","82F/L","85V","89T")
      drv3_out <- map(drv3, mut)
      drv3_len <- sum(sapply(drv3_out, any))*0.5
      
      drv4 <- c("33I","35N","89I")
      drv4_len <- sum(drv4%in%mut)*0.25
      
      drv5 <- c("50L","88S")
      drv5_len <- sum(drv5%in%mut)*-0.25
      
      scoredrv <- drv1_len + drv2_len + drv3_len + drv4_len + drv5_len
      
      if(scoredrv >= 3.5){
        drvclass <- "Resistant"
      } else if (scoredrv >= 2 & scoredrv < 3.5){
        drvclass <- "Intermediate Resistant"
      } else {
        drvclass <- "Susceptible"
      }
      
      #FPV
      
      fpvclass <- NA
      
      fpv1 <- c("50V")
      fpv1_len <- sum(fpv1%in%mut)*2
      
      fpv2 <- c("76V","84A/C/V")
      fpv2_out <- map(fpv2, mut)
      fpv2_len <- sum(sapply(fpv2_out, any))*1.5
      
      fpv3 <- c("47A/V", "48M","54M/T/V", "82F")
      fpv3_out <- map(fpv3, mut)
      fpv3_len <- sum(sapply(fpv3_out, any))*1
      
      fpv4 <- c("10F","20T","24I","32I","33F","43T","46I/L", "48A","54A/L","58E","82A/M/S/T", 
                "89T/V")
      fpv4_out <- map(fpv4, mut)
      fpv4_len <- sum(sapply(fpv4_out, any))*0.5
      
      fpv5 <- c("10I/V", "20I/R/M/V", "33I","43R","48V","89I","90M")
      fpv5_out <- map(fpv5, mut)
      fpv5_len <- sum(sapply(fpv5_out, any))*0.25
      
      fpv6 <- c("50L","88S")
      fpv6_len <- sum(fpv6%in%mut)*-0.25
      
      scorefpv <- fpv1_len + fpv2_len + fpv3_len + fpv4_len + fpv5_len + fpv6_len
      
      if(scorefpv >= 3.5){
        fpvclass <- "Resistant"
      } else if (scorefpv >= 2 & scorefpv < 3.5){
        fpvclass <- "Intermediate Resistant"
      } else {
        fpvclass <- "Susceptible"
      }
      
      #IDV
      
      idvclass <- NA
      
      idv1 <- c("82A/F/S/T","84A/V")
      idv1_out <- map(idv1, mut)
      idv1_len <- sum(sapply(idv1_out, any))*2
      
      idv2 <- c("46I/L", "54A/S/T")
      idv2_out <- map(idv2, mut)
      idv2_len <- sum(sapply(idv2_out, any))*1.5
      
      idv3 <- c("10F","24I","32I","48M/V","54V","76V","82M","88S","90M")
      idv3_out <- map(idv3, mut)
      idv3_len <- sum(sapply(idv3_out, any))*1
      
      idv4 <- c("20T","43T","48A","54L/M","66F","71I", "73S/T/C","74P","84C","88D","89T/V","95F")
      idv4_out <- map(idv4, mut)
      idv4_len <- sum(sapply(idv4_out, any))*0.5
      
      idv5 <- c("10I/V","20I/M/R/V","35G/N","43R","62V","71T/V","74A/S","89I")
      idv5_out <- map(idv5, mut)
      idv5_len <- sum(sapply(idv5_out, any))*0.25
      
      idv6 <- c("50L")
      idv6_len <- sum(idv6%in%mut)*-0.25
      
      scoreidv <- idv1_len + idv2_len + idv3_len + idv4_len + idv5_len + idv6_len  
      
      if(scoreidv >= 3){
        idvclass <- "Resistant"
      } else if (scoreidv >= 2 & scoreidv < 3){
        idvclass <- "Intermediate Resistant"
      } else {
        idvclass <- "Susceptible"
      }
      
      #LPV
      
      lpvclass <- NA
      
      lpv1 <- c("47A")
      lpv1_len <- sum(lpv1%in%mut)*2
      
      lpv2 <- c("50V","54A/T","76V")
      lpv2_out <- map(lpv2, mut)
      lpv2_len <- sum(sapply(lpv2_out, any))*1.5
      
      lpv3 <- c("47V","48M","54S/V","82F","84S/A")
      lpv3_out <- map(lpv3, mut)
      lpv3_len <- sum(sapply(lpv3_out, any))*1
      
      lpv4 <- c("10F","20T","24F/I","32I","33F","43T","46I/L","48A/V", "53L","54L/M","71I",
                "73S/T/C","82A/L/M/T","84V","88D","90M")
      lpv4_out <- map(lpv4, mut)
      lpv4_len <- sum(sapply(lpv4_out, any))*0.5
      
      lpv5 <- c("10I/V","20I/M/R/V","33I","43R","64M/V","71T/V", "77A/T/V")
      lpv5_out <- map(lpv5, mut)
      lpv5_len <- sum(sapply(lpv5_out, any))*0.25
      
      lpv6 <- c("50L")
      lpv6_len <- sum(lpv6%in%mut)*-0.25
      
      scorelpv <- lpv1_len + lpv2_len + lpv3_len + lpv4_len + lpv5_len + lpv6_len
      
      if(scorelpv >= 3.5){
        lpvclass <- "Resistant"
      } else if (scorelpv >= 2 & scorelpv < 3.5){
        lpvclass <- "Intermediate Resistant"
      } else {
        lpvclass <- "Susceptible"
      }
      
      #SQV
      
      sqvclass <- NA
      
      sqv1 <- c("48M/V", "90M")
      sqv1_out <- map(sqv1, mut)
      sqv1_len <- sum(sapply(sqv1_out, any))*2
      
      sqv2 <- c("54A/S/T", "84A/C")
      sqv2_out <- map(sqv2, mut)
      sqv2_len <- sum(sapply(sqv2_out, any))*1.5
      
      sqv3 <- c("24I","53L","54V","71V","84V","88D/S", "89V")
      sqv3_out <- map(sqv3, mut)
      sqv3_len <- sum(sapply(sqv3_out, any))*1
      
      sqv4 <- c("10F","20T","46I/L", "48A","50V","54L/M", "58E","71I","73S/T/C", "74S", "75P","89T")
      sqv4_out <- map(sqv4, mut)
      sqv4_len <- sum(sapply(sqv4_out, any))*0.5
      
      sqv5 <- c("10I/V", "11I","20I/M/R/V", "62V","71T","74A","82A/F/L/M/S/T", "89I")
      sqv5_out <- map(sqv5, mut)
      sqv5_len <- sum(sapply(sqv5_out, any))*0.25
      
      sqv6 <- c("50L")
      sqv6_len <- sum(sqv6%in%mut)*-0.25
      
      scoresqv <- sqv1_len + sqv2_len + sqv3_len + sqv4_len + sqv5_len + sqv6_len 
      
      if(scoresqv >= 3.5){
        sqvclass <- "Resistant"
      } else if (scoresqv >= 2 & scoresqv < 3.5){
        sqvclass <- "Intermediate Resistant"
      } else {
        sqvclass <- "Susceptible"
      }
      
      #TPV
      
      tpvclass <- NA
      
      tpv1 <- c("47V", "74P","82T")
      tpv1_len <- sum(tpv1%in%mut)*2
      
      tpv2 <- c("58E","82L/S", "83D")
      tpv2_out <- map(tpv2, mut)
      tpv2_len <- sum(sapply(tpv2_out, any))*1.5
      
      tpv3 <- c("41T","43T","54A/M/V", "84A/C/V")
      tpv3_out <- map(tpv3, mut)
      tpv3_len <- sum(sapply(tpv3_out, any))*1
      
      tpv4 <- c("33F","47A","46L","54S/T")
      tpv4_out <- map(tpv4, mut)
      tpv4_len <- sum(sapply(tpv4_out, any))*0.5
      
      tpv5 <- c("11L","32I","38W","45I","71L","73T","82A/F/M", "89T/V", "90M")
      tpv5_out <- map(tpv5, mut)
      tpv5_len <- sum(sapply(tpv5_out, any))*0.25
      
      tpv6 <- c("50L")
      tpv6_len <- sum(tpv6%in%mut)*-0.25
      
      scoretpv <- tpv1_len + tpv2_len + tpv3_len + tpv4_len + tpv5_len + tpv6_len 
      
      if(scoretpv >= 3.5){
        tpvclass <- "Resistant"
      } else if (scoretpv >= 2 & scoretpv < 3.5){
        tpvclass <- "Intermediate Resistant"
      } else {
        tpvclass <- "Susceptible"
      }
      
      #NFV
      
      nfvclass <- NA
      
      nfv1 <- c("30N","90M")
      nfv1_len <- sum(nfv1%in%mut)*2
      
      nfv2 <- c("54A/S/T", "84A/C", "88S")
      nfv2_out <- map(nfv2, mut)
      nfv2_len <- sum(sapply(nfv2_out, any))*1.5
      
      nfv3 <- c("54V","82A/T/F","88D")
      nfv3_out <- map(nfv3, mut)
      nfv3_len <- sum(sapply(nfv3_out, any))*1
      
      nfv4 <- c("10F","20T","23I","24I","32I","43T","46I/L", "48A/V", "54L/M","58E","66F","71I",
                "73S/T/C", "74P","76V","82L/M/S", "84V","89T/V", "93M")
      nfv4_out <- map(nfv4, mut)
      nfv4_len <- sum(sapply(nfv4_out, any))*0.5
      
      nfv5 <- c("10I/V", "20I/M/R/V", "33F/I", "35G/N", "43R","62V","64V","71T/V","74A/S", "89I")
      nfv5_out <- map(nfv5, mut)
      nfv5_len <- sum(sapply(nfv5_out, any))*0.25
      
      nfv6 <- c("50L")
      nfv6_len <- sum(nfv6%in%mut)*-0.25
      
      scorenfv <- nfv1_len + nfv2_len + nfv3_len + nfv4_len + nfv5_len + nfv6_len 
      
      if(scorenfv >= 2){
        nfvclass <- "Resistant"
      } else if (scorenfv >= 1.25 & scorenfv < 2){
        nfvclass <- "Intermediate Resistant"
      } else {
        nfvclass <- "Susceptible"
      }
      
      res_rega <- rbind.data.frame(c("atazanavir/r (ATV/r)", atvclass),
                                   c("darunavir/r (DRV/r)", drvclass),
                                   c("fosamprenavir/r (FPV/r)", fpvclass),
                                   c("indinavir/r (IDV/r)", idvclass),
                                   c("lopinavir/r (LPV/r)", lpvclass),
                                   c("nelfinavir (NFV)", nfvclass),
                                   c("saquinavir/r (SQV/r)", sqvclass),
                                   c("tipranavir/r (TPV/r)", tpvclass))
      
      
      colnames(res_rega) <- c("drug", "Rega")
      
      ##### Brazilian Algorithm #####
      
      mut <- paste0(posi, aa)
      
      #ATV
      
      atvclass <- NA
      
      atv1 <- c("50L","88S")
      atv1_len <- sum(atv1%in%mut)
      
      atv2 <- c("10I/V/F/R", "16E", "20R/M/T/I", "24I", "32I", "33F/I/V", 
                "36I/L/V", "45R", "46I/L/V", "48M/V", "53L", "54L/M/V", 
                "60E", "63P", "71V/T/I", "73C/S/T/A", "82A/F/I/S/T/M/L/C", 
                "84A/C/V", "85V", "88D/T/S", "89M/V/Q/T", "90M")
      atv2_out <- map(atv2, mut)
      atv2_len <- sum(sapply(atv2_out, any))
      
      atv3 <- c("10F/I/V","16E","33F/I/V","46I/L","60E","84V","85V","90M")
      atv3_out <- map(atv3, mut)
      atv3_len <- sum(sapply(atv3_out, any))
      
      atv4 <- c("90M","84M/V/Q/T","82A/F/I/S/T/M/L/C","54L/M/V","88D")
      atv4_out <- map(atv4, mut)
      atv4_len <- sum(sapply(atv4_out, any))
      
      atv5 <- c("10I/V/F/R","16E","20R/M/T/I","24I","32I","33F/I/V","36I/L/V",
                "45R","46I/L/V","48M/V","53L","60E","63P","71V/T/I","73C/S/T/A",
                "85V","89M/V/Q/T")
      atv5_out <- map(atv5, mut)
      atv5_len <- sum(sapply(atv5_out, any))
      
      if(atv1_len >= 1 | atv2_len >= 8 | atv3_len >= 3 | (atv4_len >= 1 & atv5_len >= 3)){
        atvclass <- "Resistant"
      } else if (atv2_len == 6 | atv2_len == 7 | (atv4_len >= 1 & atv5_len >= 2)){
        atvclass <- "Intermediate Resistant"
      } else {
        atvclass <- "Susceptible"
      }
      
      #DRV
      
      drvclass <- NA
      
      drv1 <- c("10F","11L/I","15V", "32I", "33F","35G", "41I/T","46I/L",
                "47A/V", "48M", "50V", "54L/M/V/A", "70E",
                "73S/T/A/C", "74P","76V","84V","85V", "89V")
      drv1_out <- map(drv1, mut)
      drv1_len <- sum(sapply(drv1_out, any))
      
      drv2 <- c("11I", "32I", "33F", "47F/V", "50V", "54L/M", "73S", "74P",
                "76V", "84V", "89V")
      drv2_out <- map(drv2, mut)
      drv2_len <- sum(sapply(drv2_out, any))
      
      if(drv1_len >= 8 | drv2_len >= 5){
        drvclass <- "Resistant"
      } else if (drv1_len == 6 | drv1_len == 7 | drv2_len == 3 | drv2_len == 4){
        drvclass <- "Intermediate Resistant"
      } else {
        drvclass <- "Susceptible"
      }
      
      #FPV
      
      fpvclass <- NA
      
      fpv1 <- c("10F/I/V/R","11I","20M/R/T/I","24I","32I","33F","35D","41K",
                "43R","46I/L","47V","48M","50V","54L/V/M/A/T/S","58E","62V", 
                "63P","71V/T","73S/T/C/A","76V","82A/F/I/T/S/M/C","84V/A/C",
                "85V","89V/T","90M","93L")
      fpv1_out <- map(fpv1, mut)
      fpv1_len <- sum(sapply(fpv1_out, any))
      
      if(fpv1_len >= 8){
        fpvclass <- "Resistant"
      } else if (fpv1_len == 6 | fpv1_len ==7){
        fpvclass <- "Intermediate Resistant"
      } else {
        fpvclass <- "Susceptible"
      }
      
      #IDV
      
      idvclass <- NA
      
      idv1 <- c("10F/I/R/V", "20M/R/T/I", "24I", "32I" ,"36L/V", "46I/L",
                "48V", "54L/T/V", "58E", "71/T/V/I", "73A/S/T/C",
                "76V", "77I","82A/F/I/S/T", "84V/A/C", "88S/D", "90M","93L")
      idv1_out <- map(idv1, mut)
      idv1_len <- sum(sapply(idv1_out, any))
      
      if(idv1_len >= 5){
        idvclass <- "Resistant"
      } else if ((idv1_len == 4)){
        idvclass <- "Intermediate Resistant"
      } else {
        idvclass <- "Susceptible"
      }
      
      #LPV
      
      lpvclass <- NA
      
      lpv1 <- c("47A")
      lpv1_len <- sum(lpv1%in%mut)
      
      lpv2 <- c("10F/I/R/V", "16A/E", "20M/R/T/I", "24I/V", "32I",
                "33F","34Q", "36I/V", "43T", "46I/L", "48V/M", "50V",
                "53L","54A/M/L/S/T/V","55E/R", "58E","63P",
                "71I/L/T/V", "73S/P/T/C/A","74S", "76V",
                "82A/F/I/S/T/M/L/C", "84V/A/C","89I/M/V", "90M",
                "91S")
      lpv2_out <- map(lpv2, mut)
      lpv2_len <- sum(sapply(lpv2_out, any))
      
      lpv3 <- c("10F/I/R/V", "20M/R", "24I", "33F", "46I/L","50V",
                "53L","54M/L/T/V","63P", "71I/L/T/V",
                "76V","82A/F/S/T", "84V", "90M")
      lpv3_out <- map(lpv3, mut)
      lpv3_len <- sum(sapply(lpv3_out, any))
      
      lpv4 <- c("50V","54A/M/L/S/T/V","82A/F/I/S/T/M/L/C","84V/A/C","90M")
      lpv4_out <- map(lpv4, mut)
      lpv4_len <- sum(sapply(lpv4_out, any))
      
      lpv5 <- c("10F/I/R/V","20M/R/T/I","24I/V","33F","46I/L","47A/V",
                "48V/M","53L","71I/L/T/V","73S/P/T/C/A")
      lpv5_out <- map(lpv5, mut)
      lpv5_len <- sum(sapply(lpv5_out, any))
      
      if(lpv1_len == 1 | lpv3_len >= 6 | lpv4_len >= 3 | (lpv4_len >= 2 & lpv5_len >= 3)){
        lpvclass <- "Resistant"
      } else if (lpv2_len == 6 | lpv2_len == 7 | lpv4_len >= 2 | lpv3_len == 4 | lpv3_len == 5){
        lpvclass <- "Intermediate Resistant"
      } else {
        lpvclass <- "Susceptible"
      }
      
      #SQV
      
      sqvclass <- NA
      
      sqv1 <- c("10F/I/R/V", "20M/R/T/I", "24I","36V/L", "46I/L",
                "48V","53L","54L/T/V/M", "58E", "62V",
                "71V/T/I","73S/T/C/A", "74A", "82A/F/T/I/S",
                "84V/A/C", "88D/S", "90M")
      sqv1_out <- map(sqv1, mut)
      sqv1_len <- sum(sapply(sqv1_out, any))
      
      if(sqv1_len >= 5){
        sqvclass <- "Resistant"
      } else if (sqv1_len == 4){
        sqvclass <- "Intermediate Resistant"
      } else {
        sqvclass <- "Susceptible"
      }
      
      #TPV
      
      tpvclass <- NA
      
      tpv1 <- c("10V/F/I", "13V", "20M/R/V","24M", "32I", "33F", "35G",
                "36I/L/V", "37D", "43T", "45I", "46I/L", "47A/V", "48V",
                "54A/M/V/S/T", "58E", "69I/K/N/Q/R/Y","71V",
                "73C/S/T/A","74P", "77I", "82AL/T", "83D", "84A/V/C",
                "89V/I/M/R/T")
      tpv1_out <- map(tpv1, mut)
      tpv1_len <- sum(sapply(tpv1_out, any))
      
      if(tpv1_len >= 8){
        tpvclass <- "Resistant"
      } else if (tpv1_len == 6 | tpv1_len == 7){
        tpvclass <- "Intermediate Resistant"
      } else {
        tpvclass <- "Susceptible"
      }
      
      nfvclass <- "Not present"
      
      res_bra <- rbind.data.frame(c("atazanavir/r (ATV/r)", atvclass),
                                  c("darunavir/r (DRV/r)", drvclass),
                                  c("fosamprenavir/r (FPV/r)", fpvclass),
                                  c("indinavir/r (IDV/r)", idvclass),
                                  c("lopinavir/r (LPV/r)", lpvclass),
                                  c("nelfinavir (NFV)", nfvclass),
                                  c("saquinavir/r (SQV/r)", sqvclass),
                                  c("tipranavir/r (TPV/r)", tpvclass))
      
      colnames(res_bra) <- c("drug", "Brazilian Algorithm")
      
      fim <- cbind.data.frame(res_hivdb[1], res_hivdb[3], res_rega[2], res_anrs[2], res_bra[2])
      
      clasfim <- fim
      
    }
    
    
    if(input$region == "Reverse Transcriptase"){
      
      Scores_NNRTI <- read.table(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/HIVdb/Scores_NNRTI_2018.txt"), header = T)
      
      Scores_NRTI <- read.table(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/HIVdb/Scores_NRTI_2018.txt"), header = T)
      
      drugName <- c("DOR","DOR","DOR","DOR","DOR","DOR","DOR","DOR","DOR","DOR","DOR",
                    "EFV","EFV","EFV","EFV",
                    "ETR","ETR","ETR","ETR","ETR","ETR","ETR","ETR",
                    "NVP","NVP","NVP",   
                    "RPV","RPV","RPV","RPV","RPV","RPV","RPV")
      
      p1 <- c(100,101,101,103,106,108,108,181,181,98,98,
              101,103,106,98,  
              101,101,101,101,179,179,181,98,  
              101,103,98,  
              101,103,138,179,179,181,98)
      
      a1 <- c("I","E","E","N","A","I","I","C","C","G","G",
              "E","R","A","G",  
              "E","E","E","E","F","T","C","G",  
              "E","R","G",  
              "E","R","K","F","T","C","G")
      
      p2 <- c(103,190,190,181,227,181,234,190,190,181,227, 
              181,179,227,181,  
              181,188,190,190,181,181,190,181,  
              181,179,181,  
              184,179,184,181,181,190,181)
      
      a2 <- c("N","A","S","C","L","C","I","A","CSTV","C","C",
              "C","D","L","C",  
              "C","L","A","S","C","C","ACSTV","C",  
              "c","D","C",  
              "I","D","I","C","C","ACSTV","C")
      
      score <- c(15,5,5,10,15,10,15,20,10,5,15,
                 5,20,15,5,  
                 5,5,5,5,15,10,10,5,  
                 5,20,5,  
                 15,15,15,15,10,10,5)
      
      comb_NNRTI <- cbind.data.frame(drugName, p1, a1, p2, a2, score)
      
      drugName2 <- c("ABC","ABC","ABC","ABC","ABC","ABC","ABC","ABC","ABC","ABC","ABC","ABC","ABC","ABC","ABC",
                     "DDI","DDI","DDI","DDI","DDI","DDI","DDI","DDI","DDI","DDI","DDI","DDI",
                     "FTC","FTC","FTC","FTC","FTC","FTC",
                     "X3TC","X3TC","X3TC","X3TC","X3TC","X3TC",
                     "D4T","D4T","D4T","D4T","D4T","D4T","D4T","D4T","D4T","D4T","D4T","D4T","D4T","D4T",
                     "TDF","TDF","TDF","TDF","TDF","TDF","TDF","TDF","TDF","TDF","TDF","TDF","TDF","TDF","TDF","TDF",
                     "AZT","AZT","AZT","AZT","AZT","AZT","AZT","AZT","AZT","AZT","AZT","AZT","AZT","AZT")
      
      p12<- c(210,210,40,41,41,41,41,41,41,67,67,67,70,74,77,
              210,210,40,41,41,41,41,41,67,67,70,77,
              210,41,41,41,67,77,
              210,41,41,41,67,77,
              151,210,210,40,41,41,41,41,41,67,67,70,70,77,
              151,210,210,40,41,41,41,41,41,41,65,67,67,70,70,77,
              151,210,210,40,41,41,41,41,41,65,67,67,70,77)
      
      a12 <- c("W","W","F","L","L","L","L","L","L","EGN","EGN","EGN","R","IV","L", 
               "W","W","F","L","L","L","L","L","EGN","EGN","R","L",
               "W","L","L","L","EGN","L",
               "W","L","L","L","EGN","L",
               "M","W","W","F","L","L","L","L","L","EGN","EGN","EGNQST","R","L",
               "M","W","W","F","L","L","L","L","L","L","R","EGN","EGN","EGNQST","R","L",
               "M","W","W","F","L","L","L","L","L","R","EGN","EGN","R","L")
      
      p22 <- c(215,215,41,210,210,215,215,44,67,215,70,70,215,184,116,  
               215,215,41,210,215,215,44,67,215,70,215,116,
               215,210,215,67,70,116,
               215,210,215,67,70,116,
               184,215,215,41,210,215,215,44,67,215,70,184,215,116,
               184,215,215,41,210,210,215,215,44,67,151,215,70,184,215,116,
               184,215,215,41,210,215,215,44,67,151,215,70,215,116)
      
      a22 <- c("ACDEILNSV","FY","L","W","W","ACDEILNSV","FY","AD","EGN","FY","R","R","FY","IV","Y",  
               "ACDEILNSV","FY","L","W","ACDEILNSV","FY","AD","EGN","FY","R","FY","Y",
               "FY","W","FY","EGN","R","Y",
               "FY","W","FY","EGN","R","Y",
               "IV","ACDEILNSV","FY","L","W","ACDEILNSV","FY","AD","EGN","FY","R","IV","FY","Y",
               "IV","ACDEILNSV","FY","L","W","W","ACDEILNSV","FY","AD","EGN","M","FY","R","IV","FY","Y",
               "IV","ACDEILNSV","FY","L","W","ACDEILNSV","FY","AD","EGN","M","FY","R","FY","Y")
      
      p32 <- c(0,0,210,0,215,0,0,210,215,219,184,219,0,0,151, 
               0,0,210,0,0,0,210,215,219,219,0,151,
               0,215,0,215,219,151,
               0,215,0,215,219,151,
               0,0,0,210,0,0,0,210,215,219,219,0,0,151,
               0,0,0,210,0,215,0,0,210,215,0,219,219,0,0,151,
               0,0,0,210,0,0,0,210,215,0,219,219,0,151)
      
      a32 <- c("O","O","W","O","FY","O","O","W","FY","ENQR","IV","ENQR","O","O","M",
               "O","O","W","O","O","O","W","FY","ENQR","ENQR","O","M",
               "O","FY","O","FY","ENQR","M",
               "O","FY","O","FY","ENQR","M",
               "O","O","O","W","O","O","O","W","FY","ENQR","ENQR","O","O","M",
               "O","O","O","W","O","FY","O","O","W","FY","O","ENQR","ENQR","O","O","M",
               "O","O","O","W","O","O","O","W","FY","O","ENQR","ENQR","O","M")	
      
      p42 <- c(0,0,215,0,0,0,0,215,0,0,219,0,0,0,0, 
               0,0,215,0,0,0,215,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,215,0,0,0,215,0,0,0,0,0,0,
               0,0,0,215,0,0,0,0,215,0,0,0,0,0,0,0,
               0,0,0,215,0,0,0,215,0,0,0,0,0,0)
      
      a42 <- c("O","O","FY","O","O","O","O","FY","O","O","ENQR","O","O","O","O", 
               "O","O","FY","O","O","O","FY","O","O","O","O","O",
               "O","O","O","O","O","O",
               "O","O","O","O","O","O",
               "O","O","O","FY","O","O","O","FY","O","O","O","O","O","O",
               "O","O","O","FY","O","O","O","O","FY","O","O","O","O","O","O","O",
               "O","O","O","FY","O","O","O","FY","O","O","O","O","O","O")
      
      score2 <- c(5,10,5,10,5,5,15,5,5,5,20,10,5,15,10,
                  5,10,5,10,5,10,5,5,5,10,5,10,
                  5,5,5,5,10,15,
                  5,5,5,5,10,15,
                  10,5,10,5,10,5,10,5,5,5,10,10,5,10,
                  10,5,10,5,10,5,5,10,5,5,10,5,10,10,5,15,
                  10,5,10,5,10,5,10,5,5,10,5,10,5,10)
      
      comb_NRTI <- cbind.data.frame(drugName2, p12, a12, p22, a22, p32, a32, score2)
      
      if (input$hivdb == ">=20%"){
        
        cutoff <- subset(dposmut, freqf >= 0.2)
        
      }
      
      if (input$hivdb == ">=1%"){
        
        cutoff <- subset(dposmut, freqf >= 0.01)
      }
      
      ##### HIVdb #####
      
      ## NNRTI ##
      
      pa <- with(Scores_NNRTI, paste0(Position, AA))
      pro <- cbind(pa, Scores_NNRTI)
      aa <- cutoff$aa
      posi <- cutoff$pos
      mut <- paste0(posi, aa)
      new <- pro[pro$pa %in% mut,]
      res <- colSums(new[,4:dim(new)[2]]) #punctuation simple resistance
      z <- data.frame(names(res), res) #reorganizing the simple answer
      row.names(z) <- NULL
      names(z) <- c("drugName", "score")
      cc <- data.frame()
      tam <- length(new$pa) #items number found in table
      
      #punctuation combined resistance
      
      #positions should be in ascending order
      
      library(gtools)
      
      if (tam >= 3){
        ind <- permutations(tam, 2, seq(1:tam))
        
        for(k in 1:dim(ind)[1]){
          ind <- permutations(tam, 2, seq(1:tam))[k,] #to generate the combinations between the positions
          #combination
          a <- subset(comb_NNRTI, comb_NNRTI[,2] == new$Position[ind[1]])
          b <- a[grep(new$AA[ind[1]], a[,3]),] 
          c <- subset(b, b[,4] == new$Position[ind[2]])
          d <- c[grep(new$AA[ind[2]], c[,5]),]
          
          indNA <- apply(d, 1, function(x) all(is.na(x)))
          d <- d[!indNA, ]
          
          cc <- rbind(cc,d)
          cc <- unique(cc)
        }
        
        if (dim(cc)[1] == 0){ #in case there is no combination, the data set will be empty
          cfim <- z
          
        } else {
          cc1 <- cc[with(cc, order(drugName)), ] #orders the drugs to stay close to each other
          cc2 <- cc1[,c("drugName", "score")] 
          cc3 <- aggregate(. ~ drugName, data=cc2, FUN=sum) #sum the values of the same drugs
          
          cfim <- aggregate(. ~ drugName, data=rbind(cc3, z), FUN=sum)
        }
      }
      
      
      if (tam == 2){
        
        ind <- permutations(tam, 2, seq(1:tam))
        
        for(k in 1:dim(ind)[1]){
          ind <- permutations(tam, 2, seq(1:tam))[k,]
          a <- subset(comb_NNRTI, comb_NNRTI[,2] == new$Position[1])
          b <- a[grep(new$AA[1], a[,3]),] 
          c <- subset(b, b[,4] == new$Position[2])
          d <- c[grep(new$AA[2], c[,5]),]
          
          indNA <- apply(d, 1, function(x) all(is.na(x)))
          d <- d[!indNA, ]
          
          cc <- d
          cc <- unique(cc)
        }
        
        if (dim(cc)[1] == 0){ #in case there is no combination, the data set will be empty
          cfim <- z
          
        } else {
          cc1 <- cc[with(cc, order(drugName)), ] #orders the drugs to stay close to each other
          cc2 <- cc1[,c("drugName", "score")] 
          cc3 <- aggregate(. ~ drugName, data=cc2, FUN=sum) #sum the values of the same drugs
          
          cfim <- aggregate(. ~ drugName, data=rbind(cc2, z), FUN=sum)
        }
      }
      
      if (tam <= 1){
        cfim <- z
      }
      
      cfim$class <- rep(NA, nrow(cfim)) #added a column in dataframe "cfim"
      
      #Stanford (HIVdb) classifications
      
      for (f in 1:dim(cfim)[1]){
        cfim$score[f]
        if (cfim$score[f] <= 9){cfim$class[f] <- "Susceptible"}
        if (cfim$score[f] >= 10 & cfim$score[f] <= 14){cfim$class[f] <- "Potential low-level resistance"}
        if (cfim$score[f] >= 15 & cfim$score[f] <= 29){cfim$class[f] <- "Low-level resistance"}
        if (cfim$score[f] >= 30 & cfim$score[f] <= 59){cfim$class[f] <- "Intermediate resistance"}
        if (cfim$score[f] >= 60){cfim$class[f] <- "High-level resistance"}
      }
      
      cfim$drugName <- gsub("DOR", "doravirine (DOR)", cfim$drugName)
      cfim$drugName <- gsub("EFV", "efavirenz (EFV)", cfim$drugName)
      cfim$drugName <- gsub("ETR", "etravirine (ETR)", cfim$drugName)
      cfim$drugName <- gsub("NVP", "nevirapine (NVP)", cfim$drugName)
      cfim$drugName <- gsub("RPV", "rilpivirine (RPV)", cfim$drugName)
      
      names(cfim) <- c("Reverse Transcriptase Inhibitors", "Mutation Scoring", "HIVdb")
      
      nnrtifim <- cfim
      
      rm(cfim)
      
      ## NRTI ##
      
      pa <- with(Scores_NRTI, paste0(Position, AA))
      pro <- cbind(pa, Scores_NRTI)
      aa <- cutoff$aa
      posi <- cutoff$pos
      mut <- paste0(posi, aa)
      new <- pro[pro$pa %in% mut,]
      res <- colSums(new[,4:dim(new)[2]]) #punctuation simple resistance
      z <- data.frame(names(res), res) #reorganizing the simple answer
      row.names(z) <- NULL
      names(z) <- c("drugName2", "score2")
      cc <- data.frame()
      tam <- length(new$pa) #items number found in table
      
      #punctuation combined resistance
      
      #### positions should be in ascending order
      
      library(gtools)
      
      if (tam >= 3){ #in this case we can have 3 combinations, so the combinations of 3 in 3 only occur when it is greater than or equal to 4
        ind <- permutations(tam, 3, seq(1:tam)) #number of combinations
        
        for(k in 1:dim(ind)[1]){
          ind <- permutations(tam, 3, seq(1:tam))[k,] #to generate the combinations between the positions
          #combination
          a <- subset(comb_NRTI, comb_NRTI[,2] == new$Position[ind[1]])
          b <- a[grep(new$AA[ind[1]], a[,3]),] 
          c <- subset(b, b[,4] == new$Position[ind[2]])
          d <- c[grep(new$AA[ind[2]], c[,5]),]
          e <- subset(d, d[,6] == new$Position[ind[3]] | d[,6] == 0)
          if(dim(e)[1] != 0 & d[1,6] == 0){
            f <- e[grepl("O", e[,7]),]
          } else {
            f <- e[grepl(new$AA[ind[3]], e[,7]),]
          }
          
          indNA <- apply(f, 1, function(x) all(is.na(x)))
          f <- f[!indNA, ]
          
          cc <- rbind(cc,f)
          cc <- unique(cc) #remove repeated lines
          
        }
        
        if (dim(cc)[1] == 0){ #in case there is no combination, the data set will be empty
          cfim <- z
          
        } else {
          cc1 <- cc[with(cc, order(drugName2)), ] #orders the drugs to stay close to each other
          cc2 <- cc1[,c("drugName2", "score2")] 
          cc3 <- aggregate(. ~ drugName2, data=cc2, FUN=sum) #sum the values of the same drugs
          
          cfim <- aggregate(. ~ drugName2, data=rbind(cc3, z), FUN=sum)
        }
      }
      
      
      if (tam == 2){ #in this case we can have 3 combinations, so the combinations of 3 in 3 only occur when it is greater than or equal to 4
        ind <- permutations(tam, 2, seq(1:tam)) #number of combinations
        
        for(k in 1:dim(ind)[1]){
          ind <- permutations(tam, 2, seq(1:tam))[k,] #to generate the combinations between the positions
          #combination
          a <- subset(comb_NRTI, comb_NRTI[,2] == new$Position[ind[1]])
          b <- a[grep(new$AA[ind[1]], a[,3]),] 
          c <- subset(b, b[,4] == new$Position[ind[2]])
          d <- c[grep(new$AA[ind[2]], c[,5]),]
          e <- subset(d, d[,6] == new$Position[ind[3]] | d[,6] == 0)
          if(dim(e)[1] != 0 & d[1,6] == 0){
            f <- e[grepl("O", e[,7]),]
          } else {
            f <- e[grepl(new$AA[ind[3]], e[,7]),]
          }
          
          indNA <- apply(f, 1, function(x) all(is.na(x)))
          f <- f[!indNA, ]
          
          cc <- rbind(cc,f)
          cc <- unique(cc) #remove repeated lines
          
        }
        
        if (dim(cc)[1] == 0){ #in case there is no combination, the data set will be empty
          cfim <- z
          
        } else {
          cc1 <- cc[with(cc, order(drugName2)), ] #orders the drugs to stay close to each other
          cc2 <- cc1[,c("drugName2", "score2")] 
          cc3 <- aggregate(. ~ drugName2, data=cc2, FUN=sum) #sum the values of the same drugs
          
          cfim <- aggregate(. ~ drugName2, data=rbind(cc3, z), FUN=sum)
        }
      }
      
      if (tam <= 1){
        cfim <- z
      }
      
      cfim$class <- rep(NA, nrow(cfim)) #added a column in dataframe "cfim"
      #Stanford (HIVdb) classifications
      for (f in 1:dim(cfim)[1]){
        cfim$score[f]
        if (cfim$score[f] <= 9){cfim$class[f] <- "Susceptible"}
        if (cfim$score[f] >= 10 & cfim$score[f] <= 14){cfim$class[f] <- "Potential low-level resistance"}
        if (cfim$score[f] >= 15 & cfim$score[f] <= 29){cfim$class[f] <- "Low-level resistance"}
        if (cfim$score[f] >= 30 & cfim$score[f] <= 59){cfim$class[f] <- "Intermediate resistance"}
        if (cfim$score[f] >= 60){cfim$class[f] <- "High-level resistance"}
      }
      
      cfim$drugName2 <- gsub("ABC", "abacavir (ABC)", cfim$drugName2)
      cfim$drugName2 <- gsub("DDI", "didanosine (DDI)", cfim$drugName2)
      cfim$drugName2 <- gsub("FTC", "emtricitabine (FTC)", cfim$drugName2)
      cfim$drugName2 <- gsub("X3TC", "lamivudine (3TC)", cfim$drugName2)
      cfim$drugName2 <- gsub("D4T", "stavudine (D4T)", cfim$drugName2)
      cfim$drugName2 <- gsub("TDF", "tenofovir (TDF)", cfim$drugName2)
      cfim$drugName2 <- gsub("AZT", "zidovudine (AZT)", cfim$drugName2)
      
      names(cfim) <- c("Reverse Transcriptase Inhibitors", "Mutation Scoring", "HIVdb")
      
      nrtifim <- cfim
      
      totfim <- rbind.data.frame(nnrtifim, nrtifim)
      
      ##### ANRS #####
      
      mut <- paste0(posi, aa)
      
      #3TC = FTC
      
      tc3class <- NA
      
      tc31 <- c("69i","184V/I", "65R")
      tc31_out <- map(tc31, mut)
      tc31_len <- sum(sapply(tc31_out, any))
      
      tc32 <- c("151M")
      tc32_len <- sum(tc32%in%mut)
      
      if(tc31_len >= 1){
        tc3class <- "Resistant"
      } else if (tc32_len == 1){
        tc3class <- "Possible Resistant"
      } else {
        tc3class <- "Susceptible"
      }
      
      #ABC
      
      abcclass <- NA
      
      abc1 <- c("65R", "69i","74V/I", "151M", "115F")
      abc1_out <- map(abc1, mut)
      abc1_len <- sum(sapply(abc1_out, any))
      
      abc2 <- c("41L", "67N", "184V/I", "210W", "215A/C/D/G/E/H/I/L/N/S/V/Y/F")
      abc2_out <- map(abc2, mut)
      abc2_len <- sum(sapply(abc2_out, any))
      
      abc3 <- c("41L", "67N", "210W", "215A/C/D/G/E/H/I/L/N/S/V/Y/F")
      abc3_out <- map(abc3, mut)
      abc3_len <- sum(sapply(abc3_out, any))
      
      abc4 <- c("184V/I")
      abc4_out <- map(abc4, mut)
      abc4_len <- ifelse(sum(sapply(abc4_out, any))>=1, 1,0)
      
      if(abc1_len >= 1 | abc2_len >= 3){
        abcclass <- "Resistant"
      } else if (abc3_len == 2 | abc4_len == 1){
        abcclass <- "Possible Resistant"
      } else {
        abcclass <- "Susceptible"
      }
      
      #AZT
      
      aztclass <- NA
      
      azt1 <- c("69i", "215A/C/D/G/E/H/I/L/N/S/V/Y/F", "151M")
      azt1_out <- map(azt1, mut)
      azt1_len <- sum(sapply(azt1_out, any))
      
      azt2 <- c("41L", "67N", "70R", "210W", "219Q/E")
      azt2_out <- map(azt2, mut)
      azt2_len <- sum(sapply(azt2_out, any))
      
      if(azt1_len >= 1 | azt2_len >= 3){
        aztclass <- "Resistant"
      } else {
        aztclass <- "Susceptible"
      }
      
      #ddI
      
      ddiclass <- NA
      
      ddi1 <- c("65R", "74V/I", "151M", "69i")
      ddi1_out <- map(ddi1, mut)
      ddi1_len <- sum(sapply(ddi1_out, any))
      
      ddi2 <- c("184V/I")
      ddi2_out <- map(ddi2, mut)
      ddi2_len <- ifelse(sum(sapply(ddi2_out, any))>=1, 1,0)
      
      ddi3 <- c("41L")
      ddi3_len <- sum(ddi3%in%mut)
      
      ddi4 <- c("69D")
      ddi4_len <- sum(ddi4%in%mut)
      
      ddi5 <- c("215Y/F")
      ddi5_out <- map(ddi5, mut)
      ddi5_len <- ifelse(sum(sapply(ddi5_out, any))>=1, 1,0)
      
      ddi6 <- c("219Q/E")
      ddi6_out <- map(ddi6, mut)
      ddi6_len <- ifelse(sum(sapply(ddi6_out, any))>=1, 1,0)
      
      ddi7 <- c("70R")
      ddi7_len <- sum(ddi7%in%mut)
      
      scoreddi <- ddi3_len + ddi4_len + ddi5_len + ddi6_len - ddi7_len - ddi2_len
      
      if(ddi1_len >= 1 | scoreddi >= 2){
        ddiclass <- "Resistant"
      } else {
        ddiclass <- "Susceptible"
      }
      
      #TDF
      
      tdfclass <- NA
      
      tdf1 <- c("65R/E/N", "69i", "70E")
      tdf1_out <- map(tdf1, mut)
      tdf1_len <- sum(sapply(tdf1_out, any))
      
      tdf2 <- c("41L", "44D", "67N", "69D/N/S", "74V/I", "210W", "215A/C/D/G/E/H/I/L/N/S/V/Y/F")
      tdf2_out <- map(tdf2, mut)
      tdf2_len <- sum(sapply(tdf2_out, any))
      
      if(tdf1_len >= 1 | tdf2_len >= 4){
        tdfclass <- "Resistant"
      } else if (tdf2_len == 3){
        tdfclass <- "Possible Resistant"
      } else {
        tdfclass <- "Susceptible"
      }
      
      #d4T
      
      d4tclass <- NA
      
      d4t1 <- c("75A/M/S/T", "215A/C/D/G/E/H/I/L/N/S/V/Y/F", "K65R", "151M", "69i")
      d4t1_out <- map(d4t1, mut)
      d4t1_len <- sum(sapply(d4t1_out, any))
      
      d4t2 <- c("41L", "67N", "70R", "210W", "219Q/E")
      d4t2_out <- map(d4t2, mut)
      d4t2_len <- sum(sapply(d4t2_out, any))
      
      if(d4t1_len >= 1 | d4t2_len >= 3){
        d4tclass <- "Resistant"
      } else {
        d4tclass <- "Susceptible"
      }
      
      #EFV
      
      efvclass <- NA
      
      efv1 <- c("100I", "101E", "103N/S/T/H", "106M", "138K", "181C/I","188C/L","190A/S/E/Q/C/T/V","225H", "230L")
      efv1_out <- map(efv1, mut)
      efv1_len <- sum(sapply(efv1_out, any))
      
      if(efv1_len >= 1){
        efvclass <- "Resistant"
      } else {
        efvclass <- "Susceptible"
      }
      
      #NVP
      
      nvpclass <- NA
      
      nvp1 <- c("100I", "101E", "103N/S/T/H", "106A/M", "181C/I", "188C/H/L", "190A/C/S/E/Q/T/V", "230L")
      nvp1_out <- map(nvp1, mut)
      nvp1_len <- sum(sapply(nvp1_out, any))
      
      nvp2 <- c("138K")
      nvp2_len <- sum(nvp2%in%mut)
      
      if(nvp1_len >= 1){
        nvpclass <- "Resistant"
      } else if (nvp2_len == 1){
        nvpclass <- "Possible Resistant"
      } else {
        nvpclass <- "Susceptible"
      }
      
      #ETR
      
      etrclass <- NA
      
      etr1 <- c("138K", "181C/I/V", "221Y")
      etr1_out <- map(etr1, mut)
      etr1_len <- sum(sapply(etr1_out, any))
      
      etr2 <- c("90I","98G","100I","101E/H/I/P/R", "106I", "179D/F/I/L/M/T", "190A/S", "230L")
      etr2_out <- map(etr2, mut)
      etr2_len <- sum(sapply(etr2_out, any))
      
      etr3 <- c("138A/G/Q/R/S")
      etr3_out <- map(etr3, mut)
      etr3_len <- ifelse(sum(sapply(etr3_out, any))>=1,1,0)
      
      if(etr1_len >= 1 | etr2_len >= 3){
        etrclass <- "Resistant"
      } else if (etr2_len == 2 | etr3_len == 1){
        etrclass <- "Possible Resistant"
      } else {
        etrclass <- "Susceptible"
      }
      
      #RPV
      
      rpvclass <- NA
      
      rpv1 <- c("101E/P", "138A/G/K/Q/R/S", "179L", "181C/I/V", "188L", "221Y", "227C","230I/L/V")
      rpv1_out <- map(rpv1, mut)
      rpv1_len <- sum(sapply(rpv1_out, any))
      
      rpv2 <- c("100I")
      rpv2_len <- sum(rpv2%in%mut)
      
      rpv3 <- c("103N/S")
      rpv3_out <- map(rpv3, mut)
      rpv3_len <- ifelse(sum(sapply(rpv3_out, any))>=1,1,0)
      
      rpv4 <- c("103R")
      rpv4_len <- sum(rpv4%in%mut)
      
      rpv5 <- c("179D")
      rpv5_len <- sum(rpv5%in%mut)
      
      scorerpv1 <- rpv2_len + rpv3_len
      scorerpv2 <- rpv2_len + rpv4_len + rpv5_len
      
      if(rpv1_len >= 1 | scorerpv1 == 2 | scorerpv2 == 3){
        rpvclass <- "Resistant"
      } else if (rpv5_len == 1){
        rpvclass <- "Possible Resistant"
      } else {
        rpvclass <- "Susceptible"
      }
      
      #DOR
      
      dorclass <- NA
      
      dor1 <- c("106A/M", "188L", "190S", "230L", "227C")
      dor1_out <- map(dor1, mut)
      dor1_len <- sum(sapply(dor1_out, any))
      
      dor2 <- c("100I")
      dor2_len <- sum(dor2%in%mut)
      
      dor3 <- c("103N")
      dor3_len <- sum(dor3%in%mut)
      
      dor4 <- c("181C")
      dor4_len <- sum(dor4%in%mut)
      
      dor5 <- c("225H")
      dor5_len <- sum(dor5%in%mut)
      
      scoredor1 <- dor2_len + dor3_len
      scoredor2 <- dor3_len + dor4_len 
      scoredor3 <- dor3_len + dor5_len 
      
      if(dor1_len >= 1 | scoredor1 == 2 | scoredor2 == 2 | scoredor3 == 2){
        dorclass <- "Resistant"
      } else {
        dorclass <- "Susceptible"
      }
      
      res_anrs <- rbind.data.frame(c("doravirine (DOR)", dorclass),
                                   c("efavirenz (EFV)", efvclass),
                                   c("etravirine (ETR)", etrclass),
                                   c("nevirapine (NVP)", nvpclass),
                                   c("rilpivirine (RPV)", rpvclass),
                                   c("abacavir (ABC)", abcclass),
                                   c("didanosine (DDI)", ddiclass),
                                   c("emtricitabine (FTC)", tc3class),
                                   c("lamivudine (3TC)", tc3class),
                                   c("stavudine (D4T)", d4tclass),
                                   c("tenofovir (TDF)", tdfclass),
                                   c("zidovudine (AZT)", aztclass))
      
      colnames(res_anrs) <- c("drug", "ANRS")
      
      ##### Rega #####
      
      mut <- paste0(posi, aa)
      
      #3TC
      
      tc3class <- NA
      
      tc31 <- c("184V/I")
      tc31_out <- map(tc31, mut)
      tc31_len <- ifelse(sum(sapply(tc31_out, any))>=1,1,0)
      
      tc32 <- c("151M")
      tc32_len <- sum(tc32%in%mut)
      
      tc33 <- c("65N/R")
      tc33_out <- map(tc31, mut)
      tc33_len <- ifelse(sum(sapply(tc33_out, any))>=1,1,0)
      
      tc34 <- c("65N/R","67d","69i","70E/G")
      tc34_out <- map(tc34, mut)
      tc34_len <- sum(sapply(tc34_out, any))
      
      tc35 <- c("44A/D","118I)")
      tc35_out <- map(tc35, mut)
      tc35_len <- sum(sapply(tc35_out, any))
      
      tc36 <- c("41L","67N","69A/N","70R","210W","215A/C/D/E/G/H/I/L/N/S/V/F/Y","219E/H/N/Q/R")
      tc36_out <- map(tc36, mut)
      tc36_len <- sum(sapply(tc36_out, any))
      
      tc37 <- c("75I","77L","116Y")
      tc37_len <- sum(tc37%in%mut)
      
      if(tc31_len == 1 | (tc32_len == 1 & tc33_len == 1)){
        tc3class <- "Resistant"
      } else if (tc34_len >= 1 | (tc35_len >= 1 & tc36_len >= 3) | (tc32_len == 1 & tc37_len >= 3)){
        tc3class <- "Intermediate Resistant"
      } else {
        tc3class <- "Susceptible"
      }
      
      #ABC
      
      abcclass <- NA
      
      abc1 <- c("67d", "69G")
      abc1_len <- sum(abc1%in%mut)
      
      abc2 <- c("65N/R","69i","74I/V","115F","151M","184I/V")
      abc2_out <- map(abc2, mut)
      abc2_len <- sum(sapply(abc2_out, any))
      
      abc3 <- c("69i")
      abc3_len <- sum(abc3%in%mut)
      
      abc4 <- c("41L","67N","70R","210W","215A/C/D/E/G/H/I/L/N/S/V/F/Y","219E/H/N/Q/R")
      abc4_out <- map(abc4, mut)
      abc4_len <- sum(sapply(abc4_out, any))
      
      abc5 <- c("65N/R","74I/V","115F","184I/V")
      abc5_out <- map(abc5, mut)
      abc5_len <- sum(sapply(abc5_out, any))
      
      abc6 <- c("151M")
      abc6_len <- sum(abc6%in%mut)
      
      abc7 <- c("75I","77L","116Y")
      abc7_len <- sum(abc7%in%mut)
      
      abc8 <- c("69i","70E/G","151M")
      abc8_out <- map(abc8, mut)
      abc8_len <- sum(sapply(abc8_out, any))
      
      if(abc1_len >= 1 | abc2_len >= 2 | (abc3_len == 1 & abc4_len >= 3) | (abc5_len >= 1 & abc4_len >= 4) | (abc6_len == 1 & abc7_len >= 2)){
        abcclass <- "Resistant"
      } else if (abc8_len == 1 | (abc5_len == 1 & (abc4_len >= 2 & abc4_len <= 3)) | abc4_len >= 3){
        abcclass <- "Intermediate Resistant"
      } else {
        abcclass <- "Susceptible"
      }
      
      #AZT
      
      aztclass <- NA
      
      azt1 <- c("151M", "69i")
      azt1_len <- sum(azt1%in%mut)
      
      azt2 <- c("41L","67G/N","69A/N","70R","210W", "215A/C/D/E/G/H/I/L/N/S/V/F/Y","219E/H/N/Q/R")
      azt2_out <- map(azt2, mut)
      azt2_len <- sum(sapply(azt2_out, any))
      
      azt3 <- c("41L","210W","215Y")
      azt3_len <- sum(azt3%in%mut)
      
      azt4 <- c("184I/V")
      azt4_out <- map(azt4, mut)
      azt4_len <- ifelse(sum(sapply(azt4_out, any))>=1,1,0)
      
      azt5 <- c("74V")
      azt5_len <- sum(azt5%in%mut)
      
      azt6 <- c("67G/N","70R","215F","219E/Q")
      azt6_out <- map(azt6, mut)
      azt6_len <- sum(sapply(azt6_out, any))
      
      azt7 <- c("215A/C/D/E/G/H/I/L/N/S/V/F/Y")
      azt7_out <- map(azt7, mut)
      azt7_len <- ifelse(sum(sapply(azt7_out, any))>=1,1,0)
      
      if(azt1_len >= 1 | azt2_len >= 4 | (azt3_len >= 3 & azt4_len == 0 & azt5_len == 0) | azt6_len >= 3){
        aztclass <- "Resistant"
      } else if ((azt2_len >= 2 & azt2_len <= 3) | (azt7_len == 1 & azt4_len == 0)){
        aztclass <- "Intermediate Resistant"
      } else {
        aztclass <- "Susceptible"
      }
      
      #ddI
      
      ddiclass <- NA
      
      ddi1 <- c("69D/G/N", "69i", "151M")
      ddi1_out <- map(ddi1, mut)
      ddi1_len <- sum(sapply(ddi1_out, any))
      
      ddi2 <- c("184V/I")
      ddi2_out <- map(ddi2, mut)
      ddi2_len <- ifelse(sum(sapply(ddi2_out, any))>=1,1,0)
      
      ddi3 <- c("65R/N","74I/V")
      ddi3_out <- map(ddi3, mut)
      ddi3_len <- sum(sapply(ddi3_out, any))
      
      ddi4 <- c("41L","67N","70R","74I/V","210W","215A/C/D/E/G/H/I/L/N/S/V/F/Y","219E/H/N/Q/R")
      ddi4_out <- map(ddi4, mut)
      ddi4_len <- sum(sapply(ddi4_out, any))
      
      ddi5 <- c("65N/R","70E/G","74I/V","75T")
      ddi5_out <- map(ddi5, mut)
      ddi5_len <- sum(sapply(ddi5_out, any))
      
      ddi6 <- c("41L","215F/Y")
      ddi6_out <- map(ddi6, mut)
      ddi6_len <- sum(sapply(ddi6_out, any))
      
      ddi7 <- c("41L","67N","70R","210W","215A/C/D/E/G/H/I/L/N/S/V/F/Y","219E/H/N/Q/R")
      ddi7_out <- map(ddi7, mut)
      ddi7_len <- sum(sapply(ddi7_out, any))
      
      if(ddi1_len >= 1 | (ddi2_len == 1 & ddi3_len >= 1) | ddi4_len >= 5){
        ddiclass <- "Resistant"
      } else if (ddi5_len >= 1 | ddi6_len == 2 | ddi7_len >= 3){
        ddiclass <- "Intermediate Resistant"
      } else {
        ddiclass <- "Susceptible"
      }
      
      #TDF
      
      tdfclass <- NA
      
      tdf1 <- c("69i")
      tdf1_len <- sum(tdf1%in%mut)
      
      tdf2 <- c("65N/R")
      tdf2_out <- map(tdf2, mut)
      tdf2_len <- ifelse(sum(sapply(tdf2_out, any))>=1,1,0)
      
      tdf3 <- c("184I/V")
      tdf3_out <- map(tdf3, mut)
      tdf3_len <- ifelse(sum(sapply(tdf3_out, any))>=1,1,0)
      
      tdf4 <- c("41L", "67N", "70R", "210W", "215A/C/D/E/G/H/I/L/N/S/V/F/Y","219E/H/N/Q/R")
      tdf4_out <- map(tdf4, mut)
      tdf4_len <- sum(sapply(tdf4_out, any))
      
      tdf5 <- c("70E/G")
      tdf5_out <- map(tdf5, mut)
      tdf5_len <- ifelse(sum(sapply(tdf5_out, any))>=1,1,0)
      
      tdf6 <- c("41L", "210W", "215Y")
      tdf6_len <- sum(tdf6%in%mut)
      
      tdf7 <- c("151M")
      tdf7_len <- sum(tdf7%in%mut)
      
      tdf8 <- c("75I","77L","116Y")
      tdf8_len <- sum(tdf8%in%mut)
      
      if(tdf1_len == 1 | (tdf2_len == 1 & tdf3_len == 0) | tdf4_len >= 5){
        tdfclass <- "Resistant"
      } else if (tdf5_len == 1 | (tdf2_len == 1 & tdf3_len == 1) | tdf6_len >= 3 | tdf4_len >= 4 | (tdf7_len == 1 & tdf8_len >= 3)){
        tdfclass <- "Intermediate Resistant"
      } else {
        tdfclass <- "Susceptible"
      }
      
      #d4T
      
      d4tclass <- NA
      
      d4t1 <- c("67d","69i","75A/M/S/T","151M")
      d4t1_out <- map(d4t1, mut)
      d4t1_len <- sum(sapply(d4t1_out, any))
      
      d4t2 <- c("41L","67N","69A/D/G/N","70R","210W","215A/C/D/E/G/H/I/L/N/S/V/F/Y","219E/H/N/Q/R")
      d4t2_out <- map(d4t2, mut)
      d4t2_len <- sum(sapply(d4t2_out, any))
      
      d4t3 <- c("41L","210W","215Y")
      d4t3_len <- sum(d4t3%in%mut)
      
      d4t4 <- c("215A/C/D/E/G/H/I/L/N/S/V/F/Y")
      d4t4_out <- map(d4t4, mut)
      d4t4_len <- ifelse(sum(sapply(d4t4_out, any))>=1,1,0)
      
      d4t5 <- c("184I/V")
      d4t5_out <- map(d4t5, mut)
      d4t5_len <- ifelse(sum(sapply(d4t5_out, any))>=1,1,0)
      
      if(d4t1_len >= 1 | d4t2_len >= 4 | d4t3_len == 3){
        d4tclass <- "Resistant"
      } else if ((d4t2_len >= 2 & d4t2_len <= 3)| (d4t4_len == 1 & d4t5_len == 0)){
        d4tclass <- "Intermediate Resistant"
      } else {
        d4tclass <- "Susceptible"
      }
      
      #FTC
      
      ftcclass <- NA
      
      ftc1 <- c("184I/V")
      ftc1_out <- map(ftc1, mut)
      ftc1_len <- ifelse(sum(sapply(ftc1_out, any))>=1,1,0)
      
      ftc2 <- c("65N/R", "151M")
      ftc2_out <- map(ftc2, mut)
      ftc2_len <- sum(sapply(ftc2_out, any))
      
      ftc3 <- c("65N/R","67d", "69i", "70E/G")
      ftc3_out <- map(ftc3, mut)
      ftc3_len <- sum(sapply(ftc3_out, any))
      
      ftc4 <- c("44A/D","118I")
      ftc4_out <- map(ftc4, mut)
      ftc4_len <- sum(sapply(ftc4_out, any))
      
      ftc5 <- c("41L","67N", "69A/N", "70R","210W","215A/C/D/E/G/H/I/L/N/S/V/F/Y","219E/H/N/Q/R")
      ftc5_out <- map(ftc5, mut)
      ftc5_len <- sum(sapply(ftc5_out, any))
      
      ftc6 <- c("151M")
      ftc6_len <- sum(ftc6%in%mut)
      
      ftc7 <- c("75I","77L","116Y")
      ftc7_len <- sum(ftc7%in%mut)
      
      if(ftc1_len == 1 | ftc2_len == 2){
        ftcclass <- "Resistant"
      } else if (ftc3_len >= 1| (ftc4_len >= 1 & ftc5_len >= 3) | (ftc6_len == 1 & ftc7_len >= 3)){
        ftcclass <- "Intermediate Resistant"
      } else {
        ftcclass <- "Susceptible"
      }
      
      #EFV
      
      efvclass <- NA
      
      efv1 <- c("100I", "101P", "103N/S/T/H","106M", "181C/T/V","188C/H/F/L", "190A/S/E/Q/C/T/V","230I/L")
      efv1_out <- map(efv1, mut)
      efv1_len <- sum(sapply(efv1_out, any))*2
      
      efv2 <- c("106A","138Q","225H")
      efv2_len <- sum(efv2%in%mut)*1.5
      
      efv3 <- c("90I", "98G", "101E/H/N/Q", "103R", "106I", "108I", "138K", 
                "179D/E", "221Y", "227C/L", "238N/T", "318F")
      efv3_out <- map(efv3, mut)
      efv3_len <- sum(sapply(efv3_out, any))*0.5
      
      scoreefv <- efv1_len + efv2_len + efv3_len
      
      
      if(scoreefv >= 2){
        efvclass <- "Resistant"
      } else if (scoreefv >= 1 & scoreefv < 2){
        efvclass <- "Intermediate Resistant" 
      } else {
        efvclass <- "Susceptible"
      }
      
      #NVP
      
      nvpclass <- NA
      
      nvp1 <- c("100I", "101P", "103N/S/T/H", "106A/M", "181C/I/V", "188C/H/L", "190A/S/E/Q/C/T/V", "230I/L")
      nvp1_out <- map(nvp1, mut)
      nvp1_len <- sum(sapply(nvp1_out, any))*2
      
      nvp2 <- c("138K/Q", "227C/L","238N/T","318F")
      nvp2_out <- map(nvp2, mut)
      nvp2_len <- sum(sapply(nvp2_out, any))*1.5
      
      nvp3 <- c("98G","101E/H/N/Q","103R","106I","108I","179D/E","221Y")
      nvp3_out <- map(nvp3, mut)
      nvp3_len <- sum(sapply(nvp3_out, any))*0.5
      
      scorenvp <- nvp1_len + nvp2_len + nvp3_len
      
      if(scorenvp >= 2){
        nvpclass <- "Resistant"
      } else if (scorenvp >= 1 & scorenvp < 2){
        nvpclass <- "Intermediate Resistant"
      } else {
        nvpclass <- "Susceptible"
      }
      
      #ETR
      
      etrclass <- NA
      
      etr1 <- c("100I","101P","181C/I/V","227C","230I/L")
      etr1_out <- map(etr1, mut)
      etr1_len <- sum(sapply(etr1_out, any))*1.5
      
      etr2 <- c("138K/Q","179F","190C/E/Q/S/T/V")
      etr2_out <- map(etr2, mut)
      etr2_len <- sum(sapply(etr2_out, any))*0.75
      
      etr3 <- c("101E/H/N/Q","106I","138A","190A")
      etr3_out <- map(etr3, mut)
      etr3_len <- sum(sapply(etr3_out, any))*0.5
      
      etr4 <- c("90I","98G","101R","103H/N/T","106A/M","108I","179D/E/T","188C/H/F/L", "221Y","225H","227L","234I","236L","238N/T", "318F")
      etr4_out <- map(etr4, mut)
      etr4_len <- sum(sapply(etr4_out, any))*0.25
      
      scoreetr <- etr1_len + etr2_len + etr3_len + etr4_len
      
      if(scoreetr >= 3){
        etrclass <- "Resistant"
      } else if (scoreetr >= 1.5 & scoreetr < 3){
        etrclass <- "Intermediate Resistant"
      } else {
        etrclass <- "Susceptible"
      }
      
      #RPV
      
      rpvclass <- NA
      
      rpv1 <- c("101P", "181I/V", "188L")
      rpv1_out <- map(rpv1, mut)
      rpv1_len <- sum(sapply(rpv1_out, any))
      
      rpv2 <- c("101E")
      rpv2_len <- sum(rpv2%in%mut)
      
      rpv3 <- c("184I/V")
      rpv3_out <- map(rpv3, mut)
      rpv3_len <- ifelse(sum(sapply(rpv3_out, any))>=1,1,0)
      
      rpv4 <- c("100I")
      rpv4_len <- sum(rpv4%in%mut)
      
      rpv5 <- c("103N")
      rpv5_len <- sum(rpv5%in%mut)
      
      rpv6 <- c("138K")
      rpv6_len <- sum(rpv6%in%mut)
      
      rpv7 <- c("101E","138A/G/K/Q/R/S", "179L","181C","221Y","227C","230I/L")
      rpv7_out <- map(rpv7, mut)
      rpv7_len <- sum(sapply(rpv7_out, any))*1.5
      
      rpv8 <- c("90I","101Q/T", "103S","106A/I", "108I","179D/I/T", "189I", "190E")
      rpv8_out <- map(rpv8, mut)
      rpv8_len <- sum(sapply(rpv8_out, any))*0.5
      
      if(rpv1_len>= 1 | ((rpv2_len + rpv3_len) == 2) | ((rpv4_len + rpv5_len) == 2) | ((rpv3_len + rpv6_len) == 2)){
        s_rpv1 <- 2.5
      } else {
        s_rpv1 <- 0
      }
      
      scorerpv <- s_rpv1 + rpv7_len + rpv8_len
      
      if(scorerpv >= 2.5){
        rpvclass <- "Resistant"
      } else if (scorerpv >= 1.5 & scorerpv < 2.5){
        rpvclass <- "Intermediate Resistant"
      } else {
        rpvclass <- "Susceptible"
      }
      
      dorclass <- "Not Present"
      
      res_rega <- rbind.data.frame(c("doravirine (DOR)", dorclass),
                                   c("efavirenz (EFV)", efvclass),
                                   c("etravirine (ETR)", etrclass),
                                   c("nevirapine (NVP)", nvpclass),
                                   c("rilpivirine (RPV)", rpvclass),
                                   c("abacavir (ABC)", abcclass),
                                   c("didanosine (DDI)", ddiclass),
                                   c("emtricitabine (FTC)", tc3class),
                                   c("lamivudine (3TC)", tc3class),
                                   c("stavudine (D4T)", d4tclass),
                                   c("tenofovir (TDF)", tdfclass),
                                   c("zidovudine (AZT)", aztclass))
      
      colnames(res_rega) <- c("drug", "Rega")
      
      ##### Brazilian Algorithm #####
      
      mut <- paste0(posi, aa)
      
      #EFV
      
      efvclass <- NA
      
      efv1 <- c("100I", "101E/P/Q", "103N/A/S/T/Q/H", "106A/M", "181C/I/V,188C/H/L", "190A/S/E/Q/C/T/V",
                "225H", "230L")
      efv1_out <- map(efv1, mut)
      efv1_len <- sum(sapply(efv1_out, any))
      
      if(efv1_len >= 1){
        efvclass <- "Resistant"
      } else {
        efvclass <- "Susceptible"
      }
      
      #NVP
      
      nvpclass <- NA
      
      nvp1 <- c("98G", "100I", "101E/P/Q", "103N/A/S/T/Q/H",
                "106A/M", "108I", "179D/E/F", "181C/I/V",
                "188C/H/L", "190A/S/E/Q/C/T/V", "227L/C", "230L")
      nvp1_out <- map(nvp1, mut)
      nvp1_len <- sum(sapply(nvp1_out, any))
      
      if(nvp1_len >= 1){
        nvpclass <- "Resistant"
      } else {
        nvpclass <- "Susceptible"
      }
      
      #ETR
      
      etrclass <- NA
      
      etr1 <- c("181I/V/C", "100I","101P","230L")
      etr1_out <- map(etr1, mut)
      etr1_len <- sum(sapply(etr1_out, any))
      
      etr2 <- c("90I","98G","101E/H","106I","138A/G/K/Q","179F/T/D","190A/S")
      etr2_out <- map(etr2, mut)
      etr2_len <- sum(sapply(etr2_out, any))
      
      if(etr1_len >= 2 | (etr1_len >= 1 & etr2_len >= 1) | etr2_len >= 3){
        etrclass <- "Resistant"
      } else if (etr2_len == 2 | etr1_len == 1){
        etrclass <- "Intermediate Resistant"
      } else {
        etrclass <- "Susceptible"
      }
      
      #3TC
      
      tc3class <- NA
      
      tc31 <- c("69i","151M/L", "67d")
      tc31_out <- map(tc31, mut)
      tc31_len <- sum(sapply(tc31_out, any))
      
      tc32 <- c("184V/I")
      tc32_out <- map(tc32, mut)
      tc32_len <- ifelse(sum(sapply(tc32_out, any))>1, 1,0)
      
      tc33 <- c("65R/N")
      tc33_out <- map(tc33, mut)
      tc33_len <- ifelse(sum(sapply(tc33_out, any))>1, 1,0)
      
      tc34 <- c("44A/D","118I")
      tc34_out <- map(tc34, mut)
      tc34_len <- sum(sapply(tc34_out, any))
      
      tc35 <- c("41L", "67N/E/G", "70R/G/N", "210W", "215Y/F","219E/Q/N/R")
      tc35_out <- map(tc35, mut)
      tc35_len <- sum(sapply(tc35_out, any))
      
      if(tc31_len >= 1 | tc32_len == 1){
        tc3class <- "Resistant"
      } else if ((tc33_len == 1 | tc34_len == 2) | tc35_len >= 3){
        tc3class <- "Intermediate Resistant"
      } else {
        tc3class <- "Susceptible"
      }
      
      #ABC
      
      abcclass <- NA
      
      abc1 <- c("69i","151M/L", "67d")
      abc1_out <- map(abc1, mut)
      abc1_len <- sum(sapply(abc1_out, any))
      
      abc2 <- c("184V/I")
      abc2_out <- map(abc2, mut)
      abc2_len <- ifelse(sum(sapply(abc2_out, any))>1,1,0)
      
      abc3 <- c("65R","74V/I","115F")
      abc3_out <- map(abc3, mut)
      abc3_len <- sum(sapply(abc3_out, any))
      
      abc4 <- c("41L", "67N/E/G", "70R/G/N","210W", "215Y/F","219E/Q/N/R")
      abc4_out <- map(abc4, mut)
      abc4_len <- sum(sapply(abc4_out, any))
      
      if(abc1_len >= 1 | (abc2_len == 1 & abc3_len >= 1) | (abc2_len == 1 & abc4_len >= 3)){
        abcclass <- "Resistant"
      } else if (abc3_len >= 1 | abc4_len >= 4 | abc2_len == 1){
        abcclass <- "Intermediate Resistant"
      } else {
        abcclass <- "Susceptible"
      }
      
      #AZT
      
      aztclass <- NA
      
      azt1 <- c("69i","151M/L", "67d")
      azt1_out <- map(azt1, mut)
      azt1_len <- sum(sapply(azt1_out, any))
      
      azt2 <- c("215Y/F/C/D/S/I/E/N/V")
      azt2_out <- map(azt2, mut)
      azt2_len <- sum(sapply(azt2_out, any))
      
      azt3 <- c("40F", "41L", "67N/E/G", "70R/G/N","210W","219Q/E/N/R")
      azt3_out <- map(azt3, mut)
      azt3_len <- sum(sapply(azt3_out, any))
      
      if(azt1_len >= 1 | (azt2_len >= 1 & azt3_len >= 1)){
        aztclass <- "Resistant"
      } else if (azt3_len >= 2 | azt2_len >= 1){
        aztclass <- "Intermediate Resistant"
      } else {
        aztclass <- "Susceptible"
      }
      
      #ddI
      
      ddiclass <- NA
      
      ddi1 <- c("65R", "74V/I")
      ddi1_out <- map(ddi1, mut)
      ddi1_len <- sum(sapply(ddi1_out, any))
      
      ddi2 <- c("41L", "67N/E/G", "70R/G/N", "210W", "215Y/F","219E/Q/N/R")
      ddi2_out <- map(ddi2, mut)
      ddi2_len <- sum(sapply(ddi2_out, any))
      
      ddi3 <- c("184V/I")
      ddi3_out <- map(ddi3, mut)
      ddi3_len <- ifelse(sum(sapply(ddi3_out, any))>1,1,0)
      
      ddi4 <- c("69D/N","75T","184V/I", "70E")
      ddi4_out <- map(ddi4, mut)
      ddi4_len <- sum(sapply(ddi4_out, any))
      
      ddi5 <- c("69i","151M/L","67d")
      ddi5_out <- map(ddi5, mut)
      ddi5_len <- sum(sapply(ddi5_out, any))
      
      if(ddi5_len >= 1 | ddi1_len >= 1 | (ddi2_len >= 3 & ddi3_len == 1)){
        ddiclass <- "Resistant"
      } else if (ddi4_len >= 1 | ddi2_len >= 3){
        ddiclass <- "Intermediate Resistant"
      } else {
        ddiclass <- "Susceptible"
      }
      
      #TDF
      
      tdfclass <- NA
      
      tdf1 <- c("65R", "69i", "151M/L", "67d", "70E")
      tdf1_out <- map(tdf1, mut)
      tdf1_len <- sum(sapply(tdf1_out, any))
      
      tdf2 <- c("41L", "67N", "70R/G/N", "210W", "215Y/F","219E/Q/N/R")
      tdf2_out <- map(tdf2, mut)
      tdf2_len <- sum(sapply(tdf2_out, any))
      
      tdf3 <- c("41L","210W")
      tdf3_len <- sum(tdf3%in%mut)
      
      tdf4 <- c("115F")
      tdf4_len <- sum(tdf4%in%mut)
      
      if(tdf1_len >= 1 | (tdf2_len >= 3 & tdf3_len >= 1)){
        tdfclass <- "Resistant"
      } else if (tdf4_len == 1){
        tdfclass <- "Intermediate Resistant"
      } else {
        tdfclass <- "Susceptible"
      }
      
      rvpclass <- "Not Present"
      ftcclass <- "Not Present"
      d4tclass <- "Not Present"
      dorclass <- "Not Present"
      
      res_bra <- rbind.data.frame(c("doravirine (DOR)", dorclass),
                                  c("efavirenz (EFV)", efvclass),
                                  c("etravirine (ETR)", etrclass),
                                  c("nevirapine (NVP)", nvpclass),
                                  c("rilpivirine (RPV)", rpvclass),
                                  c("abacavir (ABC)", abcclass),
                                  c("didanosine (DDI)", ddiclass),
                                  c("emtricitabine (FTC)", tc3class),
                                  c("lamivudine (3TC)", tc3class),
                                  c("stavudine (D4T)", d4tclass),
                                  c("tenofovir (TDF)", tdfclass),
                                  c("zidovudine (AZT)", aztclass))
      
      colnames(res_bra) <- c("drug", "Brazilian Algorithm")
      
      
      fim <- cbind.data.frame(totfim[1], totfim[3], res_rega[2], res_anrs[2], res_bra[2]) 
      
      clasfim <- fim  
      
    }
    
    if(input$region == "Integrase"){
      
      Scores_INI <- read.table(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/HIVdb/Scores_INI_2018.txt"), header = T)
      
      drugName <- c("BIC","BIC","BIC","BIC","BIC","BIC","BIC","BIC","BIC","BIC","BIC","BIC","BIC",
                    "DTG","DTG","DTG","DTG","DTG","DTG","DTG","DTG","DTG","DTG","DTG","DTG","DTG",
                    "EVG","EVG","EVG","EVG","EVG","EVG","EVG", 
                    "RAL","RAL")
      
      p1 <- c(138,138,140,143,143,148,148,157,51,74,74,92,97,
              138,138,140,143,143,148,148,157,51,74,74,92,97, 
              138,143,143,51,74,74,97, 
              138,74)
      
      a1 <- c("AKT","AKT","ACS","ACGHRS","ACGHRS","HKR","HRK","Q","Y","FM","FM","Q","A",
              "AKT","AKT","ACS","ACGHRS","ACGHRS","HKR","HRK","Q","Y","FM","FM","Q","A", 
              "AKT","ACGHRS","ACGHRS","Y","FM","FM","A", 
              "AKT","FM")
      
      p2 <- c(140,148,148,163,230,155,163,263,263,143,148,155,148,
              140,148,148,163,230,155,163,263,263,143,148,155,148, 
              140,163,230,263,143,148,143, 
              140,148)
      
      a2 <- c("ACS","HKR","HKR","R","R","H","KR","K","K","ACGHRS","HKR","H","HKR",
              "ACS","HKR","HKR","R","R","H","KR","K","K","ACGHRS","HKR","H","HKR", 
              "ACS","R","R","K","ACGHRS","HKR","ACGHRS", 
              "ACS","HKR")
      
      score <- c(10,10,10,5,5,10,5,10,10,5,10,5,10,
                 10,10,10,5,5,10,5,10,10,5,10,5,10, 
                 15,5,5,15,5,10,5, 
                 15,10)           
      
      comb_INI <- cbind.data.frame(drugName, p1, a1, p2, a2, score)
      
      if (input$hivdb == ">=20%"){
        
        cutoff <- subset(dposmut, freqf >= 0.2)
        
      }
      
      if (input$hivdb == ">=1%"){
        
        cutoff <- subset(dposmut, freqf >= 0.01)
      }
      
      ##### HIVdb #####
      
      pa <- with(Scores_INI, paste0(Position, AA))
      pro <- cbind(pa, Scores_INI)
      aa <- cutoff$aa
      posi <- cutoff$pos
      mut <- paste0(posi, aa)
      new <- pro[pro$pa %in% mut,]
      res <- colSums(new[,4:dim(new)[2]]) #punctuation simple resistance
      z <- data.frame(names(res), res) #reorganizing the simple answer
      row.names(z) <- NULL
      names(z) <- c("drugName", "score")
      
      cc <- data.frame()
      
      tam <- length(new$pa) #items number found in table
      
      #punctuation combined resistance
      library(gtools)
      if (tam >= 3){
        for(k in 1:tam){
          ind <- combn(seq(1:tam), 2)[,k] #to generate the combinations between the positions
          #combination
          a <- subset(comb_INI, comb_INI[,2] == new$Position[ind[1]])
          b <- a[grep(new$AA[ind[1]], a[,3]),] 
          c <- subset(b, b[,4] == new$Position[ind[2]])
          d <- c[grep(new$AA[ind[2]], c[,5]),]
          
          cc <- rbind(cc, d)
        }
        
        if (dim(cc)[1] == 0){ #in case there is no combination, the data set will be empty
          cfim <- z
          
        } else {
          cc1 <- cc[with(cc, order(drugName)), ] #orders the drugs to stay close to each other
          cc2 <- cc1[,c("drugName", "score")] 
          cc3 <- aggregate(. ~ drugName, data=cc2, FUN=sum) #sum the values of the same drugs
          
          cfim <- aggregate(. ~ drugName, data=rbind(cc3, z), FUN=sum)
        }
      }
      
      
      if (tam == 2){
        a <- subset(comb_INI, comb_INI[,2] == new$Position[1])
        b <- a[grep(new$AA[1], a[,3]),] 
        c <- subset(b, b[,4] == new$Position[2])
        d <- c[grep(new$AA[2], c[,5]),]
        
        cc <- d
        
        if (dim(cc)[1] == 0){ #in case there is no combination, the data set will be empty
          cfim <- z
          
        } else {
          cc1 <- cc[with(cc, order(drugName)), ] #orders the drugs to stay close to each other
          cc2 <- cc1[,c("drugName", "score")] 
          cc3 <- aggregate(. ~ drugName, data=cc2, FUN=sum) #sum the values of the same drugs
          
          cfim <- aggregate(. ~ drugName, data=rbind(cc3, z), FUN=sum)
        }
      }
      
      if (tam <= 1){
        cfim <- z
      }
      
      cfim$class <- rep(NA, nrow(cfim)) #added a column in dataframe "cfim"
      
      #Stanford (HIVdb) classifications
      
      for (f in 1:dim(cfim)[1]){
        cfim$score[f]
        if (cfim$score[f] <= 9){cfim$class[f] <- "Susceptible"}
        if (cfim$score[f] >= 10 & cfim$score[f] <= 14){cfim$class[f] <- "Potential low-level resistance"}
        if (cfim$score[f] >= 15 & cfim$score[f] <= 29){cfim$class[f] <- "Low-level resistance"}
        if (cfim$score[f] >= 30 & cfim$score[f] <= 59){cfim$class[f] <- "Intermediate resistance"}
        if (cfim$score[f] >= 60){cfim$class[f] <- "High-level resistance"}
      }
      
      cfim$drugName <- gsub("RAL", "raltegravir (RAL)", cfim$drugName)
      cfim$drugName <- gsub("EVG", "elvitegravir (EVG)", cfim$drugName)
      cfim$drugName <- gsub("DTG", "dolutegravir (DTG)", cfim$drugName)
      cfim$drugName <- gsub("BIC", "bictegravir (BIC)", cfim$drugName)
      
      names(cfim) <- c("Integrase Inhibitors", "Mutation Scoring", "HIVdb")
      
      res_hivdb <- cfim
      
      ##### ANRS #####
      
      mut <- paste0(posi, aa)
      
      #RAL
      
      ralclass <- NA
      
      ral1 <- c("66K", "92Q", "118R", "121Y", "140A/S", "143A/C/G/H/R/S", "148E/G/H/K/R", "151L", "155H/S/T", "157Q", "230R", "263K")
      ral1_out <- map(ral1, mut)
      ral1_len <- sum(sapply(ral1_out, any))
      
      ral2 <- c("74F")
      ral2_len <- sum(sapply(mut, function(x) grepl(x, ral2)))
      
      ral3 <- c("75I")
      ral3_len <- sum(sapply(mut, function(x) grepl(x, ral3)))
      
      scoreral <- ral2_len + ral3_len 
      
      if(ral1_len >= 1 | scoreral == 2){
        ralclass <- "Resistant"
      } else {
        ralclass <- "Susceptible"
      }
      
      #EVG
      
      evgclass <- NA
      
      evg1 <- c("66I/A/K", "92Q", "97A","118R", "121Y", "138K", "140C/S", "143A/C/G/H/R/S", "145S", "147G", 
                "148H/R/K", "145S", "147G", "148H/R/K","151L", "155H/S/T","157Q", "230R","263K")
      evg1_out <- map(evg1, mut)
      evg1_len <- sum(sapply(evg1_out, any))
      
      evg2 <- c("74F")
      evg2_len <- sum(sapply(mut, function(x) grepl(x, evg2)))
      
      evg3 <- c("75I")
      evg3_len <- sum(sapply(mut, function(x) grepl(x, evg3)))
      
      scoreevg <- evg2_len + evg3_len 
      
      if(evg1_len >= 1 | scoreevg == 2){
        evgclass <- "Resistant"
      } else {
        evgclass <- "Susceptible"
      }
      
      
      #DTG
      
      dtgclass <- NA
      
      dtg1 <- c("118R","121Y", "151L", "153F/Y", "263K")
      dtg1_out <- map(dtg1, mut)
      dtg1_len <- sum(sapply(dtg1_out, any))
      
      dtg2 <- c("66K")
      dtg2_len <- sum(sapply(mut, function(x) grepl(x, dtg2)))
      
      dtg3 <- c("74M")
      dtg3_len <- sum(sapply(mut, function(x) grepl(x, dtg3)))
      
      dtg4 <- c("92Q")
      dtg4_len <- sum(sapply(mut, function(x) grepl(x, dtg4)))
      
      dtg5 <- c("155H")
      dtg5_len <- sum(sapply(mut, function(x) grepl(x, dtg5)))
      
      dtg6 <- c("148H/K/R")
      dtg6_out <- map(dtg6, mut)
      dtg6_len <- ifelse(sum(sapply(dtg6_out, any))>=1,1,0)
      
      dtg7 <- c("74I", "138A/K/T", "140A/C/S")
      dtg7_out <- map(dtg7, mut)
      dtg7_len <- sum(sapply(dtg7_out, any))
      
      dtg8 <- c("157Q")
      dtg8_len <- sum(sapply(mut, function(x) grepl(x, dtg8)))
      
      scoredtg1 <- dtg2_len + dtg3_len
      scoredtg2 <- dtg4_len + dtg5_len
      scoredtg3 <- dtg6_len + dtg5_len
      
      if(dtg1_len >= 1 | (dtg6_len == 1 & dtg7_len >= 2) | scoredtg1 == 2 | scoredtg2 == 2 | scoredtg3 == 2){
        dtgclass <- "Resistant"
      } else if (dtg2_len == 1 | dtg8_len == 1 | (dtg6_len == 1 & dtg7_len == 1)){
        dtgclass <- "Possible Resistant"
      } else {
        dtgclass <- "Susceptible"
      }
      
      #BIC
      
      bicclass <- NA
      
      bic1 <- c("118R","121Y","138A/K/T", "140A/C/S", "148H/K/R", "151L", "153F/Y", "155H","230R","263K")
      bic1_out <- map(bic1, mut)
      bic1_len <- sum(sapply(mut, function(x) grepl(x, bic1)))
      
      bic2 <- c("66K")
      bic2_len <- sum(sapply(mut, function(x) grepl(x, bic2)))
      
      bic3 <- c("74M")
      bic3_len <- sum(sapply(mut, function(x) grepl(x, bic3)))
      
      bic4 <- c("74I")
      bic4_len <- sum(sapply(mut, function(x) grepl(x, bic4)))
      
      bic5 <- c("92Q")
      bic5_len <- sum(sapply(mut, function(x) grepl(x, bic5)))
      
      bic6 <- c("157Q")
      bic6_len <- sum(sapply(mut, function(x) grepl(x, bic6)))
      
      scorebic1 <- bic2_len + bic3_len
      scorebic2 <- bic4_len + bic5_len
      
      if(bic1_len >= 1 | scorebic1 == 2 | scorebic2 == 2){
        bicclass <- "Resistant"
      } else if (bic2_len == 1 | bic6_len == 1){
        bicclass <- "Possible Resistant"
      } else {
        bicclass <- "Susceptible"
      }
      
      res_anrs <- rbind.data.frame(c("bictegravir (BIC)", bicclass),
                                   c("dolutegravir (DTG)", dtgclass),
                                   c("elvitegravir (EVG)", evgclass),
                                   c("raltegravir (RAL)", ralclass))
      
      colnames(res_anrs) <- c("drug", "ANRS")
      
      ###### Rega ######
      
      mut <- paste0(posi, aa)
      
      #RAL
      
      ralclass <- NA
      
      ral1 <- c("143R", "148H/K/R", "155H/S")
      ral1_out <- map(ral1, mut)
      ral1_len <- sum(sapply(ral1_out, any))*20
      
      ral2 <- c("74M","92Q","138A","140A/S", "143C/H", "163K","232N")
      ral2_out <- map(ral2, mut)
      ral2_len <- sum(sapply(ral2_out, any))*1
      
      ral3 <- c("97A","138K","151I","156N","163R","206S","230N")
      ral3_len <- sum(sapply(mut, function(x) grepl(x, ral3)))*0.25
      
      scoreral <- ral1_len + ral2_len + ral3_len
      
      if(scoreral >= 20){
        ralclass <- "Resistant"
      } else if (scoreral >= 1 & scoreral < 20){
        ralclass <- "Intermediate Resistant"
      } else {
        ralclass <- "Susceptible"
      }
      
      #EVG
      
      evgclass <- NA
      
      evg1 <- c("66A/I/K","92Q/G","145S","146I/L/P/R","147G","148H/K/R","155H/S/T")
      evg1_out <- map(evg1, mut)
      evg1_len <- sum(sapply(evg1_out, any))*2
      
      evg2 <- c("121Y","140A/C/S")
      evg2_out <- map(evg2, mut)
      evg2_len <- sum(sapply(evg2_out, any))*1
      
      evg3 <- c("51Y","97A","114Y","146K","153F/Y","263K")
      evg3_out <- map(evg3, mut)
      evg3_len <- sum(sapply(evg3_out, any))*0.5
      
      scoreevg <- evg1_len + evg2_len + evg3_len
      
      if(scoreevg >= 2){
        evgclass <- "Resistant"
      } else if (scoreevg >= 1 & scoreevg < 2){
        evgclass <- "Intermediate Resistant"
      } else {
        evgclass <- "Susceptible"
      }
      
      #DTG
      
      dtgclass <- NA
      
      dtg1 <- c("118R","263K")
      dtg1_len <- sum(sapply(mut, function(x) grepl(x, dtg1)))*1.5
      
      dtg2 <- c("148H/K/R")
      dtg2_out <- map(dtg2, mut)
      dtg2_len <- ifelse(sum(sapply(dtg2_out, any))>=1,1,0)*1
      
      dtg3 <- c("51Y","74M","92Q","97A","124A","138K","140A/C/S","143C/H/R","153F","155H")
      dtg3_out <- map(dtg3, mut)
      dtg3_len <- sum(sapply(dtg3_out, any))*0.5
      
      scoredtg <- dtg1_len + dtg2_len + dtg3_len
      
      if(scoredtg >= 2){
        dtgclass <- "Resistant"
      } else if (scoredtg >= 1.5 & scoredtg < 2){
        dtgclass <- "Intermediate Resistant"
      } else {
        dtgclass <- "Susceptible"
      }
      
      bicclass <- "Not present"
      
      res_rega <- rbind.data.frame(c("bictegravir (BIC)", bicclass),
                                   c("dolutegravir (DTG)", dtgclass),
                                   c("elvitegravir (EVG)", evgclass),
                                   c("raltegravir (RAL)", ralclass))
      
      colnames(res_rega) <- c("drug", "Rega")
      
      ##### Brazilian Algorithm #####
      
      dtgclass <- "Not Present"
      evgclass <- "Not Present"
      ralclass <- "Not present"
      bicclass <- "Not present"
      
      res_bra <- rbind.data.frame(c("bictegravir (BIC)", bicclass),
                                  c("dolutegravir (DTG)", dtgclass),
                                  c("elvitegravir (EVG)", evgclass),
                                  c("raltegravir (RAL)", ralclass))
      
      colnames(res_bra) <- c("drug", "Brazilian Algorithm")
      
      fim <- cbind.data.frame(res_hivdb[1], res_hivdb[3], res_rega[2], res_anrs[2], res_bra[2])
      
      clasfim <- fim
      
    }
    
    clasfim
    
  })
  
  output$table2 <- DT::renderDataTable({
    d <- datahivdb()
    dd <- cbind(d[, 1], d[, input$show_vars, drop = F])
    names(dd)[1] <- paste(input$region, "Inhibitors")
    
    table2 <-
      datatable(
        dd,
        extensions = 'Buttons',
        options = list(
          pageLength = nrow(dd),
          searching = T,
          dom = 'lfBrtip',
          buttons =
            list(list(
              extend = 'collection',
              buttons = list(
                list(
                  extend = 'csv',
                  filename = paste(
                    input$text,
                    "-",
                    'Frequency:',
                    input$hivdb,
                    '_',
                    input$region,
                    sep = ''
                  ), title = "Genotypic Resistance Interpretation Algorithms"
                ),
                list(
                  extend = 'excel',
                  filename = paste(
                    input$text,
                    "-",
                    'Frequency:',
                    input$hivdb,
                    '_',
                    input$region,
                    sep = ''
                  ), title = "Genotypic Resistance Interpretation Algorithms"
                ),
                list(
                  extend = 'pdf',
                  filename = paste(
                    input$text,
                    "-",
                    'Frequency:',
                    input$hivdb,
                    '_',
                    input$region,
                    sep = ''
                  ), title = "Genotypic Resistance Interpretation Algorithms", messageTop = paste0 ("This table shows the HIV drug susceptibility to ", input$region, " inhibitors in according to rule-based algorithms chosen."), messageBottom = "All sequence interpretations are for research use only."
                  
                ),
                list(
                  extend = 'print',
                  filename = paste(
                    input$text,
                    "-",
                    'Frequency:',
                    input$hivdb,
                    '_',
                    input$region,
                    sep = ''
                  ), title = "Genotypic Resistance Interpretation Algorithms", messageTop = paste0 ("This table shows the HIV drug susceptibility to ", input$region, " inhibitors in according to rule-based algorithms chosen."), messageBottom = "All sequence interpretations are for research use only."
                )
              ),
              text = 'Download'
            )),
          scrollX = TRUE,
          initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'background-color': '#2c3e50', 'color': '	#ffffff'});",
            "}"
          )
        )
      ) %>% formatStyle(
        colnames(dd),
        backgroundColor = styleEqual(c("Susceptible", "Potential low-level resistance", "Low-level resistance", "Intermediate resistance", "Intermediate Resistant", "Possible Resistant", "High-level resistance", "Resistant"),
                                     c('#addd8e', '#ffffcc', '#fed976', '#fd8d3c', '#fd8d3c', '#fd8d3c', '#f03b20', '#f03b20'))
      )
    
  })
}

ui <- fluidPage(theme = "bootstrap.css",
                  list(tags$head(HTML('<link rel="icon", href="sira.png", 
                                      type="image/png" />'))),
                  div(style="padding: 1px 0px; width: '100%'",
                      titlePanel(
                        title="", windowTitle="SIRA-HIV"
                      )
                  ),
                  
                  # titlePanel("Application Title"),
                  
                  navbarPage(
                    
                    title=strong("SIRA-HIV"),
                    
                    tabPanel("Home", 
                             fluidPage(
                               br(),
                               fluidRow(
                                 column(12,
                                        
                                        h4("SIRA-HIV", align = "center", style = "font-weight: 1000; font-size: 50px; line-height: 3"),
                                        h3("A web application for identifying HIV drug resistance from a next generation sequencing data", align = "center"),
                                        br(),
                                        p("- SIRA-HIV is a tool that can detect minority mutations in HIV sequences with a 1% cut-off."),
                                        p("- Segminator II is used for the characterization of viral read data. This system supports the next generation sequencing platforms including: 454 Life Sciences, Illumina, Ion Torrent and Pacific Biosciences."),
                                        p("- The output of SIRA-HIV is a drug resistance report which provides the frequencies of the amino acids found in the regions of gene pol and the susceptibility of HIV-1 variants to antiretroviral drugs by genotypic interpretation systems."),
                                        tags$head(tags$style("#text1{color: red; font-size: 20px; font-style: italic;
                                                             }"
                         )
                                        )
                                        )))
                         
                                 ),
                    
                    tabPanel("Drug Resistance Positions",
                             sidebarLayout(fluid = TRUE,
                                           sidebarPanel(width = 4,
                                                        tags$head(tags$style(type="text/css", "
                                                                             #loadmessage {
                                                                             position: fixed;
                                                                             top: 0px;
                                                                             left: 0px;
                                                                             width: 100%;
                                                                             padding: 5px 0px 5px 0px;
                                                                             text-align: center;
                                                                             font-weight: bold;
                                                                             font-size: 100%;
                                                                             color: #000000;
                                                                             background-color: #4b5a6a;
                                                                             z-index: 105;
                                                                             }
                                                                             ")),
                                                        textInput("text", label = "Analysis Name (Optional):", value = "", placeholder = "Sample_ID"),
                                                        hr(),
                                                        h4(strong("Step 1 (Optional if you already have VEMETable file): Run Segminator II")),
                                                        tags$head(tags$script(src = "message-handler.js")),
                                                        actionButton("do", "Run"),
                                                        bsTooltip("do", "Run Segminator II for the characterization of viral read data from next generation sequencing platforms","right", options = list(container = "body")),
                                                        hr(),
                                                        fileInput('file1', 'Step 2: Choose VEMETable to upload:', accept = c('text/csv','text/comma-separated-values','text/tab-separated-values','text/plain','.csv','.tsv')),
                                                        bsTooltip("file1", "VEMETable is the file saved from Segminator II",
                                                                  "right", options = list(container = "body")),
                                                        # actionButton("loading", "Load an example"),
                                                        # bsTooltip("loading", "Load an example to see how SIRA-HIV works","right", options = list(container = "body")),
                                                        hr(),
                                                        selectInput("region", "Step 3: Choose a region of the pol gene:", choices = c("Protease", "Reverse Transcriptase", "Integrase"), selectize = T),
                                                        actionButton("Load", "Run"),
                                                        bsTooltip("Load", "Display the amino acid frequency table of the region of the pol gene selected",
                                                                  "right", options = list(container = "body")),
                                                        # h4("Download the results:"),
                                                        hr(),
                                                        actionButton("Loadd", "Coverage Plot"),
                                                        bsTooltip("Loadd", "Show the coverage plot of the region of the pol gene selected",
                                                                  "right", options = list(container = "body")),
                                                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                                         tags$div("Loading...",id="loadmessage"))
                                                        ),
                                           mainPanel(
                                             DT::dataTableOutput("table1"),
                                             HTML('<style>.rChart {width: 100%; height: 500px}</style>'),
                                             plotlyOutput("plot")
                                             
                                           )
                                                        )
                                           ),
                    
                    tabPanel("Genotypic Resistance Interpretation Algorithms",
                             sidebarLayout(fluid = TRUE,
                                           sidebarPanel(width = 4,
                                                        tags$head(tags$style(type="text/css", "
                                                                             #loadmessage {
                                                                             position: fixed;
                                                                             top: 0px;
                                                                             left: 0px;
                                                                             width: 100%;
                                                                             padding: 5px 0px 5px 0px;
                                                                             text-align: center;
                                                                             font-weight: bold;
                                                                             font-size: 100%;
                                                                             color: #000000;
                                                                             background-color: #4b5a6a;
                                                                             z-index: 105;
                                                                             }
                                                                             ")),
                                                        h2("Rule-based Algorithms"),
                                                        #h4("HIVdb"),
                                                        checkboxGroupInput("show_vars", label = h4(strong("Step 1: Choose the systems:")), 
                                                                           choices = list("ANRS" = "ANRS", "HIVdb" = "HIVdb", "REGA" = "Rega", "Brazilian Algorithm" = "Brazilian Algorithm"),
                                                                           selected = "ANRS"),
                                                        selectInput("hivdb", "Step 2: Choose cutoff:", choices = c(">=20%", ">=1%"), selectize =T),
                                                        actionButton("load", "Run"),
                                                        hr(),
                                                        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                                         tags$div("Loading...",id="loadmessage"))
                                                        
                                                        ),
                                           mainPanel(
                                             DT::dataTableOutput("table2"),
                                             HTML('<style>.rChart {width: 100%; height: 500px}</style>'))
                                                        )
                                           ),
                    
                    tabPanel("About",
                             fluidPage(
                               br(),
                               fluidRow(
                                 
                                 h2(strong("SIRA-HIV Scientific Board"), align = "center"),
                                 br(),
                                 h4(""),
                                 h4(""),
                                 br(),
                                 p("For any comments, suggestions or problems, please contact us."),
                                 br(),
                                 tags$b("Contact:"),
                                 tags$ul(
                                   tags$li("Leticia Martins Raposo (raposo@peb.ufrj.br)"), 
                                   tags$li("Flavio Fonseca Nobre (flavio@peb.ufrj.br)")
                                 ),
                                 # p("Federal University of Rio de Janeiro (UFRJ)"),
                                 # p("Alberto Luiz Coimbra Institute for Graduate Studies and Research in Engineering (COPPE)"),
                                 # p("Biomedical Engineering Program (PEB)"),
                                 # p("Rio de Janeiro, Brazil"),
                                 tags$img(src='contact.png', height = 130, width = 950),
                                 br(),
                                 br(),
                                 tags$b("Funding:"),
                                 p("This work was supported by the following organizations:"),
                                 tags$img(src='agency.png', height = 100, width = 700)
                                 
                               ))
                             
                    )
                    
                    )
                             )

shinyApp(ui = ui, server = server)