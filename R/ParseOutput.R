#' @export
parse_em_output <- function(prot, model, categories, toys = FALSE){
  emfiles <- list.files(path = file.path(getwd(), "output", "wood_likelihood", model, prot, fsep = .Platform$file.sep), pattern = paste(c(prot,"_[a-z]+_", categories,".*"), collapse=""), all.files = FALSE,
                        full.names = TRUE, recursive = FALSE)
  emshortfiles <- list.files(path = file.path(getwd(), "output", "wood_likelihood", model, prot, fsep = .Platform$file.sep), pattern = paste(c(prot,"_[a-z]+_", categories,".*"), collapse=""), all.files = FALSE,
                             full.names = FALSE, recursive = FALSE)
  if (toys){
    emfiles <- list.files(path = file.path(getwd(), "output", "toys", model, prot, fsep = .Platform$file.sep), pattern = paste(c(prot,"_[a-z]+_", categories,".*"), collapse=""), all.files = FALSE,
                          full.names = TRUE, recursive = FALSE)
    emshortfiles <- list.files(path = file.path(getwd(), "output", "toys", model, prot, fsep = .Platform$file.sep), pattern = paste(c(prot,"_[a-z]+_", categories,".*"), collapse=""), all.files = FALSE,
                               full.names = FALSE, recursive = FALSE)
    
  }
  boo <- sapply(emfiles, function(ef){
    length(grep('model bic',readLines(ef),  value = TRUE, perl = TRUE)) == 1
  })
  
  emfiles <- emfiles[boo]
  emshortfiles <- emshortfiles[boo]
  
  em_output <- lapply(emfiles, function(file){
    em_file <-  readLines(file)  
    grepped_params <- tail(grep('^\\s+[0-9]', em_file,  value = TRUE, perl = TRUE), categories)
    grepped_lambdas <- tail(grep('^\\[1\\]\\s+\\"expon', em_file,  value = TRUE, perl = TRUE), categories)
    grepped_weights <- tail(grep('^\\[1\\]\\s+[0-9]\\.[0-9]', em_file,  value = TRUE, perl = TRUE), 1)
    
    my_params <-matrix(ncol = 4, nrow = categories, byrow = TRUE)
    colnames(my_params) = c("lambda_weib", "p", "weight", "lambda_exp")
    
    for (c in seq(1, categories, 1)){
      my_params[c,1] <- as.numeric(strsplit(grepped_params[c], '\\s+')[[1]][3])
      my_params[c,2] <- as.numeric(strsplit(grepped_params[c], '\\s+')[[1]][2])
      my_params[c,3] <- as.numeric(strsplit(grepped_weights[1], '\\s+')[[1]][c+1])
      my_str <- strsplit(grepped_lambdas[c], '\\s+')[[1]][4]
      my_params[c,4] <- as.numeric(substr(my_str, 0, nchar(my_str)-1))
    }
    grepped_bic <- tail(em_file, 1)
    bic <- as.numeric(strsplit(grepped_bic, '\\s+')[[1]][2])
    list(params = my_params, bic = bic)
  })
  
  names(em_output) <- emshortfiles
  em_output
}

#' @export
parse_group_output <- function(prot, model, categories){
  emfiles <- list.files(path = file.path(getwd(), "output", "group_likelihood", model, prot, fsep = .Platform$file.sep), pattern = paste(c(prot,"_.*"), collapse=""), all.files = FALSE,
                        full.names = TRUE, recursive = FALSE)
  emshortfiles <- list.files(path = file.path(getwd(), "output", "group_likelihood", model, prot, fsep = .Platform$file.sep), pattern = paste(c(prot,"_.*"), collapse=""), all.files = FALSE,
                             full.names = FALSE, recursive = FALSE)
  boo <- sapply(emfiles, function(ef){
    length(grep('model bic',readLines(ef),  value = TRUE, perl = TRUE)) == 1
  })
  print (emshortfiles)
  emfiles <- emfiles[boo]
  emshortfiles <- emshortfiles[boo]
  
  em_output <- lapply(emfiles, function(file){
    em_file <-  readLines(file)  
    grepped_params <- tail(grep('^\\s+[0-9]', em_file,  value = TRUE, perl = TRUE), categories)
    grepped_lambdas <- tail(grep('^\\[1\\]\\s+\\"expon', em_file,  value = TRUE, perl = TRUE), categories)
    grepped_weights <- tail(grep('^\\[1\\]\\s+[0-9]\\.[0-9]', em_file,  value = TRUE, perl = TRUE), 1)
    
    my_params <-matrix(ncol = 4, nrow = categories, byrow = TRUE)
    colnames(my_params) = c("lambda_weib", "p", "weight", "lambda_exp")
    
    for (c in seq(1, categories, 1)){
      my_params[c,1] <- as.numeric(strsplit(grepped_params[c], '\\s+')[[1]][3])
      my_params[c,2] <- as.numeric(strsplit(grepped_params[c], '\\s+')[[1]][2])
      my_params[c,3] <- as.numeric(strsplit(grepped_weights[1], '\\s+')[[1]][c+1])
      my_str <- strsplit(grepped_lambdas[c], '\\s+')[[1]][4]
      my_params[c,4] <- as.numeric(substr(my_str, 0, nchar(my_str)-1))
    }
    grepped_bic <- tail(em_file, 1)
    bic <- as.numeric(strsplit(grepped_bic, '\\s+')[[1]][2])
    list(params = my_params, bic = bic)
  })
  
  names(em_output) <- emshortfiles
  em_output
}

#' @export
parse_group_LRT <-function(prot){
  #prot <- "h3"
  dir <- file.path(getwd(), "output", "group_likelihood", "weibull", prot, fsep = .Platform$file.sep)
  complfiles <- list.files(path = dir, pattern = paste(c("^", prot,"_.*_complement"), collapse=""), all.files = FALSE,
                        full.names = TRUE, recursive = FALSE)
  groupfiles <- sapply(complfiles, function(f){
    temp <- gregexpr('_', f)
    tochop <- temp[[1]][length(temp[[1]])]
    substr(f,1, tochop-1)
  })

  print (groupfiles)
  print (complfiles)
  
  output <- lapply(groupfiles, function(wgfile){
    gfname <-basename(wgfile)
    
    wgroup_file <-  readLines(wgfile)  
    wg_grepped_lnL <- tail(wgroup_file, 5)[1]
    wg_lnL <- as.numeric(strsplit(wg_grepped_lnL, '\\s+')[[1]][2])
    
    egroup_file <-  readLines(file.path(getwd(), "output", "group_likelihood", "exponential", prot, gfname, fsep = .Platform$file.sep))  
    eg_grepped_lnL <- tail(egroup_file, 5)[1]
    eg_lnL <- as.numeric(strsplit(eg_grepped_lnL, '\\s+')[[1]][2])
    
    lr <- 2*(wg_lnL-eg_lnL)
    
    print (eg_lnL)
    print (wg_lnL)
    print (gfname)
    print (lr)
  })
  

}

color.gradient <- function(x, colors=c("red","yellow","springgreen","royalblue"), colsteps=15) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

#' @export
show_em_results <- function (em_output, model = c("weibull", "exponential"), filter_worst = 0, filter_best = 0){
  bics <- sapply(em_output, function (elm){
    elm$bic
  })
  
  sorted_bics <-sort(bics)
  good_names <- names(head(sorted_bics, length(sorted_bics)-filter_worst))
  bics <- bics[good_names]
  em_output<- em_output[good_names]
  sorted_bics <- sorted_bics[good_names]
  
  bad_names <- names(tail(sorted_bics, length(sorted_bics)-filter_best))
  bics <- bics[bad_names]
  em_output<- em_output[bad_names]
  sorted_bics <- sorted_bics[bad_names]
  
  if (model == "weibull"){
    myplot <- scatterplot3d(0,0,0, color= "white", type="h", xlab = "lambda", ylab = "p", zlab = "weight", xlim = c(0, 0.1), ylim = c(0,10), zlim = c(0,1), pch=19)
    
    for (c in seq(1, categories, 1)){
      firstl <- sapply(em_output, function (e){
        e$params[c,1]
      })
      firstp <- sapply(em_output, function (e){
        e$params[c,2]
      })
      firstw <- sapply(em_output, function (e){
        e$params[c,3]
      })
      
      # only first cluster for all files
      myplot$points3d(firstl, firstp, firstw, col= color.gradient(bics), type="h", pch=19)
    }
    
    cluster <- em_output[grep('clu', names(em_output))]
    
    # only first cluster for all files
    if (length(cluster) > 0){
      for (c in seq(1, categories, 1)){
        myplot$points3d(cluster[[1]]$params[c,1], cluster[[1]]$params[c,2], cluster[[1]]$params[c,3], col= color.gradient(bics), type="h", pch=17)
      }
    }
  }
  
  else {
    myplot <- plot(0,0, col= "white",  xlim = c(0, 0.1), ylim = c(0,1), xlab = "lambda", ylab = "weight", pch=19)
    for (c in seq(1,categories,1)){
      firstl <- sapply(em_output, function (e){
        e$params[c,4]
      })
      firstw <- sapply(em_output, function (e){
        e$params[c,3]
      })
      points(firstl, firstw, col= color.gradient(bics), pch=19)
    }
    
    cluster <- em_output[grep('clu', names(em_output))]
    if (length(cluster) > 0){
    # only first cluster for all files
    for (c in seq(1, categories, 1)){
      points(cluster[[1]]$params[c,4], cluster[[1]]$params[c,3], col= color.gradient(bics), pch=17)
    }
    }
  }
  ##
  
  print (paste ("best bic", names(head(sorted_bics, 1)), head(sorted_bics, 1)))
  print (paste ("worst bic", names(tail(sorted_bics, 1)),tail(sorted_bics, 1)))
  print (paste ("diff", tail(sorted_bics, 1)-head(sorted_bics, 1)))

}

show_em_all_results <- function (em_output, model = c("weibull", "exponential")){
  bics <- sapply(em_output, function (elm){
    elm$bic
  })
  
  sorted_bics <-sort(bics)
  print (length(sorted_bics))
  print (sorted_bics)
  
  if (model == "weibull"){
    myplot <- scatterplot3d(0,0,0, color= "white", type="h", xlab = "lambda", ylab = "p", zlab = "weight", xlim = c(0, 0.1), ylim = c(0,10), zlim = c(0,1), pch=19)
    
    for (c in seq(1, categories, 1)){
      firstl <- sapply(em_output, function (e){
        e$params[c,1]
      })
      firstp <- sapply(em_output, function (e){
        e$params[c,2]
      })
      firstw <- sapply(em_output, function (e){
        e$params[c,3]
      })
      
      # only first cluster for all files
      myplot$points3d(firstl, firstp, firstw, col= color.gradient(bics), type="h", pch=19)
    }
    
    cluster <- em_output[grep('clu', names(em_output))]
    print ("better than cluster")
    print(length(sorted_bics[sorted_bics<cluster[[1]]$bic]))
    print ("equal to cluster")
    print(length(sorted_bics[sorted_bics == cluster[[1]]$bic]) - 1)
    print ("worse than cluster")
    print(length(sorted_bics[sorted_bics > cluster[[1]]$bic]))
    
    
    # only first cluster for all files
    if (length(cluster) > 0){
      for (c in seq(1, categories, 1)){
        myplot$points3d(cluster[[1]]$params[c,1], cluster[[1]]$params[c,2], cluster[[1]]$params[c,3], col= color.gradient(bics), type="h", pch=17)
      }
    }
  }
  
  else {
    myplot <- plot(0,0, col= "white",  xlim = c(0, 0.1), ylim = c(0,1), xlab = "lambda", ylab = "weight", pch=19)
    for (c in seq(1,categories,1)){
      firstl <- sapply(em_output, function (e){
        e$params[c,4]
      })
      firstw <- sapply(em_output, function (e){
        e$params[c,3]
      })
      points(firstl, firstw, col= color.gradient(bics), pch=19)
    }
    
    cluster <- em_output[grep('clu', names(em_output))]
    print ("better than cluster")
    print(length(sorted_bics[sorted_bics<cluster[[1]]$bic]))
    print ("equal to cluster")
    print(length(sorted_bics[sorted_bics == cluster[[1]]$bic])- 1)
    print ("worse than cluster")
    print(length(sorted_bics[sorted_bics > cluster[[1]]$bic]))
    if (length(cluster) > 0){
      # only first cluster for all files
      for (c in seq(1, categories, 1)){
        points(cluster[[1]]$params[c,4], cluster[[1]]$params[c,3], col= color.gradient(bics), pch=17)
      }
    }
  }
  ##
  
  print (paste ("best bic", names(head(sorted_bics, 1)), head(sorted_bics, 1)))
  print (paste ("worst bic", names(tail(sorted_bics, 1)),tail(sorted_bics, 1)))
  print (paste ("diff", tail(sorted_bics, 1)-head(sorted_bics, 1)))
  
}




