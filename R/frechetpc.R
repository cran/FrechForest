#' Ordonne time measurements
#'
#' Used to put the data to the right form at a given time
#'
#'
#' @param X [vector]: matrix of covariables, each variable must be in column
#' @param time [vector]: time measurements of the observations
#' @param id [vector]: IDs of the measurements to their belonging curves
#'
#' @keywords internal
#'
#'
ordonne <- function(X , time , id){
  mat  <- matrix(NA, length(unique(id)), length(unique(time)))
  for( i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    t_w <- time[w]
    w_time <- NULL
    for (j in 1:length(w)){
      w_time <- c(w_time, which(unique(time)==t_w[j]))
    }
    mat[i,w_time] <- X[w]
  }
  return(mat)
}

#' Impurity
#'
#' Computes the Frechet Variance from a group of curves
#'
#' @param traj [vector]: Vector of curves
#' @param time [vector]: Time measurements of the different curves (same length as traj)
#' @param id [vector]: IDs of the different measurements (same length as traj)
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#'
#' @keywords internal
#'
#'
impurity <- function(traj,time, id, timeScale=0.1){
  imp <- 0
  trajLong <- data.frame(id=id,time=time,traj=traj)
  meanF <- kmlShape::meanFrechet(trajLong = trajLong, timeScale = timeScale)
  for (i in unique(id)){
    imp <- imp + kmlShape::distFrechet(meanF$times, meanF$traj, time[which(id==i)], traj[which(id==i)], timeScale = timeScale)
  }
  imp <- imp/length(unique(id))
  return(imp)
}


#' Split Impurity
#'
#' Computes the Frechet variance reduction from a split
#'
#' @param traj [vector]: Trajectories
#' @param time [vector]: Times at which measures are made
#' @param id [vector]: IDs of the different measurements
#' @param split [numeric]: Variable of split
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#'
#' @keywords internal
#'
impurity_split <- function(traj, time, id,split,timeScale=0.1){
  impur <- 0
  imp <- list()
  for (i in unique(split)){
    fils <- unique(id)[which(split==i)]
    prop <- length(fils)/length(unique(id))
    w <- NULL
    for (j in 1:length(fils)){
      w <- c(w, which(id==fils[j]))
    }
    imp[[i]] <- impurity(traj[w], time[w], id[w])
    impur <- impur + imp[[i]]*prop
  }
  return(list(impur=impur, imp_list=imp))
}

#' Title
#'
#' @param traj
#' @param time
#' @param id
#' @param split
#' @param timeScale
#'
#' @keywords internal
#'
ERimpurity_split <- function(traj, time, id,split,timeScale=0.1){
  impur <- 0
  imp <- list()
  for (i in unique(split$sp)){
    prop <- length(unique(id[which(split$sp==i)]))/length(unique(id))
    w <- which(split$sp==i)
    imp[[i]] <- impurity(traj[w], time[w], id[w])
    impur <- impur + imp[[i]]*prop
  }
  return(list(impur=impur, imp_list=imp))
}

#' Variable split
#'
#' Determines the split from a set of input and output variables which are curves:
#'
#' @param X [matrix]: Matrix of explanatory variables, each column codes for a variable.
#' @param Y [vector]: Output curves
#' @param time [vector]: Time at which measures are made
#' @param id [vector]: IDs of trajectories
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the euclidiean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#' @param ... :
#'
#' @keywords internal
#'
#'
var_split <- function(X,Y,time,id,timeScale=0.1,...){
  impur <- rep(0,dim(X)[2])
  toutes_imp <- list()
  split <- list()
  for (i in 1:dim(X)[2]){
    mclds <- kmlShape::cldsWide(ordonne(X[,i], time, id), unique(time), unique(id))
    crit <- kmlShape::kmlShape(mclds, nbClusters = 2, timeScale = timeScale, ...)
    att <- attributes(crit)
    split[[i]] <- att$clusters
    impurete <- impurity_split(Y,time,id,split[[i]], timeScale)
    impur[i] <- impurete$impur
    toutes_imp[[i]] <- impurete$imp_list
  }
  true_split <- which.min(impur)
  split <- split[[true_split]]
  return(list(split=split, impurete=min(impur),impur_list = toutes_imp[[true_split]], variable=which.min(impur)))
}


#' Title
#'
#' @param X
#' @param Y
#' @param time
#' @param id
#' @param timeScale
#' @param ...
#'
#' @keywords internal
#'
ERvar_split <- function(X,Y,time,id,timeScale=0.1,...){
  impur <- rep(0,dim(X)[2])
  toutes_imp <- list()
  split <- list()
  for (i in 1:dim(X)[2]){
    ## Il nous faut tirer les centres ::
    centers <- sample(unique(id),2)
    Split <- rep(NA,dim(X)[2])
    # Il nous faut ensuite calculer toutes les distances afin de dire qui est dans quelle classe::
    for (j in 1:length(unique(id))){
      w <- which(id==unique(id)[j])
      Split[w] <- 1*(kmlShape::distFrechet(time[w],X[w,i],time[which(id==centers[1])], X[which(id==centers[1]),i])<=kmlShape::distFrechet(time[w],X[w,i],time[which(id==centers[2])], X[which(id==centers[2]),i]))+2*(kmlShape::distFrechet(time[w],X[w,i],time[which(id==centers[1])], X[which(id==centers[1]),i])>kmlShape::distFrechet(time[w],X[w,i],time[which(id==centers[2])], X[which(id==centers[2]),i]))
    }
    split[[i]] <- list(sp = Split, centers=centers)

    ### Maintenant il nous faut de quoi calculer l'erreur commise :::

    impurete <- ERimpurity_split(Y,time,id,split[[i]], timeScale)
    impur[i] <- impurete$impur
    toutes_imp[[i]] <- impurete$imp_list
  }
  true_split <- which.min(impur)
  centers <- split[[true_split]]$centers
  split <- split[[true_split]]$sp
  return(list(split=split, impurete=min(impur),impur_list = toutes_imp[[true_split]], variable=which.min(impur), centers=centers))
}



#' Maximal Frechet tree
#'
#' builds the maximal Frechet tree
#'
#' @param X [matrix]: Matrix of explanatory variables, each column codes for a variable
#' @param Y [vector]: Output curves
#' @param id [vector]: IDs of trajectories
#' @param time [vector]: Time at which measures are made
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#' @param ... :
#'
#' @keywords internal
#'
#'
Tmax <- function(X,Y,id,time,timeScale=0.1 , ...){
  id_feuille <- rep(1,length(id))
  id_feuille_prime <- id_feuille
  V_split <- NULL
  impurete <- impurity(Y, time, id, timeScale)
  hist_nodes <- list()
  hist_imp_nodes <- NULL
  imp_nodes <- list()
  imp_nodes[[1]] <- impurity(Y, time, id, timeScale)
  Y_curves <- list()
  Y_curves_imputation <- list()
  hist_imp_nodes <- as.matrix(cbind(1, impurete,length(unique(id))))
  for (p in 1:(length(unique(id))/2-1)){
    count_split <- 0
    for (i in 1:length(unique(id_feuille))){

      w <- which(id_feuille==unique(id_feuille)[i])
      if (length(unique(id[w]))>1){
        feuille_split <- var_split(X[w,],Y[w],time[w],id[w],timeScale = timeScale,...) #### te renvoie l'impurete ainsi que le split ? gauche et a droite

        imp_avant_split <- imp_nodes[[unique(id_feuille)[i]]]
        imp_apres_split <- feuille_split$impurete

        if (imp_apres_split<imp_avant_split){


          Y_curves[[unique(id_feuille)[i]]] <- kmlShape::meanFrechet(as.data.frame(cbind(id[w], time[w], Y[w])))
          #Y_curves[[unique(id_feuille)[i]]] <- DouglasPeuckerNbPoints(Y_curves_imputation[[unique(id_feuille)[i]]]$times, Y_curves_imputation[[unique(id_feuille)[i]]]$traj, nbPoints = length(unique(time)))

          imp_nodes[[2*unique(id_feuille)[i]]] <- feuille_split$impur_list[[1]]
          imp_nodes[[2*unique(id_feuille)[i]+1]] <- feuille_split$impur_list[[2]]

          hist_imp_nodes <- rbind(hist_imp_nodes, c(2*unique(id_feuille)[i],feuille_split$impur_list[[1]], length(which(feuille_split$split==1))))
          hist_imp_nodes <- rbind(hist_imp_nodes, c(2*unique(id_feuille)[i]+1,feuille_split$impur_list[[2]], length(which(feuille_split$split==2))))
          gauche_id <- unique(id[w])[which(feuille_split$split==1)]
          droit_id <- unique(id[w])[which(feuille_split$split==2)]
          V_split <- rbind(V_split,c(unique(id_feuille)[i],feuille_split$variable))
          w_gauche <- NULL
          w_droit <- NULL
          print(paste("Split on the variable", feuille_split$variable))
          for (k in 1:length(gauche_id)){
            w_gauche <- c(w_gauche, which(id==gauche_id[k]))
          }
          for (k in 1:length(droit_id)){
            w_droit <- c(w_droit, which(id==droit_id[k]))
          }
          id_feuille_prime[w_gauche] <- 2*(unique(id_feuille)[i])
          id_feuille_prime[w_droit] <- 2*(unique(id_feuille)[i])+1
          trajG <- as.data.frame(cbind(id[w_gauche], time[w_gauche], X[w_gauche,feuille_split$variable]))
          trajD <- as.data.frame(cbind(id[w_droit], time[w_droit], X[w_droit,feuille_split$variable]))
          meanFg <- as.matrix(kmlShape::meanFrechet(trajG))
          meanFd <- as.matrix(kmlShape::meanFrechet(trajD))
          hist_nodes[[2*(unique(id_feuille)[i])]] <- meanFg
          hist_nodes[[2*(unique(id_feuille)[i])+1]] <- meanFd
          count_split <- count_split+1
        }
      }
    }

    id_feuille <- id_feuille_prime
    imp <- 0
    imp1 <- NULL
    for (i in 1:length(unique(id_feuille))){
      w <- which(id_feuille==unique(id_feuille)[i])
      prop <- length(unique(id[w]))/length(unique(id))
      imp <- imp+ impurity(Y[w], time[w], id[w])*prop
    }
    impurete <- c(impurete,imp)

    if (count_split ==0 ){

      impurete <- c(impurete,imp)
      V_split <- data.frame(V_split)
      names(V_split) <- c("num_noeud", "var_split")
      for (q in unique(id_feuille)){
        w <- which(id_feuille == q)
        Y_curves[[q]] <- kmlShape::meanFrechet(data.frame(id[w], time[w], Y[w]))
      }
      return(list(feuilles = id_feuille, id=id, V_split=V_split, impurity=impurete, hist_nodes=hist_nodes, Y_curves = Y_curves, time = time, Y=Y, hist_imp_nodes=hist_imp_nodes, Alpha =0))
    }
  }

  V_split <- data.frame(V_split)
  names(V_split) <- c("num_noeud", "var_split")
  for (q in unique(id_feuille)){
    w <- which(id_feuille == q)
    Y_curves[[q]] <- kmlShape::meanFrechet(data.frame(id[w], time[w], Y[w]))
    #Y_curves[[q]] <- DouglasPeuckerNbPoints(Y_curves_imputation[[q]]$times, Y_curves_imputation[[q]]$traj, nbPoints = length(unique(time)))
  }
  return(list(feuilles = id_feuille, id=id, V_split=V_split, impuity=impurete, hist_nodes=hist_nodes, Y_curves= Y_curves, time=time, Y=Y, hist_imp_nodes=hist_imp_nodes, Alpha=0))
}


#' Frechet Tree prediction
#'
#' Given a Frechet tree and new input trajectories predictors \code{X}, this function returns the identifier of the leaf in which each observation falls.
#'
#' @param tree : Frechet tree obtained with the function \code{\link{FrechTree}}.
#' @param X [matrix]: a data frame or a matrix of trajectories predictors.
#' @param time [vector]: time measurements of the new trajectories to predict.
#' @param id [vector]: identifier, one for each trajectory to attribute each measurement of \code{X} to one of the trajectories.
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#'
#' @return a vector of the identifier of the leaf in which each observation falls.
#' @keywords internal
#'
pred.FT <- function(tree, X, time, id, timeScale=0.1){
  pred <- rep(NA,length(unique(id)))
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    noeud_courant <- 1
    while (is.element(noeud_courant, tree$feuilles)==FALSE){
      var_split <- tree$V_split[which(tree$V_split[,1]==noeud_courant),2]
      meanG <- tree$hist_nodes[[2*noeud_courant]]
      meanD <- tree$hist_nodes[[2*noeud_courant+1]]
      distG <- kmlShape::distFrechet(meanG[,1], meanG[,2], time[w], X[w,var_split], timeScale = timeScale)
      distD <- kmlShape::distFrechet(meanD[,1], meanD[,2], time[w], X[w,var_split], timeScale = timeScale)
      if (distG < distD) { noeud_courant <- 2*noeud_courant}
      if (distD <= distG) {noeud_courant <- 2*noeud_courant +1}
    }
    pred[i] <- noeud_courant
  }
  return(pred)
}



#' Frechet Tree prediction
#'
#' Given a Frechet tree and new input trajectories predictors \code{X}, this function returns the identifier of the leaf in which each observation falls.
#'
#' @param object : Frechet tree obtained with the function \code{FrechTree}.
#' @param X [matrix]: a data frame or a matrix of trajectories predictors.
#' @param time [vector]: time measurements of the new trajectories to predict.
#' @param id [vector]: identifier, one for each trajectory to attribute each measurement of \code{X} to one of the trajectories.
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#' @param ... : optional parameters to be passed to the low level function.
#'
#' @return a matrix of  the identifier of the leaf (second column) in which each individual (first column) falls.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(10)
#' data <- DataGenCurves(40)
#' Ft <- FrechTree(data$X,data$Y, data$id,data$time,select = "Hubert", toPlot="none")
#' }
predict.FrechTree <- function(object, X, time, id, timeScale=0.1, ...){
  tree <- object
  pred <- matrix(NA,length(unique(id)),2)
  colnames(pred) <- c("ID","leaf")
  for (i in 1:length(unique(id))){
    w <- which(id==unique(id)[i])
    noeud_courant <- 1
    while (is.element(noeud_courant, tree$feuilles)==FALSE){
      var_split <- tree$V_split[which(tree$V_split[,1]==noeud_courant),2]
      meanG <- tree$hist_nodes[[2*noeud_courant]]
      meanD <- tree$hist_nodes[[2*noeud_courant+1]]
      distG <- kmlShape::distFrechet(meanG[,1], meanG[,2], time[w], X[w,var_split], timeScale = timeScale)
      distD <- kmlShape::distFrechet(meanD[,1], meanD[,2], time[w], X[w,var_split], timeScale = timeScale)
      if (distG < distD) { noeud_courant <- 2*noeud_courant}
      if (distD <= distG) {noeud_courant <- 2*noeud_courant +1}
    }
    pred[i,] <- c(unique(id)[i],noeud_courant)
  }
  return(pred)
}

#' Branche
#'
#' Used to prune a maximal Frechet tree
#'
#' @param tree : Frechet tree
#' @param t [numeric]: a node of tree
#'
#' @keywords internal
#'
#'
branche <- function(tree, t){
  Y_curves <- list()
  f <- unique(tree$feuilles)
  sous_split <- tree$V_split[which(tree$V_split[,1]==t),]
  N <- 2
  g <- which(tree$V_split[,1]==2*t)
  d <- which(tree$V_split[,1]==2*t+1)
  noeuds_courants <- tree$V_split[c(g,d),1]
  noeuds_courants1 <- noeuds_courants
  sous_split <- rbind(sous_split, tree$V_split[c(g,d),])
  sous_feuilles <- NULL
  hist_nodes <- list()
  if (length(g)>0) {hist_nodes[[2*t]] <- tree$hist_nodes[[2*t]]}
  if (length(d)>0) {hist_nodes[[2*t+1]] <- tree$hist_nodes[[2*t+1]]}
  if (length(d)== 0) {sous_feuilles <- c(sous_feuilles, 2*t+1)
  Y_curves[[2*t+1]] <- tree$Y_curves[[2*t+1]]}
  if (length(g)== 0) {sous_feuilles <- c(sous_feuilles, 2*t)
  Y_curves[[2*t]] <- tree$Y_curves[[2*t]]}
  racine <- t
  if (length(noeuds_courants)>0) {
    while(N>0){
      p <- 0
      courant_prime <- NULL
      for (t in noeuds_courants){
        g <- which(tree$V_split[,1]==2*t)
        d <- which(tree$V_split[,1]==2*t+1)

        if (length(g)>0){ p <- p+2
        courant_prime <- c(courant_prime, tree$V_split[g,1])
        sous_split <- rbind(sous_split, tree$V_split[g,])
        hist_nodes[[2*t]] <- tree$hist_nodes[[2*t]]}

        if (length(d)>0){ p <- p+2
        courant_prime <- c(courant_prime, tree$V_split[d,1])
        sous_split <- rbind(sous_split, tree$V_split[d,])
        hist_nodes[[2*p+1]] <- tree$hist_nodes[[2*t+1]]}

        if(length(g)==0) {sous_feuilles <- c(sous_feuilles,2*t)
        Y_curves[[2*t]] <- tree$Y_curves[[2*t]]}

        if (length(d)==0) { sous_feuilles <- c(sous_feuilles, 2*t+1)
        Y_curves[[2*t+1]] <- tree$Y_curves[[2*t+1]]}
      }
      noeuds_courants <- courant_prime
      N <-p
    }
  }

  if (length(noeuds_courants1)==0) {sous_feuilles <- c(2*t, 2*t+1)}
  s_feuilles <- NULL
  s_id <- NULL
  s_time <- NULL
  s_Y <- NULL
  for(f in unique(sous_feuilles)){
    w <- which(tree$feuilles==f)
    s_feuilles <- c(s_feuilles, tree$feuilles[w])
    s_id <- c(s_id, tree$id[w])
    s_time <- c(s_time, tree$time[w])
    s_Y <- c(s_Y, tree$Y[w])
  }
  #### il faut maintenant calculer l'impurete de la branche ainsi que celle du noeud t
  #### impurete dans le noeud racine :::
  impurity_racine <- tree$hist_imp_nodes[which(tree$hist_imp_nodes[,1]==racine),2]
  n_racine <- tree$hist_imp_nodes[which(tree$hist_imp_nodes[,1]==racine),3]
  n_base <- tree$hist_imp_nodes[1,3]
  impurity_racine <- impurity_racine*(n_racine/n_base)

  impurity_T <- 0
  for (i in unique(s_feuilles)){
    w <- which(tree$hist_imp_nodes[,1]==i)
    prop <- tree$hist_imp_nodes[w,3]/n_base
    impurity_T <- impurity_T + tree$hist_imp_nodes[w,2]*prop
  }

  return(list(feuilles=s_feuilles, V_split = sous_split, id= s_id, time=s_time,hist_nodes=hist_nodes, Y=s_Y, impurity_T = impurity_T, impurity_racine = impurity_racine, n_racine=n_racine, Y_curves=tree$Y_curves))
}

#' Noeud deg
#'
#' This function is used to detect the nodes to prune at each step
#'
#' @param tree : Frechet tree
#'
#' @keywords internal
#'
#'
noeuds_deg <- function(tree){
  noeuds <- tree$V_split[,1]
  deg <- NULL
  alpha <- rep()
  mat_pen <- matrix(0, length(noeuds), 5)
  mat_pen[,1] <- noeuds
  for (t in noeuds){
    b <- branche(tree,t) ### on recupère la branche associee à t
    if (length(unique(b$feuilles))>1){
      mat_pen[which(noeuds==t), 2] <- b$impurity_racine
      mat_pen[which(noeuds==t), 3] <- b$impurity_T
      mat_pen[which(noeuds==t), 4] <- length(unique(b$feuilles))
      mat_pen[which(noeuds==t), 5] <- (b$impurity_racine-b$impurity_T)/(length(unique(b$feuilles))-1)}
    #pen <- mat_pen[which(noeuds==t), 5]
    #err <- b$impurity_T + pen*length(unique(b$feuilles)) - b$impurity_racine - pen
    #print(err)
  }
  alpha <- min(mat_pen[,5])
  err <- rep(0, length(noeuds))
  for (i in  1:dim(mat_pen)[1]){
    err[i] <- round(mat_pen[i,3] + alpha*mat_pen[i,4] - mat_pen[i,2] - alpha, 5)
    if (err[i]==0){
      deg <- rbind(deg, c(mat_pen[i,1], alpha))
    }
  }
  return(deg)
}

#' Pruning Frechet tree
#'
#' To prune a Frechet tree for curves
#'
#' @param tree : Frechet tree
#'
#' @keywords internal
#'
#'
elagage <- function(tree){

  t_feuilles <- NULL
  t_id <- NULL
  t_time <- NULL
  t_split <- NULL
  t_hist <- NULL
  t_Y <- tree$Y


  tree_courant <- tree
  nb_feuilles <- length(unique(tree$feuilles))
  n_max <- nb_feuilles
  courant <- 2
  TREES <- list()
  TREES[[1]] <- tree
  ##### il faut aussi trouver les d?coupe superficielles :::: on garde un historique des d?coupes :::::
  while(nb_feuilles >1){
    deg <- noeuds_deg(tree_courant)
    if (dim(deg)[1]>1) deg <- apply(deg, 2, sort, decreasing=TRUE)
    t_feuilles_courant <- tree_courant$feuilles
    for (t in deg[,1]){
      b <- branche(tree_courant, t)
      feuilles_b <- unique(b$feuilles)
      w_feuilles <- NULL
      for (f in feuilles_b ){
        w_feuilles <- c(w_feuilles, which(tree_courant$feuilles==f))
      }
      t_feuilles_courant[w_feuilles] <- t
      #### il faut maintenant retirer toute la branche de t

      nodes <- b$V_split[,1]
      w_nodes <- NULL
      for (node in nodes){
        w_nodes <- c(w_nodes, which(tree_courant$V_split[,1]==node))
      }

      t_split_courant <- as.matrix(tree_courant$V_split[-w_nodes,, drop = FALSE])
      ##### il faut alors recalculer l'importance dans les feuilles
      tree_courant <- list(feuilles=t_feuilles_courant, V_split = t_split_courant, id= tree$id, time=tree$time,hist_nodes=tree$hist_nodes, Y=tree$Y, hist_imp_nodes=tree$hist_imp_nodes, Y_curves = tree$Y_curves, Alpha=unique(deg[,2]))
      TREES[[courant]] <- tree_courant
    }
    courant <- courant+1
    nb_feuilles <- length(unique(tree_courant$feuilles))
  }

  return(TREES)
}


#' Cross validated Frechet tree
#'
#' Builds a Frechet tree that is selected with cross validation
#'
#' @param X [matrix]: a data frame or a matrix of trajectories predictors. Each colunm codes for a trajectory predictor.
#' @param Y [vector]: a vector containing the output trajectories (same length as \code{nrow(X)}).
#' @param id [vector (factor)]: identifier, one for each trajectory to attribute each measurement of \code{X} and \code{Y} to one of the trajectories (same length as \code{Y}).
#' @param time [vector]: time measurements for the observations of both \code{X} and \code{Y} (same length as \code{Y}).
#' @param ... : optional parameters to be passed to the low level function
#'
#' @return
#' @keywords internal
#'
TreeShapeCV <- function(X,Y,id,time, ...){

  TMAX <- Tmax(X,Y,id,time, ...)

  elag_max <- elagage(TMAX)
  ALPHA <- rep(NA, length(elag_max))
  for (i in 1:length(ALPHA)){
    ALPHA[i] <- elag_max[[i]]$Alpha
  }
  #### on transforme le tout en beta
  beta <- rep(NA, length(ALPHA))
  beta[length(ALPHA)] <- ALPHA[length(ALPHA)]
  for (i in 1:(length(ALPHA)-1)){
    beta[i] <- sqrt(abs(ALPHA[i]*ALPHA[i+1]))
  }

  #### Il faut faire faire les sous ensemble de validation crois?e::::
  ELAG <- list()
  n_folds <- 10
  VC <- sample(rep(1:n_folds, length.out= length(unique(id))))
  tmax <- list()
  APP <- list()
  err <- matrix(0, length(beta), n_folds)
  for (p in 1:n_folds){
    app <- unique(id)[which(VC!=p)] ### on r?cup?re les identifiants
    w <- NULL
    for (a in app){
      w <- c(w, which(id==a))
    }
    APP[[p]] <- w
    X.app <- X[w,, drop=FALSE]
    Y.app <- Y[w]
    time.app <- time[w]
    id.app <- id[w]

    X.val <- X[-w,, drop=FALSE]
    Y.val <- Y[-w]
    time.val <- time[-w]
    id.val <- id[-w]

    tmax[[p]] <- Tmax(X.app,Y.app,id.app, time.app, ...)

    ELAG[[p]] <- elagage(tmax[[p]])

    pen <- rep(NA,length(ELAG[[p]]))
    #pen[length(ELAG[[p]])] <- ELAG[[p]][[length(ELAG[[p]])]]$Alpha
    for (l in 1:length(pen)){
      pen[l] <- ELAG[[p]][[l]]$Alpha
    }

    for (k in 1:length(beta)){
      sous_arbre <- ELAG[[p]][[which.min(abs(pen-beta[k]))]]
      where <- pred.FT(sous_arbre,X.val,time.val , id.val ) #### on doit trouver les feuilles de pr?diction :::
      ##### il nous faut maintenant pr?dire les diff?rentes courbes ::::
      err_courant <- rep(0, length(where))
      for (j in 1:length(where)){
        ww <- which(id.val == unique(id.val)[j])
        #mean_courant <- DouglasPeuckerEpsilon(sous_arbre$Y_curves[[where[j]]][,1],sous_arbre$Y_curves[[where[j]]][,2], 0.01)
        err_courant[j] <-  kmlShape::distFrechet(time.val[ww], Y.val[ww],sous_arbre$Y_curves[[where[j]]][,1], sous_arbre$Y_curves[[where[j]]][,2])
      }
      err[k,p] <- mean(err_courant)
    }
  }

  SD <- apply(err, 1, "sd")
  err_M <- apply(err, 1, "mean")

  ### On prend le meilleur mod?le
  seuil <- min(err_M) + SD[which.min(err_M)]
  ### il faut s?lectionner le meilleur arbre
  optimal.tree <- max(which(err_M<=seuil))
  beta.opt <- beta[optimal.tree]
  final_tree <- elag_max[[max(which(err_M<=seuil))]]

  #### On va s?lectionner l'arbre optimal pour chaque ensemble d'apprentissage puis calculer l'importance des variables sur ceux-ci::
  Importance <- matrix(0, n_folds, dim(X)[2])
  err_arbres_select <- rep(NA, n_folds)
  for (k in 1:n_folds){
    ### on r?cup?re les ?l?ments de validation::::
    X.val <- X[-APP[[k]],]
    Y.val <- Y[-APP[[k]]]
    time.val <- time[-APP[[k]]]
    id.val <- id[-APP[[k]]]

    pen <- rep(NA,length(ELAG[[k]]))
    #pen[length(ELAG[[p]])] <- ELAG[[p]][[length(ELAG[[p]])]]$Alpha
    for (l in 1:length(pen)){
      pen[l] <- ELAG[[k]][[l]]$Alpha
    }
  }

  m_leafs <- max(unique(final_tree$feuilles))
  return(list(feuilles = final_tree$feuilles, V_split=final_tree$V_split, hist_nodes=final_tree$hist_nodes[1:m_leafs], Y_curves=final_tree$Y_curves[1:m_leafs], err_elag = err_M, seuil=seuil))
}


#' Hubert's statistic for Frechet tree
#'
#' Computes the Huberts statistic for a given Frechet tree
#'
#' @param tree : Frechet tree
#'
#' @keywords internal
#'
#'
Huberts <- function(tree){
  hub_stat <- 0
  n <- length(unique(tree$id))
  M <- (n*(n-1))/2
  for (i in 1:(length(unique(tree$id))-1)){
    for (j in (i+1):length(unique(tree$id))){
      w_i <- which(tree$id==unique(tree$id)[i])
      w_j <- which(tree$id==unique(tree$id)[j])
      leaf_i <- tree$Y_curves[[unique(tree$feuilles[w_i])]]
      leaf_j <- tree$Y_curves[[unique(tree$feuilles[w_j])]]
      U <- kmlShape::distFrechet(tree$time[w_i], tree$Y[w_i], tree$time[w_j], tree$Y[w_j])
      V <- kmlShape::distFrechet(leaf_i$times,leaf_i$traj, leaf_j$times, leaf_j$traj)
      hub_stat <- hub_stat+ (U*V)/M
    }
  }
  hub_stat <- hub_stat
  return(hub_stat)
}


#' Pruning function for Hubert Frechet tree
#'
#' To prune a tree and selects the optimal sub Frechet tree according to the Huberts statistic
#'
#' @param tree : Frechet tree
#'
#' @keywords internal
#'
#'
elagageHuberts <- function(tree){

  tree_courant <- tree
  nb_feuilles <- length(unique(tree$feuilles))
  n_max <- nb_feuilles
  courant <- 2
  TREES <- list()
  TREES[[1]] <- tree ### on stocke la premi?re partition tmax

  huberts_stat <- Huberts(tree_courant)
  while(nb_feuilles >1){
    deg <- noeuds_deg(tree_courant)
    if (dim(deg)[1]>1) deg <- apply(deg, 2, sort, decreasing=TRUE)
    t_feuilles_courant <- tree_courant$feuilles
    for (t in deg[,1]){
      b <- branche(tree_courant, t)
      feuilles_b <- unique(b$feuilles)
      w_feuilles <- NULL
      for (f in feuilles_b ){
        w_feuilles <- c(w_feuilles, which(tree_courant$feuilles==f))
      }
      t_feuilles_courant[w_feuilles] <- t
      #### il faut maintenant retirer toute la branche de t

      nodes <- b$V_split[,1]
      w_nodes <- NULL
      for (node in nodes){
        w_nodes <- c(w_nodes, which(tree_courant$V_split[,1]==node))
      }

      t_split_courant <- as.matrix(tree_courant$V_split[-w_nodes,, drop = FALSE])
      ##### il faut alors recalculer l'importance dans les feuilles
      tree_courant <- list(feuilles=t_feuilles_courant, V_split = t_split_courant, id= tree$id, time=tree$time,hist_nodes=tree$hist_nodes, Y=tree$Y, hist_imp_nodes=tree$hist_imp_nodes, Y_curves = tree$Y_curves, Alpha=unique(deg[,2]))
      TREES[[courant]] <- tree_courant

      #### calcul de la statistique de Huberts
    }

    huberts_stat <- c(huberts_stat, Huberts(tree_courant))
    courant <- courant+1
    nb_feuilles <- length(unique(tree_courant$feuilles))
  }

  return(list(TREES=TREES, huberts_stat=huberts_stat))
}


#' FrechTree
#'
#' This function runs Frechet tree for longitudinal data using some shape respecting distance and mean.
#'
#' @inheritParams TreeShapeCV
#' @param select [character]: a character string indicating which criteria to select the optimal pruned sub-tree, either "CV" for cross validation on the prediction error or "Hubert" for the Hubert's statistics. Default is "CV".
#'
#' @return a Frechet treee which is a list of the following elements :\itemize{
#' \item \code{feuilles}: a vector indicating in which leaf is each measurement.
#' \item \code{V_split}: a matrix of two columns that describes the structure of the tree. Each row codes for a split, the first column indicates the node and the second column gives the associated splitting variable.
#' \item \code{Y_curves}: a list of the predicted trajectories (Frechet mean) for each leaf of the optimal Frechet tree.
#' \item \code{hist_nodes}: a list of the representative trajectories for each node according to the associated splitting variable.
#' }
#' @export
#'
#'
#' @examples
#' \dontrun{
#' set.seed(10)
#' data <- DataGenCurves(40)
#' Ft <- FrechTree(data$X,data$Y, data$id,data$time,select = "Hubert", toPlot="none")
#'}
FrechTree <- function(X,Y,id,time,select="CV",...){

  if (select=="CV"){
    treeshape <- TreeShapeCV(X,Y,id,time,...)
    class(treeshape) <- c("FrechTree")
    return(treeshape)
  }

  if (select=="Hubert"){
    tmax <- Tmax(X,Y,id,time, ...)
    elag <- elagageHuberts(tmax)
    tree.opt <- elag$TREES[[which.max(elag$huberts_stat)]]
    m_leafs <- max(unique(tree.opt$feuilles))
    treeshape <- list(feuilles = tree.opt$feuilles, V_split=tree.opt$V_split, hist_nodes=tree.opt$hist_nodes[1:m_leafs], Y_curves=tree.opt$Y_curves[1:m_leafs],huberts=elag$huberts_stat)
    class(treeshape) <- c("FrechTree")
    return(treeshape)
  }

}


#' Permutation
#'
#' This function is used to permute curves in order to compute the OOB error from a Frechet random forest
#'
#' @param Courbes [vector]: Curves, ordered by time measurements
#' @param id : IDs of measurements
#'
#' @keywords internal
#'
#'
permutation <- function(Courbes,id){
  perm <- sample(unique(id), length(unique(id))) #### on change les identifiants ::
  new <- NULL
  for (i in perm){
    w <- which(id==i)
    new <- c(new,Courbes[w])
  }
  return(new)
}



#'
#' @param X [matrix]: Matrix of explanatory variables, each column codes for a variable
#' @param Y [vector]: Output curves
#' @param id [vector]: IDs of measurements
#' @param time [vector]: time measurements
#' @param mtry [numeric]: number of variables randomly chosen at each split
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#' @param ... :
#'
#' @keywords internal
#'
ERtmax <- function(X,Y,id, time, mtry, timeScale=0.1, ...){
  impurity_feuilles <- NULL
  V_split <- NULL
  hist_nodes <- list()
  id_boot <- unique(sample(unique(id), length(unique(id)), replace=TRUE)) #### on tire l'?chantillon bootstrap
  boot <- id_boot
  W <- NULL
  for (k in id_boot){
    W <- c(W, which(id==k))
  }
  X_boot <- X[W,, drop=FALSE] ### bootstrap pour X
  Y_boot <- Y[W] ### idem pour Y
  time_boot <- time[W]
  id_boot <- id[W]
  impurete <- impurity(Y_boot, time_boot, id_boot, timeScale=timeScale) #### impurete dans l'ech boot au depart

  ##### creation d'une matrice de stockage
  id_feuille <- rep(1,length(id_boot)) #### localisation des feuilles de l'arbre
  id_feuille_prime <- id_feuille
  hist_imp_nodes <- NULL

  for (p in 1:(length(unique(id_boot))/2-1)){
    count_split <- 0
    for (i in 1:length(unique(id_feuille))){
      variables <- sample(c(1:dim(X[,,drop=FALSE])[2]),mtry)
      w <- which(id_feuille==unique(id_feuille)[i])
      if (length(unique(id_boot[w]))>1){

        if (length(which(hist_imp_nodes[,1]==unique(id_feuille)[i]))==0){
          hist_imp_nodes <- rbind(hist_imp_nodes, c(unique(id_feuille)[i],impurity(Y_boot[w],time_boot[w],id_boot[w],timeScale=timeScale), length(unique(id_boot[w]))))
        }

        feuille_split <- ERvar_split(X_boot[w,variables, drop=FALSE],Y_boot[w],time_boot[w],id_boot[w], timeScale=timeScale)

        gauche_id <- id_boot[w][which(feuille_split$split==1)]
        droit_id <- id_boot[w][which(feuille_split$split==2)]

        V_split <- rbind(V_split,c(unique(id_feuille)[i],variables[feuille_split$variable]))

        w_gauche <- NULL
        w_droit <- NULL
        #print(paste("d?coupage sur la variable",variables[feuille_split$variable]))
        for (k in 1:length(gauche_id)){
          w_gauche <- c(w_gauche, which(id_boot==gauche_id[k]))
        }
        for (k in 1:length(droit_id)){
          w_droit <- c(w_droit, which(id_boot==droit_id[k]))
        }

        id_feuille_prime[w_gauche] <- 2*(unique(id_feuille)[i])
        id_feuille_prime[w_droit] <- 2*(unique(id_feuille)[i])+1

        trajG <- as.data.frame(cbind(id_boot[w_gauche], time_boot[w_gauche], X_boot[w_gauche,variables[feuille_split$variable], drop=FALSE]))
        trajD <- as.data.frame(cbind(id_boot[w_droit], time_boot[w_droit], X_boot[w_droit,variables[feuille_split$variable], drop=FALSE]))

        meanFg <- data.frame(time_boot[which(id_boot==feuille_split$centers[1])],X_boot[which(id_boot==feuille_split$centers[1]),variables[feuille_split$variable]])
        meanFd <- data.frame(time_boot[which(id_boot==feuille_split$centers[2])],X_boot[which(id_boot==feuille_split$centers[2]),variables[feuille_split$variable]])

        names(meanFg) <- c("times", "traj")
        names(meanFd) <- c("times", "traj")

        hist_nodes[[2*(unique(id_feuille)[i])]] <- meanFg
        hist_nodes[[2*(unique(id_feuille)[i])+1]] <- meanFd
        count_split <- count_split+1
      }
    }

    id_feuille <- id_feuille_prime
    imp <- 0
    for (i in 1:length(unique(id_feuille))){
      w <- which(id_feuille==unique(id_feuille)[i])
      prop <- length(unique(id_boot[w]))/length(unique(id_boot))
      imp <- imp+ impurity(Y_boot[w], time_boot[w], id_boot[w] ,timeScale = timeScale)*prop
    }
    impurete <- c(impurete,imp)

    if (count_split ==0 ){
      imp <- 0
      Y_curves <- list()
      for (q in unique(id_feuille)){
        w <- which(id_feuille==q)
        #imp1 <- impurity(Y[w], time[w], id[w])
        Y_curves[[q]] <- kmlShape::meanFrechet(data.frame(id_boot[w], time_boot[w], Y_boot[w]), timeScale = timeScale)
        #impurity_feuilles <- rbind(impurity_feuilles, c(q, imp1, length(unique(id[w]))))
      }
      V_split <- data.frame(V_split)
      names(V_split) <- c("num_noeud", "var_split")
      return(list(feuilles = id_feuille, id=id_boot, V_split=V_split, impurity=impurete, hist_nodes=hist_nodes, Y_curves = Y_curves, time=time_boot , boot=boot, hist_imp_nodes=hist_imp_nodes))
    }
  }
  V_split <- data.frame(V_split)
  names(V_split) <- c("num_noeud", "var_split")
  Y_curves <- list()
  for (q in unique(id_feuille)){
    w <- which(id_feuille==q)
    #imp1 <- impurity(Y[w], time[w], id[w])
    Y_curves[[q]] <- kmlShape::meanFrechet(data.frame(id_boot[w], time_boot[w], Y_boot[w]), timeScale = timeScale)
    #impurity_feuilles <- rbind(impurity_feuilles, c(q, imp1, length(unique(id[w]))))
  }
  return(list(feuilles = id_feuille, id=id_boot, V_split=V_split, impuity=impurete, hist_nodes=hist_nodes, Y_curves= Y_curves, time=time_boot , boot=boot, hist_imp_nodes=hist_imp_nodes))
}


#' Maximal randomized Fréchet tree
#'
#' Build a maximal randomized Fréchet tree
#'
#' @param X [matrix]: Matrix of explanatory variables, each column codes for a variable
#' @param Y [vector]: Output curves
#' @param id [vector]: IDs of measurements
#' @param time [vector]: time measurements
#' @param mtry [numeric]: number of variables randomly chosen at each split
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#' @param ... :
#'
#' @keywords internal
#'
Rtmax <- function(X,Y,id, time, mtry, timeScale=0.1, ...){
  impurity_feuilles <- NULL
  V_split <- NULL
  hist_nodes <- list()
  id_boot <- unique(sample(unique(id), length(unique(id)), replace=TRUE)) #### on tire l'?chantillon bootstrap
  boot <- id_boot
  W <- NULL
  for (k in id_boot){
    W <- c(W, which(id==k))
  }
  X_boot <- X[W,, drop=FALSE] ### bootstrap pour X
  Y_boot <- Y[W] ### idem pour Y
  time_boot <- time[W]
  id_boot <- id[W]
  impurete <- impurity(Y_boot, time_boot, id_boot, timeScale=timeScale) #### impureté dans l'ech boot au départ

  ##### création d'une matrice de stockage
  id_feuille <- rep(1,length(id_boot)) #### localisation des feuilles de l'arbre
  id_feuille_prime <- id_feuille
  hist_imp_nodes <- NULL

  for (p in 1:(length(unique(id_boot))/2-1)){
    count_split <- 0
    for (i in 1:length(unique(id_feuille))){
      variables <- sample(c(1:dim(X[,,drop=FALSE])[2]),mtry)
      w <- which(id_feuille==unique(id_feuille)[i])
      if (length(unique(id_boot[w]))>1){

        if (length(which(hist_imp_nodes[,1]==unique(id_feuille)[i]))==0){
          hist_imp_nodes <- rbind(hist_imp_nodes, c(unique(id_feuille)[i],impurity(Y_boot[w],time_boot[w],id_boot[w],timeScale=timeScale), length(unique(id_boot[w]))))
        }

        feuille_split <- var_split(X_boot[w,variables, drop=FALSE],Y_boot[w],time_boot[w],id_boot[w], timeScale=timeScale,...)
        gauche_id <- unique(id_boot[w])[which(feuille_split$split==1)]
        droit_id <- unique(id_boot[w])[which(feuille_split$split==2)]
        V_split <- rbind(V_split,c(unique(id_feuille)[i],variables[feuille_split$variable]))
        w_gauche <- NULL
        w_droit <- NULL
        #print(paste("d?coupage sur la variable",variables[feuille_split$variable]))
        for (k in 1:length(gauche_id)){
          w_gauche <- c(w_gauche, which(id_boot==gauche_id[k]))
        }
        for (k in 1:length(droit_id)){
          w_droit <- c(w_droit, which(id_boot==droit_id[k]))
        }
        id_feuille_prime[w_gauche] <- 2*(unique(id_feuille)[i])
        id_feuille_prime[w_droit] <- 2*(unique(id_feuille)[i])+1
        trajG <- as.data.frame(cbind(id_boot[w_gauche], time_boot[w_gauche], X_boot[w_gauche,variables[feuille_split$variable], drop=FALSE]))
        trajD <- as.data.frame(cbind(id_boot[w_droit], time_boot[w_droit], X_boot[w_droit,variables[feuille_split$variable], drop=FALSE]))
        meanFg <- as.matrix(kmlShape::meanFrechet(trajG, timeScale = timeScale))
        meanFd <- as.matrix(kmlShape::meanFrechet(trajD, timeScale = timeScale))
        hist_nodes[[2*(unique(id_feuille)[i])]] <- meanFg
        hist_nodes[[2*(unique(id_feuille)[i])+1]] <- meanFd
        count_split <- count_split+1
      }
    }

    id_feuille <- id_feuille_prime
    imp <- 0
    for (i in 1:length(unique(id_feuille))){
      w <- which(id_feuille==unique(id_feuille)[i])
      prop <- length(unique(id_boot[w]))/length(unique(id_boot))
      imp <- imp+ impurity(Y_boot[w], time_boot[w], id_boot[w] ,timeScale = timeScale)*prop
    }
    impurete <- c(impurete,imp)

    if (count_split ==0 ){
      imp <- 0
      Y_curves <- list()
      for (q in unique(id_feuille)){
        w <- which(id_feuille==q)
        #imp1 <- impurity(Y[w], time[w], id[w])
        Y_curves[[q]] <- kmlShape::meanFrechet(data.frame(id_boot[w], time_boot[w], Y_boot[w]), timeScale = timeScale)
        #impurity_feuilles <- rbind(impurity_feuilles, c(q, imp1, length(unique(id[w]))))
      }
      V_split <- data.frame(V_split)
      names(V_split) <- c("num_noeud", "var_split")
      return(list(feuilles = id_feuille, id=id_boot, V_split=V_split, impurity=impurete, hist_nodes=hist_nodes, Y_curves = Y_curves, time=time_boot , boot=boot, hist_imp_nodes=hist_imp_nodes))
    }
  }
  V_split <- data.frame(V_split)
  names(V_split) <- c("num_noeud", "var_split")
  Y_curves <- list()
  for (q in unique(id_feuille)){
    w <- which(id_feuille==q)
    #imp1 <- impurity(Y[w], time[w], id[w])
    Y_curves[[q]] <- kmlShape::meanFrechet(data.frame(id_boot[w], time_boot[w], Y_boot[w]), timeScale = timeScale)
    #impurity_feuilles <- rbind(impurity_feuilles, c(q, imp1, length(unique(id[w]))))
  }
  return(list(feuilles = id_feuille, id=id_boot, V_split=V_split, impuity=impurete, hist_nodes=hist_nodes, Y_curves= Y_curves, time=time_boot , boot=boot, hist_imp_nodes=hist_imp_nodes))
}



#' Frechet random forest
#'
#' Build a Frechet random Forest
#'
#' @param X [matrix]: Matrix of explanatory variables, each column codes for a variable
#' @param Y [vector]: Output curves
#' @param id [vector]: IDs of measurements
#' @param time [vector]: time of measurements
#' @param mtry [numeric]: number of variables randomly chosen at each split
#' @param ntree [numeric]: number of randomized Frechet trees composing the Frechet random forest
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#' @param ... :
#'
#' @keywords internal
#'
#'
rf_shape <- function(X,Y,id,time,mtry,ntree, timeScale=0.1, ERtrees=FALSE, ...){
  trees <- list()
  if (ERtrees==TRUE){
    for (k in 1:ntree){
      trees[[k]] <- ERtmax(X,Y,id,time, mtry,timeScale=timeScale, ...)
    }
    return(trees)
  }
  for (k in 1:ntree){
    trees[[k]] <- Rtmax(X,Y,id,time, mtry,timeScale=timeScale, ...)
  }
  return(trees)
}


#' Frechet random forests prediction function
#'
#' Given a Frechet tree and new input trajectories predictors \code{X}, this function allows to predict output trajectories \code{Y} associated with the new \code{X} given a Frechet random forest
#'
#' @param object [list]: a list of the randomized Frechet tree composing the Frechet random forest.
#' @inheritParams predict.FrechTree
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(10)
#' data <- DataGenCurves(30)
#' FRF <- Frechforest(data$X,data$Y, data$id,data$time, ntree=40,ncores=2, toPlot="none")
#' pred <- predict(FRF, data$X, data$time, data$id)
#' }
predict.Frechforest <- function(object, X, time, id, timeScale=0.1,...){
  pred.feuille <- matrix(0, length(object$rf), length(unique(id)))
  for (t in 1:length(object$rf)){
    pred.feuille[t,] <- pred.FT(object$rf[[t]], X, time, id, timeScale)
  }
  pred <- NULL
  for (l in 1:dim(pred.feuille)[2]){
    pred_courant <- NULL
    for(k in 1:dim(pred.feuille)[1]){
      pred_courant <- rbind(pred_courant, cbind(rep(k,dim(object$rf[[k]]$Y_curves[[pred.feuille[k,l]]])[1]),object$rf[[k]]$Y_curves[[pred.feuille[k,l]]]))
    }
    predouille <- kmlShape::meanFrechet(pred_courant, timeScale = timeScale)
    #pred.red <- DouglasPeuckerEpsilon(predouille[,1], predouille[,2], 0.05)
    predouille <- cbind(predouille, rep(unique(id)[l],dim(predouille)[1]))
    pred <- rbind(pred, predouille)
  }
  names(pred) <- c("times", "traj", "ID")
  return(pred)
}

#' Parallelized Frechet random forest
#'
#' A parallelized version of "rf_Shape" function
#'
#' @param X [matrix]: Matrix of explanatory variables, each column codes for a variable
#' @param Y [vector]: Output curves
#' @param id [vector]: IDs of measurements
#' @param time [vector]: time of measurements
#' @param mtry [numeric]: number of variables randomly chosen at each split
#' @param ntree [numeric]: number of randomized Frechet trees composing the Frechet random forest
#' @param ncores [numeric]: number of cores used in parallel
#' @param ERT [boolean]: if TRUE uses Extremly randomized Frechet trees to build the Frechet forest
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#' @param ... :
#'
#'
#' @keywords internal
#' @import doRNG
#' @import foreach
#' @import doParallel
#'
rf_shape_para <- function(X,Y,id,time,mtry,ntree, ncores,timeScale=0.1, ERT=FALSE,...){
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  if (ERT==TRUE){
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    trees <- foreach::foreach(k=1:ntree, .packages = "kmlShape") %dorng% {
      #source("forest_shape_utils2.R")
      ERtmax(X,Y,id,time, mtry,timeScale, ...)
    }
    parallel::stopCluster(cl)
    return(trees)
  }

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  trees <- foreach::foreach(k=1:ntree, .packages = "kmlShape") %dorng% {
    #source("forest_shape_utils2.R")
    Rtmax(X,Y,id,time, mtry,timeScale, ...)
  }
  parallel::stopCluster(cl)
  return(trees)
}


#' OOB.tree
#'
#' Computes the oob error of a given tree:
#'
#' @param tree : Frechet tree
#' @param X [matrix]: matrix of the explanatory variables used to build the tree
#' @param Y [vector]: output curves used to build the tree
#' @param id [vector]: IDs of measurements
#' @param time [vector]: time measurements
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#'
#' @keywords internal
#'
OOB.tree <- function(tree, X,Y, id, time, timeScale=0.1){
  BOOT <- tree$boot
  OOB <- setdiff(unique(id), BOOT)
  xerror <- rep(NA,length(OOB))
  for (i in OOB){
    id_w <- which(id== i)
    pred_courant <- pred.FT(tree, X[id_w,,drop=FALSE], time[id_w], id[id_w], timeScale=timeScale)
    #chancla <- DouglasPeuckerNbPoints(tree$Y_curves[[pred_courant]]$times, tree$Y_curves[[pred_courant]]$traj, nbPoints = length(stats::na.omit(Y[id_w])))
    xerror[which(OOB==i)] <- kmlShape::distFrechet(tree$Y_curves[[pred_courant]]$times, tree$Y_curves[[pred_courant]]$traj, time[id_w], Y[id_w], timeScale = timeScale)
  }
  return(mean(xerror))
}


#' OOB error for a Frechet random forest
#'
#' Computes the OOB error of a given Frechet random forest
#'
#' @param rf : Frechet random forest
#' @param X [matrix]: matrix of the explanatory variables used to build the tree
#' @param Y [vector]: output curves used to build the tree
#' @param time [vector]: time measurements
#' @param id [vector]: IDs of measurements
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#'
#' @keywords internal
#'
#'
OOB.rfshape <- function(rf, X, Y,time, id, timeScale=0.1){
  err <- rep(NA,length(unique(id)))
  oob.pred <- list()
  #errdp <- rep(NA,length(unique(id)))
  for (i in 1:length(unique(id))){
    indiv <- unique(id)[i]
    w_y <- which(id==indiv)
    pred_courant <- NULL
    for (t in 1:length(rf)){
      BOOT <- rf[[t]]$boot
      oob <- setdiff(unique(id),BOOT)
      if (is.element(indiv, oob)== TRUE){
        w <- which(id== indiv)
        pred <- pred.FT(rf[[t]], X[w,, drop=FALSE], time[w],id[w], timeScale = timeScale)
        curve <- rf[[t]]$Y_curves[[pred]]
        pred_courant <- rbind(cbind(rep(t,dim(curve)[1]),curve),pred_courant)
      }
    }
    mean_pred <- kmlShape::meanFrechet(pred_courant)
    dp <- as.data.frame(curve.reduc.times(mean_pred$times, mean_pred$traj, time[w_y]))
    names(dp) <- c("x","y")
    #dp <- kmlShape::DouglasPeuckerNbPoints(mean_pred$times, mean_pred$traj, nbPoints = length(stats::na.omit(Y[w_y])))
    #err[i] <- distFrechet(mean_pred$times, mean_pred$traj, time[w_y], Y[w_y])
    oob.pred[[i]] <- dp
    err[i] <- kmlShape::distFrechet(dp$x, dp$y, time[w_y], Y[w_y])
  }
  return(list(err=err,oob.pred=oob.pred))
}


#' Frechet random forest
#'
#' This function builds Frechet random Forest introduced by Capitaine et.al, this includes the OOB predictions, OOB errors and variable importance computations.
#'
#' @inheritParams TreeShapeCV
#' @param mtry [numeric]: number of variables randomly chosen at each split. Default is \code{ncol(X)/3}
#' @param ntree [numeric]: number of randomized Frechet trees composing the Frechet random forest. Default is 100.
#' @param ncores [numeric]: number of cores used in parallel. Default is 1.
#' @param timeScale [numeric]: allow to modify the time scale, increasing or decreasing the cost of the horizontal shift. If timeScale is very big, then the Frechet mean tends to the Euclidean distance. If timeScale is very small, then it tends to the Dynamic Time Warping.
#' @param imp [logical]: TRUE to compute the variables importance FALSE otherwise (default \code{imp=}TRUE)
#' @param ERT [logical]: if TRUE uses Extremly randomized Frechet trees to build the Frechet forest.
#' @param ... : optional parameters to be passed to the low level function.
#'
#'
#' @return a Frechet random forest which is a list of the following elements: \itemize{
#' \item \code{rf:} a list of the \code{ntree} randomized Frechet trees that compose the forest.
#' \item \code{mse :} a vector containing the OOB prediction error of each randomized Frechet tree composing the forest.
#' \item \code{OOB.err: } a vector containing the OOB prediction error of each trajectory in the learning sample.
#' \item \code{OOB.pred: } a list of the OOB prediction of each trajectory in the learning set.
#' \item \code{Importance: } A vector containing the \code{p} variables importance.
#' \item \code{rsq: } “pseudo R-squared”: \code{1 - mse/Var(y)}.
#' }
#' @export
#' @import doRNG
#' @import foreach
#' @import doParallel
#'
#' @examples
#' \dontrun{
#' set.seed(10)
#' data <- DataGenCurves(50)
#' FRF <- Frechforest(data$X,data$Y, data$id,data$time, ntree=40,ncores=2, toPlot="none")
#' }
Frechforest <- function(X,Y,id, time, mtry=ceiling(ncol(X)/3), ntree=100,ncores=1,ERT=FALSE, timeScale=0.1, imp=TRUE, ...){
  print("Building the maximal Frechet trees...")
  debut <- Sys.time()
  rf <-  rf_shape_para(X,Y,id,time, mtry, ntree,timeScale = timeScale,ERT=ERT,ncores=ncores,...)
  temps <- Sys.time() - debut

  print("Forest constucted !")
  xerror <- rep(NA, ntree)
  ### il nous faut l'erreur OOB par arbre :!::

  for (i in 1:ntree){
    xerror[i] <- OOB.tree(rf[[i]], X, Y, id, time, timeScale=timeScale)
  }


  #### on calcule l'erreur OOB de pr?diction
  oob.err <- OOB.rfshape(rf,X,Y,time, id, timeScale=timeScale)

  if (imp == FALSE){
    var.ini <- impurity(Y,time, id, timeScale)
    varex <- 1 - mean(xerror)/var.ini
    frf <- list(rf=rf, error=xerror,OOB.err=oob.err$err,oob.pred= oob.err$oob.pred, rsq=varex)
    class(frf) <- c("Frechforest")
    return(frf)
  }

  print("Importance calculation...")
  p=1
  debut <- Sys.time()
  X.perm <- X
  perm.err <- matrix(NA, ntree, dim(X)[2])
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  Importance <- foreach::foreach(p=1:dim(X)[2],.packages = "kmlShape" ,.combine = "c") %dorng% {

    for (k in 1:ntree){
      BOOT <- rf[[k]]$boot
      nboot <- length(unique(id))- length(BOOT)
      id_boot <- NULL

      for (i in 1:length(BOOT)){
        id_boot <- c(id_boot,which(id==BOOT[i]))
      }

      if (nboot>30){
        X.perm[-id_boot,p] <- permutation(X[-id_boot,p], id[-id_boot])
      }

      if (nboot<=30){
        Nul <- is.na(X.perm[-id_boot,p])
        X.perm[-id_boot,p][Nul==FALSE] <- sample(X.perm[-id_boot,p][Nul==FALSE])
      }

      perm.err[k,p] <- OOB.tree(rf[[k]], X.perm, Y, id, time, timeScale=timeScale)

    }
    ##on remet la variable en place :::
    X.perm[,p] <- X[,p]
    res <- mean(perm.err[,p]- xerror)
  }
  parallel::stopCluster(cl)
  temps.imp <- Sys.time() - debut
  var.ini <- impurity(Y,time, id, timeScale)
  varex <- 1 - mean(xerror)/var.ini
  frf <- list(rf=rf, mse=xerror,OOB.err=oob.err$err,oob.pred= oob.err$oob.pred, Importance=as.vector(Importance), rsq=varex, time=temps)
  class(frf) <- c("Frechforest")
  return(frf)
}



#' Curves data generation function
#'
#' This function allows to ...
#'
#' @param n [numeric]: Number of trajectories
#'
#'
#' @return A list with the followig elements: \itemize{
#' \item \code{X:}
#' \item \code{Y:}
#' \item \code{id:}
#' \item \code{time: }
#' }
#'
#' @import stats
#' @import graphics
#' @export
#'
#' @examples
#' data <- DataGenCurves(50)
DataGenCurves <- function(n){
  tX <- seq(0,1,0.05) # Measurement times for inputs
  tY <- seq(1.1,2.2,0.05) # Measurement times for output

  idX <- NULL  # Input IDs
  idY <- NULL  # Output IDs
  for (i in 1:n){
    idX <- c(idX,rep(i,length(tX)))
    idY <- c(idY, rep(i,length(tY)))
  }

  G1 <- 1*(runif(n)<0.5)  # Groupes of trajectories for the first input variable
  G2 <- 1*(runif(n)<0.5)  # Groupes of trajectories for the second input variable

  # Define the typical temporal behavior functions f.,. for the input variables :

  f11 <- function(t){
    return(0.5*t+0.1*sin(6*t))
  }

  f12 <- function(t){
    return(0.3-0.7*(t-0.45)^2)
  }

  f21 <- function(t){
    return(2*(t-0.5)^2 - t*0.3)
  }

  f22 <- function(t){
    return(0.2 -0.3*t + 0.1*cos(8*t))
  }

  # Simulation of the Input Variables :

  X1 <- NULL
  X2 <- NULL
  beta <- rnorm(n,1, 0.1) # Dilatation terms

  for (i in 1:n){
    X1 <- c(X1, beta[i]*((G1[i]==1)*f11(tX)+(G1[i]==0)*f12(tX)))
    X2 <- c(X2, beta[i]*((G2[i]==1)*f21(tX)+(G2[i]==0)*f22(tX)))
  }

  X1 <- X1+ rnorm(length(X1), 0,0.03) # We add some noise
  X2 <- X2+ rnorm(length(X2), 0,0.03) # We add some noise

  X <- data.frame(X1,X2) # Data matrix

  # Define the typical temporal behavior functions g.,. for the Output :

  g11 <- function(t){
    return(t+ 0.3*sin(10*t))
  }

  g12<- function(t){
    return(t+2*(t-1.7)^2)
  }

  g21 <- function(t){
    return(1.5*exp(-(t-1.5)^2/0.5)- 0.1*t*cos(10*t))
  }

  g22 <- function(t){
    return(2*log(13*(t-1))/(1+t))
  }

  # Output simulation :
  Y <- NULL
  for (i in 1:n){
    Y <- c(Y, beta[i]*((G1[i]==1 & G2[i]==0)*(g11(tY))+ (G1[i]==1 & G2[i]==1)*(g12(tY )) + (G1[i]==0 & G2[i]==1)*(g21(tY)) + (G1[i]==0 & G2[i]==0)*(g22(tY))))
  }

  Y <- Y+ rnorm(length(Y), 0, 0.05) # We add some noise


  par(mfrow=c(1,3))
  w <- which(idX==1)
  plot(tX,X1[w], type="l", ylim=c(min(X1), max(X1)), col="grey",ylab="", xlab="time", main="Input variable X1")
  for (i in 2:n){
    w <- which(idX==i)
    lines(tX,X1[w],col="grey")
  }

  w <- which(idX==1)
  plot(tX,X2[w], type="l", ylim=c(min(X2), max(X2)), col="grey", ylab="", xlab="time", main="Input variable X2")
  for (i in 2:n){
    w <- which(idX==i)
    lines(tX,X2[w], col="grey")
  }


  w <- which(idY==1)
  plot(tY,Y[w], type="l", ylim=c(min(Y), max(Y)), col="grey", xlab="time", ylab="", main = "Output variable Y")
  for (i in 1:n){
    w <- which(idY==i)
    lines(tY,Y[w], col="grey")
  }
  par(mfrow=c(1,1))

  ### The Following code allows to put every element to the right format

  time <- rep(c(tX,tY),n)
  id <- NULL
  for (i in 1:n){
    id <- c(id, rep(i, length(c(tY,tX))))
  }

  M <- NULL
  y <- NULL
  X <- as.matrix(X)
  for (i in 1:n){
    wx <- which(idX==i)
    wy <- which(idY==i)
    y <- c(y, rep(NA,length(wx)),Y[wy])
    M <- rbind(M,rbind(X[wx,],matrix(NA,length(wy),ncol(X))))
  }

  Y <- y # Output
  X <- M # Input

  return(list(X=X,Y=Y,time=time,id=id))
}


#' Title
#'
#' @param time.init
#' @param traj.init
#' @param time.new
#'
#' @keywords internal
#'
curve.reduc.times <- function(time.init , traj.init, time.new){
  new.curve <- matrix(NA,length(time.new),2)
  for (j in 1:length(time.new)){
    w.time <- which.min(abs(time.new[j]-time.init))
    if (round(time.init[w.time]-time.new[j],5)==0){
      new.curve[j,] <- c(time.new[j], traj.init[w.time])
    }
    else {
      t_g <- (time.new[j]>time.init[w.time])*(time.init[w.time]) + (time.new[j]<time.init[w.time])*(time.init[w.time-1])
      t_d <- (time.new[j]<time.init[w.time])*(time.init[w.time]) + (time.new[j]>time.init[w.time])*(time.init[w.time+1])
      Y_g <- (time.new[j]>time.init[w.time])*(traj.init[w.time]) + (time.new[j]<time.init[w.time])*(traj.init[w.time-1])
      Y_d <- (time.new[j]<time.init[w.time])*(traj.init[w.time]) + (time.new[j]>time.init[w.time])*(traj.init[w.time+1])
      d1 <- time.new[j]-t_g
      d2 <- t_d - time.new[j]
      new.curve[j,] <- c(time.new[j], (1 - (d1/(d1+d2)))*Y_g + (1 - (d2/(d1+d2)))*Y_d)
    }
  }
  return(new.curve)
}


