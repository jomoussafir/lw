lw <- function (x, weights = seq(1, 3, length = nt), 
                center = TRUE,
                tol = 0.001,
                debug = FALSE)
  {

    x <- as.matrix(x)
    nassets <- ncol(x)
    nt <- nrow(x)
    
    ## weights <- weights/mean(weights)
    weights <- rep(1, nt)
    
    center <- colMeans(x * weights, na.rm = TRUE)    
    xcen <- sweep(x, 2, center, "-")

    svar <- var(xcen * sqrt(weights), use = "pairwise")
    sdiag <- sqrt(diag(svar))
    scor <- cov2cor(svar)
    cor.avg <- mean(cov2cor(svar)[lower.tri(svar, diag = FALSE)])
    prior <- svar
    prior[] <- cor.avg
    diag(prior) <- 1
    prior <- sdiag * prior * rep(sdiag, each = nassets)
    gamma <- sum((svar - prior)^2)
    
    if(debug)
      {
        cat(paste('svar\n'))
        print(svar)
        cat(paste('\n'))
        
        cat(paste('sdiag\n'))
        print(sdiag)
        cat(paste('\n'))

        cat(paste('scor\n'))
        print(scor)
        cat(paste('\n'))
        
        cat(paste('avg all cor\n'))
        print(mean(cov2cor(svar)))
        cat(paste('\n'))
        
        cat(paste('cor.avg\n'))
        print(cor.avg)
        cat(paste('\n'))
        
        cat(paste('xcen\n'))
        print(xcen)
        cat(paste('\n'))
        
        cat(paste('prior\n'))
        print(prior)
        cat(paste('\n'))

        cat(paste('gamma\n'))
        print(gamma)
        cat(paste('\n'))

        cat(paste('w\n'))
        print(weights)
        cat(paste('\n'))
      }
    
    ## compute shrinkage
    pi.mat <- array(0, c(nassets, nassets))
    theta1.mat <- array(0, c(nassets, nassets))
    theta2.mat <- array(0, c(nassets, nassets))
    
    
    for (i in 1:nt)
      {
        this.wt <- weights[i]
        this.cross <- crossprod(xcen[i, , drop = FALSE]) - svar
        pi.mat <- pi.mat + this.wt * this.cross^2
        theta1.mat <- theta1.mat + this.wt * diag(this.cross) * this.cross
        theta2.mat <- theta2.mat + rep(diag(this.cross), 
                                       each = nassets) * this.cross * this.wt
      }
    
    pi.hat <- sum(pi.mat)/nt
    theta1.mat <- theta1.mat/nt
    theta2.mat <- theta2.mat/nt
    diag(theta1.mat) <- 0
    diag(theta2.mat) <- 0
    theta1.mat <- rep(sdiag, each = nassets) * theta1.mat/sdiag
    theta2.mat <- rep(1/sdiag, each = nassets) * theta2.mat * sdiag

    if(debug)
      {
        cat(paste('Pi\n'))
        print(pi.mat)
        cat(paste('\n'))
        
        cat(paste('Theta1\n'))
        print(theta1.mat)
        cat(paste('\n'))
        
        cat(paste('Theta2\n'))
        print(theta2.mat)
        cat(paste('\n'))
      }
    
    rho <- sum(diag(pi.mat))/nt + cor.avg * 0.5 * (sum(theta1.mat) + sum(theta2.mat))
    delta <- (pi.hat - rho)/gamma/nt
    delta <- min(1, max(0, delta))

    if(debug)
      {
        cat(paste('rho\n'))
        print(rho)
        cat(paste('\n'))
        
        cat(paste('delta\n'))
        print(delta)
        cat(paste('\n'))
      }
    
    ans <- delta * prior + (1 - delta) * svar
    eig <- eigen(ans, symmetric = TRUE)
    tol <- tol * max(eig$values)
    if (min(eig$values) < tol)
      {
        ## cat(paste('small eigval.\n'))
        vals <- eig$values
        vals[vals < tol] <- tol
        ans[] <- eig$vectors %*% (vals * t(eig$vectors))
      }
    return(ans)
}
