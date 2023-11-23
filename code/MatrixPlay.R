#Lets play around and try and see how we can create a bved and amat

#lets say we are working with 10 assets, 3 are bonds and 7 are equities

b.vec <- c( 1, 0.3, 0.7, rep(0.1, 10), -rep(0.5, 10))

#Lets make a matrix of the weights of different groups that we can bind to the UB and LB
amatr <- matrix(0, nrow = 10, ncol = 3)

amatr[, 1] <- 1

amatr[1:3, 2] <- 1

amatr[4:10, 3] <- 1
#Now combine this with the constraints on the group weights
amat_full <- cbind(amatr, diag(10), -diag(10))



#Now we can even get even more technical and say that within each group the UB and LB will be different