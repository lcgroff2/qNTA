##Louis Groff - 03/10/2021:
##Custom Function list for bestNormalize transformation selections

## Define user-function for cuberoot
cuberoot_x <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^(1/3)
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('cuberoot_x', class(val))
  val
}

##predict method for cuberoot:
predict.cuberoot_x <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse-cube-root (cube)
    newdata <-  newdata^3 - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take cube root
    newdata <- (newdata + object$a)^(1/3)
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for cuberoot
print.cuberoot_x <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'cuberoot(x + a) Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- list(
  cuberoot_x = cuberoot_x,
  predict.cuberoot_x = predict.cuberoot_x,
  print.cuberoot_x = print.cuberoot_x)

## Define user-function for fourthroot
fourthroot_x <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^(1/4)
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('fourthroot_x', class(val))
  val
}

##predict method for fourthroot:
predict.fourthroot_x <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse-fourth-root
    newdata <-  newdata^4 - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take fourth root
    newdata <- (newdata + object$a)^(1/4)
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for fourthroot
print.fourthroot_x <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'fourthroot(x + a) Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  fourthroot_x = fourthroot_x,
  predict.fourthroot_x = predict.fourthroot_x,
  print.fourthroot_x = print.fourthroot_x))

## Define user-function for fifthroot
fifthroot_x <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^(1/5)
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('fifthroot_x', class(val))
  val
}

##predict method for fifthroot:
predict.fifthroot_x <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse-fifth-root
    newdata <-  newdata^5 - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take fourth root
    newdata <- (newdata + object$a)^(1/5)
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for fifthroot
print.fifthroot_x <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'fifthroot(x + a) Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  fifthroot_x = fifthroot_x,
  predict.fifthroot_x = predict.fifthroot_x,
  print.fifthroot_x = print.fifthroot_x))

## Define user-function for sixthroot
sixthroot_x <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^(1/6)
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('sixthroot_x', class(val))
  val
}

##predict method for sixthroot:
predict.sixthroot_x <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse-sixth-root
    newdata <-  newdata^6 - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take fourth root
    newdata <- (newdata + object$a)^(1/6)
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for sixthroot
print.sixthroot_x <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'sixthroot(x + a) Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  sixthroot_x = sixthroot_x,
  predict.sixthroot_x = predict.sixthroot_x,
  print.sixthroot_x = print.sixthroot_x))

## Define user-function for inverse square root
invsqrt_x <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^(-1/2)
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('invsqrt_x', class(val))
  val
}

##predict method for inverse square root:
predict.invsqrt_x <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse inverse square root
    newdata <-  newdata^-2 - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take fourth root
    newdata <- (newdata + object$a)^(-1/2)
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for inverse square root
print.invsqrt_x <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'invsqrt(x + a) Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  invsqrt_x = invsqrt_x,
  predict.invsqrt_x = predict.invsqrt_x,
  print.invsqrt_x = print.invsqrt_x))

## Define user-function for inverse cube root
invcuberoot_x <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^(-1/3)
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('invcuberoot_x', class(val))
  val
}

##predict method for inverse cube root:
predict.invcuberoot_x <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse inverse cube root
    newdata <-  newdata^-3 - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take fourth root
    newdata <- (newdata + object$a)^(-1/3)
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for inverse cube root
print.invcuberoot_x <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'invcuberoot(x + a) Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  invcuberoot_x = invcuberoot_x,
  predict.invcuberoot_x = predict.invcuberoot_x,
  print.invcuberoot_x = print.invcuberoot_x))

## Define user-function for inverse 4th root
inv4throot_x <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^(-1/4)
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('inv4throot_x', class(val))
  val
}

##predict method for inverse 4th root:
predict.inv4throot_x <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse inverse 4th root
    newdata <-  newdata^-4 - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take inverse fourth root
    newdata <- (newdata + object$a)^(-1/4)
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for inverse 4th root
print.inv4throot_x <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'inv4throot(x + a) Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  inv4throot_x = inv4throot_x,
  predict.inv4throot_x = predict.inv4throot_x,
  print.inv4throot_x = print.inv4throot_x))

## Define user-function for inverse 5th root
inv5throot_x <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^(-1/5)
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('inv5throot_x', class(val))
  val
}

##predict method for inverse 5th root:
predict.inv5throot_x <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse inverse 4th root
    newdata <-  newdata^-5 - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take inverse fifth root
    newdata <- (newdata + object$a)^(-1/5)
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for inverse 5th root
print.inv5throot_x <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'inv5throot(x + a) Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  inv5throot_x = inv5throot_x,
  predict.inv5throot_x = predict.inv5throot_x,
  print.inv5throot_x = print.inv5throot_x))

## Define user-function for inverse 6th root
inv6throot_x <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^(-1/6)
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('inv6throot_x', class(val))
  val
}

##predict method for inverse 6th root:
predict.inv6throot_x <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse inverse 6th root
    newdata <-  newdata^-6 - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take inverse 6th root
    newdata <- (newdata + object$a)^(-1/6)
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for inverse 6th root
print.inv6throot_x <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'inv6throot(x + a) Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  inv6throot_x = inv6throot_x,
  predict.inv6throot_x = predict.inv6throot_x,
  print.inv6throot_x = print.inv6throot_x))

## Define user-function for inverse x
inv_x <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^-1
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('inv_x', class(val))
  val
}

##predict method for inverse x:
predict.inv_x <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse inverse 
    newdata <-  newdata^-1 - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take inverse 
    newdata <- (newdata + object$a)^-1
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for inverse x
print.inv_x <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'inv(x + a) Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  inv_x = inv_x,
  predict.inv_x = predict.inv_x,
  print.inv_x = print.inv_x))

## Define user-function for inverse x^2
inv_x2 <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^-2
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('inv_x2', class(val))
  val
}

##predict method for inverse x^2:
predict.inv_x2 <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse inverse square
    newdata <-  newdata^(-1/2) - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take inverse square
    newdata <- (newdata + object$a)^-2
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for inverse x^2
print.inv_x2 <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'inv(x + a)^2 Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  inv_x2 = inv_x2,
  predict.inv_x2 = predict.inv_x2,
  print.inv_x2 = print.inv_x2))

## Define user-function for inverse x^3
inv_x3 <- function(x, a = NULL, standardize = FALSE, ...) {
  stopifnot(is.numeric(x))
  
  min_a <- max(0, -(min(x, na.rm = TRUE)))
  if(!length(a)) 
    a <- min_a
  if(a < min_a) {
    warning("Setting a <  max(0, -(min(x))) can lead to transformation issues",
            "Standardize set to FALSE")
    standardize <- FALSE
  }
  
  
  x.t <- (x + a)^-3
  mu <- mean(x.t, na.rm = TRUE)
  sigma <- sd(x.t, na.rm = TRUE)
  if (standardize) x.t <- (x.t - mu) / sigma
  
  # Get in-sample normality statistic results
  ptest <- nortest::pearson.test(x.t)
  
  val <- list(
    x.t = x.t,
    x = x,
    mean = mu,
    sd = sigma,
    a = a,
    n = length(x.t) - sum(is.na(x)),
    norm_stat = unname(ptest$statistic / ptest$df),
    standardize = standardize
  )
  
  # Assign class
  class(val) <- c('inv_x3', class(val))
  val
}

##predict method for inverse x^3:
predict.inv_x3 <- function(object, newdata = NULL, inverse = FALSE, ...) {
  
  # If no data supplied and not inverse
  if (is.null(newdata) & !inverse)
    newdata <- object$x
  
  # If no data supplied and inverse
  if (is.null(newdata) & inverse)
    newdata <- object$x.t
  
  # Actually performing transformations
  
  # Perform inverse transformation as estimated
  if (inverse) {
    
    # Reverse-standardize
    if (object$standardize) 
      newdata <- newdata * object$sd + object$mean
    
    # Reverse inverse cube
    newdata <-  newdata^(-1/3) - object$a
    
    
    # Otherwise, perform transformation as estimated
  } else if (!inverse) {
    # Take inverse cube
    newdata <- (newdata + object$a)^-3
    
    # Standardize to mean 0, sd 1
    if (object$standardize) 
      newdata <- (newdata - object$mean) / object$sd
  }
  
  # Return transformed data
  unname(newdata)
}

##print function for inverse x^2
print.inv_x3 <- function(x, ...) {
  cat(ifelse(x$standardize, "Standardized", "Non-Standardized"),
      'inv(x + a)^3 Transformation with', x$n, 'nonmissing obs.:\n', 
      'Relevant statistics:\n',
      '- a =', x$a, '\n',
      '- mean (before standardization) =', x$mean, '\n',
      '- sd (before standardization) =', x$sd, '\n')
}

## add to custom function list for implementation:

custom_transform <- append(custom_transform, list(
  inv_x3 = inv_x3,
  predict.inv_x3 = predict.inv_x3,
  print.inv_x3 = print.inv_x3))