#' @importFrom stats median
.make_hexbin_function <- function(x, action, cID, na.rm=FALSE) {
  if (action == "majority") {
    func_if <- !(is.factor(x) | is.character(x))

    if (func_if) {
      stop("For action 'majority' x needs to be a factor or character.")
    } else {
      res <- tapply(x, cID, FUN = function(z) {
        names(sort(table(z),
          decreasing = TRUE
        )[1])
      })
      res <- as.factor(res)
      return(res)
    }
  }

  if (action == "prop") {
    func_if <- !(is.factor(x) | is.character(x))

    if (func_if) {
      stop("For action 'prop' x needs to be a factor or character.")
    } else {
      nrows <- length(unique(cID))
      res <- vapply(unique(x),
        FUN.VALUE = rep(0, length = nrows),
        function(y) {
          tapply(x, cID, FUN = function(z) {
            sum(z == y, na.rm = na.rm) / sum(!is.na(z))
          })
        }
      )
      res <- apply(res, 2, as.numeric)
      return(res)
    }
  }

  if (action == "median") {
    func_if <- !is.numeric(x)

    if (func_if) {
      stop("For action 'median' x needs to be numeric")
    } else {
      res <- tapply(x, cID, FUN = function(z) median(z, na.rm = na.rm))
      res <- as.numeric(res)
      return(res)
    }
  }

  if (action == "mode") {
    func_if <- !is.numeric(x)

    if (func_if) {
      stop("For action 'mode' x needs to be numeric")
    } else {
      res <- tapply(x, cID, FUN = function(z) .get_mode(z))
      res <- as.numeric(res)
      return(res)
    }
  }

  if (action == "prop_0") {
    func_if <- !is.numeric(x)

    if (func_if) {
      stop("For action 'prop_0' x needs to be numeric")
    } else {
      res <- tapply(x, cID, FUN = function(z) sum(z > 0, na.rm = na.rm) / 
          sum(!is.na(z)))
      res <- as.numeric(res)
      return(res)
    }
  }
  
  if (action == "sd") {
    func_if <- !is.numeric(x)
    
    if (func_if) {
      stop("For action 'prop_0' x needs to be numeric")
    } else {
      res <- tapply(x, cID, FUN = function(z) sd(z, na.rm = na.rm))
      res <- as.numeric(res)
      return(res)
    }
  }

  if (action == "mean") {
    func_if <- !is.numeric(x)

    if (func_if) {
      stop("For action 'mean' x needs to be numeric")
    } else {
      res <- tapply(x, cID, FUN = function(z) mean(z, na.rm = na.rm))
      res <- as.numeric(res)
      return(res)
    }
  } else {
    stop("Specify valid action!")
  }
}

.make_hexbin_colnames <- function(x, name_s) {
  if (is.character(x)) {
    paste0(name_s, "_prop_", unique(x))
  } else {
    paste0(name_s, "_prop_", levels(x))
  }
}

.get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

utils::globalVariables(c("dr", "x", "assayNames", "assays", "groups", "y"))

