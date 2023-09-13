
getPermuteMatrix <- function(perm, N,  strata = NULL)
  {
    ## 'perm' is either a single number, a how() structure or a
    ## permutation matrix
    if (length(perm) == 1) {
      perm <- how(nperm = perm)
    }
    ## apply 'strata', but only if possible: ignore silently other cases
    if (!missing(strata) && !is.null(strata)) {
      if (inherits(perm, "how") && is.null(getBlocks(perm)))
        setBlocks(perm) <- strata
    }
    ## now 'perm' is either a how() or a matrix
    if (inherits(perm, "how"))
      perm <- shuffleSet(N, control = perm)
    else { # matrix: check that it *strictly* integer
      if(!is.integer(perm) && !all(perm == round(perm)))
        stop("permutation matrix must be strictly integers: use round()")
    }
    ## now 'perm' is a matrix (or always was). If it is a plain
    ## matrix, set minimal attributes for printing. This is a dirty
    ## kluge: should be handled more cleanly.
    if (is.null(attr(perm, "control")))
      attr(perm, "control") <-
        structure(list(within=list(type="supplied matrix"),
                       nperm = nrow(perm)), class = "how")
    perm
  }

simper_modify <- function (comm, group, permutations = 999, parallel = 1, ...)
{
  if (!missing(parallel))
    .NotYetUsed("parallel", error = FALSE)
  EPS <- sqrt(.Machine$double.eps)
  comm <- as.matrix(comm)
  tri <- outer(seq_len(nrow(comm)), seq_len(nrow(comm)), ">")
  rs <- rowSums(comm)
  rs <- outer(rs, rs, "+")[tri]
  spcontr <- sapply(seq_len(ncol(comm)), function(i) as.vector(vegdist(comm[,
                                                                            i, drop = FALSE], "man")))
  spcontr <- sweep(spcontr, 1, rs, "/")
  colnames(spcontr) <- colnames(comm)
  outlist <- NULL
  if (missing(group) || length(unique(group)) == 1) {
    nperm <- 0
    permat <- NULL
    average <- colMeans(spcontr)
    overall <- sum(average)
    sdi <- apply(spcontr, 2, sd)
    ord <- order(average, decreasing = TRUE)
    cusum <- cumsum(average[ord])/overall
    outlist[["total"]] <- list(species = colnames(comm),
                               average = average, overall = overall, sd = sdi,
                               ratio = average/sdi, ava = NULL, avb = NULL, ord = ord,
                               cusum = cusum, p = NULL)
  }
  else {
    comp <- t(combn(as.character(unique(group)), 2))
    spavg <- apply(comm, 2, function(x) tapply(x, group,
                                               mean))
    contrmatch <- function(X, Y, patt) X != Y & X %in% patt &
      Y %in% patt
    for (i in seq_len(nrow(comp))) {
      tmat <- outer(group, group, FUN = contrmatch, patt = comp[i,
      ])
      take <- tmat[tri]
      average <- colMeans(spcontr[take, ,drop = FALSE])
      overall <- sum(average)
      sdi <- apply(spcontr[take, , drop = FALSE], 2, sd)
      ratio <- average/sdi
      ord <- order(average, decreasing = TRUE)
      cusum <- cumsum(average[ord])/overall
      ava <- spavg[comp[i, 1], ]
      avb <- spavg[comp[i, 2], ]
      permat <- getPermuteMatrix(permutations, nrow(comm))
      nperm <- nrow(permat)
      if (nperm) {
        Pval <- rep(1, ncol(comm))
        for (k in seq_len(nperm)) {
          take <- tmat[permat[k, ], permat[k, ]][tri]
          Pval <- Pval + ((colMeans(spcontr[take, ,drop = FALSE]) -
                             EPS) >= average)
        }
        Pval <- Pval/(nperm + 1)
      }
      else {
        Pval <- NULL
      }
      outlist[[paste(comp[i, ], collapse = "_")]] <- list(species = colnames(comm),
                                                          average = average, overall = overall, sd = sdi,
                                                          ratio = ratio, ava = ava, avb = avb, ord = ord,
                                                          cusum = cusum, p = Pval)
    }
  }
  class(outlist) <- "simper"
  attr(outlist, "permutations") <- nperm
  attr(outlist, "control") <- attr(permat, "control")
  outlist
}

simper_result <- simper_modify(rarefy_V4_clean_wider, V4_meta$`adpater no`, permutations = 999)
saveRDS(object, file = "my_data.rds")
