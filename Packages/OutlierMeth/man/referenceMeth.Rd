\name{referenceMeth}
\alias{referenceMeth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reference threshold function
}
\description{
Makes a reference threshold file
}
\usage{
referenceMeth(beta)
}
\arguments{
  \item{beta}{
Matrix of Illumina 450k methylation data
}
}
\details{
Makes a reference threshold file
}
\value{
Makes a file with probe names, thresholds, and skew

}
\references{
NA
}
\author{
Bradley Downs and Leslie Cope
}

\examples{
ref=referenceMeth(example_data)

function (beta) 
{
    options(scipen = 999)
    Hyper_Norm1 = apply(beta, 1, quantile, 0.99, na.rm = T)
    Hyper_Norm2 = apply(beta, 1, quantile, 0.999, na.rm = T)
    Hyper_Norm3 = apply(beta, 1, quantile, 0.9999, na.rm = T)
    Hyper_Norm4 = apply(beta, 1, quantile, 0.99999, na.rm = T)
    Hypo_Norm1 = apply(beta, 1, quantile, 0.01, na.rm = T)
    Hypo_Norm2 = apply(beta, 1, quantile, 0.001, na.rm = T)
    Hypo_Norm3 = apply(beta, 1, quantile, 0.0001, na.rm = T)
    Hypo_Norm4 = apply(beta, 1, quantile, 0.00001, na.rm = T)
    Hyper_Norm1 = as.data.frame(Hyper_Norm1)
    Hyper_Norm2 = as.data.frame(Hyper_Norm2)
    Hyper_Norm3 = as.data.frame(Hyper_Norm3)
    Hyper_Norm4 = as.data.frame(Hyper_Norm4)
    Hypo_Norm1 = as.data.frame(Hypo_Norm1)
    Hypo_Norm2 = as.data.frame(Hypo_Norm2)
    Hypo_Norm3 = as.data.frame(Hypo_Norm3)
    Hypo_Norm4 = as.data.frame(Hypo_Norm4)
    bind_info = cbind(Hyper_Norm1, Hyper_Norm2, Hyper_Norm3, 
        Hyper_Norm4, Hypo_Norm1, Hypo_Norm2, Hypo_Norm3, Hypo_Norm4)
    bind_info = na.omit(bind_info)
    names(bind_info)[1] <- "P0.01"
    names(bind_info)[2] <- "P0.001"
    names(bind_info)[3] <- "P0.0001"
    names(bind_info)[4] <- "P0.00001"
    names(bind_info)[5] <- "N0.01"
    names(bind_info)[6] <- "N0.001"
    names(bind_info)[7] <- "N0.0001"
    names(bind_info)[8] <- "N0.00001"
    return(bind_info)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
