\name{kappas}
\alias{kappas}
\title{Condition number of the source matrix as a function of a variable index}
\usage{
kappas(mat, vars)
}
\arguments{
  \item{mat}{a compositional vector or dataframe of
  row-wise compositions}
}
\value{
A \code{l=length(vars)} vector with condition numbers of
\code{mat[vars[1:i]]} for \code{i=1:l}
}
\description{
Condition number of the source matrix as a function of a
variable index
}

