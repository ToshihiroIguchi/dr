\name{dr}
\alias{dr}
\title{Dimensionality reduction}
\usage{
dr(formula, data, normalize = TRUE, unique = TRUE,
               methods = c("pca", "kpca" ,"mds", "tsne", "lle", "som", "diffmap"),
               kpca.par = list(kernel = "rbfdot", sigma = 0.1),
               mds.par = list(k = 2),
               tsne.par = list(dims = 2, initial_dims = 50, perplexity = 30, theta = 0.5),
               lle.par = list(m = 2, k = 5, reg = 2, p = 0.5),
               som.par = list(xdim = 5, ydim = 6),
               diffmap.par = list(neigen = NULL, t = 0, maxdim = 50, delta=10^-5))
}
\description{
Wrappers of various dimension reduction algorithms.
Data after dimension reduction can also be easily visualized.
}

\arguments{

\item{formula}{
an object of class "formula" : a symbolic description of the model whose dimension is reduced.
Explanatory variables may be present or absent, and may be factors or numerical values.
}

\item{data}{
an optional data frame containing the variables in the model.
}

}

\author{
Toshihiro Iguchi
}

\examples{
result <- dr(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data = iris)
plot(result)
}
