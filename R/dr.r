#Dimensionality reduction
#次元削減

#参考
#https://www.slideshare.net/mikayoshimura50/150905-wacode-2nd
#http://d.hatena.ne.jp/hoxo_m/20120313/p1
#https://blog.albert2005.co.jp/2014/12/11/%E9%AB%98%E6%AC%A1%E5%85%83%E3%83%87%E3%83%BC%E3%82%BF%E3%81%AE%E5%8F%AF%E8%A6%96%E5%8C%96%E3%81%AE%E6%89%8B%E6%B3%95%E3%82%92swiss-roll%E3%82%92%E4%BE%8B%E3%81%AB%E8%A6%8B%E3%81%A6%E3%81%BF%E3%82%88/

library(kernlab)
library(Rtsne)
library(lle)
library(som)
library(diffusionMap)

dr <- function(formula, data, normalize = TRUE, unique = TRUE,
               methods = c("pca", "kpca" ,"mds", "tsne", "lle", "som"),
               kpca.par = list(kernel = "rbfdot", sigma = 0.1),
               mds.par = list(k = 2),
               tsne.par = list(dims = 2, initial_dims = 50, perplexity = 30, theta = 0.5),
               lle.par = list(m = 2, k = 5, reg = 2, p = 0.5),
               som.par = list(xdim = 5, ydim = 6),
               diffmap.par = list(neigen = NULL, t = 0, maxdim = 50, delta=10^-5))
  {

  #kpca.par はsigma,degree, scale, offset,orderも指定できる。
  #tsne.parのdimsは3以下がよい。

  #入力をチェック
  if(is.null(formula)){stop("formula is null.")}
  if(is.null(data)){stop("data is null.")}
  if(length(formula) == 3){stop("You have mistyped the format of formula.")}

  #データを準備
  ret <- list() #戻り値
  ret$parameter$formula <- formula
  ret$parameter$normalize <- normalize
  ret$parameter$unique <- unique

  #解析用データ
  data.original <- model.matrix(object = formula, data = data)[,-1] #解析用データ。[,-1]で切片を除く。
  ret$parameter$data.head  <- data.original[1, ]
  if(unique){data.original <- unique(data.original)} #ユニークなデータだけ抜粋。

  #規格化。data.maが解析用の元データのマトリクス
  if(normalize){
    data.mean <- apply(data.original, 2, mean)
    data.sd <- apply(data.original, 2, sd)
    #規格化。もうちょっとスマートなのがある気がする…
    data.ma <- t(apply(data.original,1,function(x){(x - data.mean)/data.sd}))
  }else{
    data.mean <- NULL
    data.sd <- NULL
    data.ma <- data.original
  }
  ret$parameter$data.mean <- data.mean
  ret$parameter$data.sd <- data.sd

  #変換
  data.df <- data.frame(data.ma) #データフレームに変換
  data.di <- dist(data.ma) #dist classの距離に変換

  #ここから解析本体-----------------------------------------------

  #PCA
  #(数値のみしか使用できない。)
  if(!is.na(match("pca", methods))){
    cat("\n[Principal Components Analysis] \n")
    ret$pca <- prcomp(~., data = data.df,
                      scale = normalize, center = normalize)
  }



  #kernel PCA
  if(!is.na(match("kpca", methods))){
    cat("\n[Kernel Principal Components Analysis] \n")
    #kpcaに渡すデータを編集。
    kpca.kernel <- kpca.par$kernel
    kpca.par$kernel <- NULL

    #kpca本体
    ret$kpca <- kpca(x = ~., data = data.df,
                     kernel = kpca.kernel, kpar = kpca.par)
  }


  #MDS
  cat("\n[Classical Multidimensional Scaling]\n")
  ret$mds <- cmdscale(d = as.matrix(data.di), k = mds.par$k)

  #t-SNE
  #https://blog.albert2005.co.jp/2015/12/02/tsne/
  #注意点として、t-SNE で圧縮する次元は、 2・3 次元が推奨されています
  cat("\n[t-Distributed Stochastic Neighbor Embedding]\n")
  ret$tsne <- Rtsne(data.df,
                    dims = tsne.par$dims, initial_dims = tsne.par$initial_dims,
                    perplexity = tsne.par$perplexity, theta = tsne.par$theta)

  #LLE
  #m,kは再考の必要あり。
  cat("\n[Locally linear embedding]\n")
  ret$lle <- lle(data.ma, m = lle.par$m,
                 k = lle.par$k,
                 reg = lle.par$reg, p = lle.par$p)

  #SOM
  #その他のパラメータもあったほうがいい？
  cat("\n[Self-Organizing Map]\n")
  ret$som <- som(data.df, xdim=som.par$xdim, ydim=som.par$ydim)

  #diffusion map
  cat("\n[Diffusion Map]\n")
  ret$diffmap <- diffuse(data.di,
                         neigen = diffmap.par$neigen,
                         t = diffmap.par$t, maxdim = diffmap.par$maxdim,
                         delta = diffmap.par$delta)




  #戻り値
  class(ret) <- "dr"
  return(ret)
}

ret <- dr(~Sepal.Length+ Sepal.Width + Petal.Length+ Petal.Width+Species, data = iris)




