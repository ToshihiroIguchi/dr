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

  #入力をチェック---
  if(is.null(formula)){stop("formula is null.")}
  if(is.null(data)){stop("data is null.")}

  #データを準備---
  ret <- list() #戻り値
  ret$parameter$formula <- formula
  ret$parameter$normalize <- normalize
  ret$parameter$unique <- unique

  #解析用データ
  if(length(formula) == 3){
    data.original <- model.matrix(object = formula[-2], data = data)[,-1] #解析用データ。[,-1]で切片を除く。
  }else{
    data.original <- model.matrix(object = formula, data = data)[,-1] #解析用データ。[,-1]で切片を除く。
  }

  #ユニークなデータを抜粋
  #https://qiita.com/weda_654/items/97c8dbba9f8198845537
  u.no <- !duplicated(data.original)
  data.original <- data.original[u.no,]

  #目的変数がある場合の処理
  if(length(formula) == 3){
    y.formula  <- as.formula(paste0(as.character(formula[2]), "~+0"))
    y.data <- model.frame(y.formula, data)[u.no,1]
    ret$y <- y.data
    formula[2] <- NULL #目的変数をさくじょ
  }else{
    ret$y <- NULL
  }




  #説明変数の項目名を抜粋
  ret$parameter$names  <- names(data.original[1, ])

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
                      scale = normalize, center = normalize,
                      retx = TRUE)
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
  if(!is.na(match("mds", methods))){
    cat("\n[Classical Multidimensional Scaling]\n")
    ret$mds <- cmdscale(d = as.matrix(data.di), k = mds.par$k)
  }

  #t-SNE
  #https://blog.albert2005.co.jp/2015/12/02/tsne/
  #注意点として、t-SNE で圧縮する次元は、 2・3 次元が推奨されています
  if(!is.na(match("tsne", methods))){
    cat("\n[t-Distributed Stochastic Neighbor Embedding]\n")
    ret$tsne <- Rtsne(data.df,
                      dims = tsne.par$dims, initial_dims = tsne.par$initial_dims,
                      perplexity = tsne.par$perplexity, theta = tsne.par$theta)
  }

  #LLE
  #m,kは再考の必要あり。
  if(!is.na(match("lle", methods))){
    cat("\n[Locally linear embedding]\n")
    ret$lle <- lle(data.ma, m = lle.par$m,
                   k = lle.par$k,
                   reg = lle.par$reg, p = lle.par$p)
  }

  #SOM
  #その他のパラメータもあったほうがいい？
  if(!is.na(match("som", methods))){
    cat("\n[Self-Organizing Map]\n")
    ret$som <- som(data.df, xdim=som.par$xdim, ydim=som.par$ydim)
  }

  #diffusion map
  if(!is.na(match("diffmap", methods))){
    cat("\n[Diffusion Map]\n")
    ret$diffmap <- diffuse(data.di,
                           neigen = diffmap.par$neigen,
                           t = diffmap.par$t, maxdim = diffmap.par$maxdim,
                           delta = diffmap.par$delta)
  }


  #戻り値
  class(ret) <- "dr"
  return(ret)
}

ret <- dr(Sepal.Length~Sepal.Length+ Sepal.Width + Petal.Length+ Petal.Width, data = iris)

predict.dr <- function(result, data){
  data.original <- model.matrix(object = result$parameter$formula, data = data)[,-1]

  #モデルで使用した説明変数と、新データで作った説明変数が一致するか確認。
  #カテゴリカル変数まわりでエラーが起きる気がする。
  if(!setequal(result$parameter$names, names(data.original[1, ]))){
    stop("Explanatory variable is different from model.")
  }

  #規格化。data.maが解析用の元データのマトリクス
  if(result$parameter$normalize){
    #規格化。もうちょっとスマートなのがある気がする…
    data.ma <- t(apply(data.original,1,
                       function(x){(x - result$parameter$data.mean)/result$parameter$data.sd}))
  }else{
    data.ma <- data.original
  }

  #変換
  data.df <- data.frame(data.ma) #データフレームに変換

  res <- list()

  res$pca <- predict(result$pca, data.df)
  res$kpca <- predict(result$kpca, data.df)


  #res$tsne <- predict(result$tsne, data.df)
  #https://stackoverflow.com/questions/43377941/t-sne-predictions-in-r
  #t-SNEでpredictできないみたい。
  #res$lle <- predict(result$lle, data.ma)
  #res$som <- predict(result$som, data.df)
}


plot.dr <- function(result, cex = 1){
  if(!is.null(result$y)){
    if(is.factor(result$y)){
      col_p <- rainbow(length(unique(result$y)))
      plot_bg <- col_p[unclass(result$y)]
    }else{
      rank_y <- round(rank(result$y))
      plot_bg <- heat.colors(length(result$y))[rank_y]
    }
  }else{
    plot_bg <- NULL; col_p <- NULL
  }

  if(!is.null(result$pca)){
    plot(result$pca$x[,1],result$pca$x[,2], xlab = "PC1", ylab = "PC2",
         type = "p", pch = 21, cex = cex, bg = plot_bg)
    #https://www1.doshisha.ac.jp/~mjin/R/Chap_07/07.html
  }

  if(!is.null(result$kpca)){
    plot(result$kpca@pcv[,1],result$kpca@pcv[,2], xlab = "K-PC1", ylab = "K-PC2",
         type = "p", pch = 21, cex = cex, bg = plot_bg)
    #https://www1.doshisha.ac.jp/~mjin/R/Chap_07/07.html
  }


}

plot(ret)





