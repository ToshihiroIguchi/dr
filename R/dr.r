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

dr <- function(formula, data, normalize = TRUE){
  #入力をチェック
  if(is.null(formula)){stop("formula is null.")}
  if(is.null(data)){stop("data is null.")}
  if(length(formula) == 3){stop("You have mistyped the format of formula.")}

  #データを準備
  ret <- list() #戻り値
  ret$formula <- formula

  #解析用データ
  #まずここをなおす。
  data.original <- model.matrix(object = formula, data = data)[,-1] #解析用データ。[,-1]で切片を除く。

  if(normalize){
    data.mean <- apply(data.original, 2, mean)
    data.sd <- apply(data.original, 2, sd)
    data.ma <- data.original
  }else{
    data.mean <- NULL
    data.sd <- NULL
    data.ma <- data.original
  }



  data.df <- data.frame(data.ma) #データフレームに変換
  data.di <- dist(data.ma) #dist classの距離に変換

  #PCA(数値のみしか使用できない。)
  ret$pca <- prcomp(~., data = data.df,
                    scale = normalize, center = normalize)

  #kernel PCA
  #kernelとkparは後で修正。
  ret$kpca <- kpca(x = ~., data = data.df,
               kernel = "rbfdot", kpar = list(sigma = 0.1))

  #MDS
  #kは後で修正。
  ret$mds <- cmdscale(d = as.matrix(data.di), k = 1)

  #t-SNE
  #https://blog.albert2005.co.jp/2015/12/02/tsne/
  #注意点として、t-SNE で圧縮する次元は、 2・3 次元が推奨されています。
  #パラメータはあとで修正

  ret$tsne <- Rtsne(unique(data.df), dims = 2, initial_dims = 50, perplexity = 30, theta = 0.5)

  #LLE
  #m,kはあとで修正
  ret$lle <- lle(data.ma, m = 2, k = 2)

  #SOM
  #normalizeの有無はオプションにしたほうがよいか？
  #パラメータはあとで修正。
  ret$som <- som(normalize(data.df), xdim=8, ydim=4)

  #diffusion map
  #scaleの有無
  #その他のパラメータ
  ret$diffmap <- diffuse(dist(scale(iris[,-5])))

  #戻り値
  class(ret) <- "dr"
  return(ret)
}

ret <- dr(~Sepal.Length+ Sepal.Width + Petal.Length+ Petal.Width+Species, data = iris)








#isomap
#isomapうまくいかない
test <- isomap(vegdist((iris[,-5])),k = 5)

## The following examples also overlay minimum spanning tree to
## the graphics in red.
data(BCI)
dis <- vegdist(BCI)
ord <- isomap(dis, k=3)


