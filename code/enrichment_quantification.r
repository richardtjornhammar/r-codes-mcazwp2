#Copyright 2020 RICHARD TJÖRNHAMMAR
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

library ( 'clusterProfiler' )
library ( 'DOSE' )
library ( 'enrichplot' )
library ( 'org.Hs.eg.db' )

doClusterProfilingRehash <- function( filename="../data/corr_de_analysis_bas.csv",
                                      cutoff = 0.05, use=".q" , addition=NULL ) {
  library ( 'clusterProfiler' )
  library ( 'org.Hs.eg.db' )
  res_tab <- read.table( filename )
  O <- list()
  for ( value in colnames(res_tab) )
  { 
      if ( grepl( use,value ) )
      {
          cluster <- res_tab[ res_tab[ c(value) ]<cutoff,c("name") ][-1]
          gsrename_ <- bitr ( cluster , fromType = "SYMBOL",
              toType = c("ENTREZID") ,
              OrgDb = org.Hs.eg.db  )
          entid <- gsrename_[,2]
          O[[value]] <- entid
      }
  }
  return ( O )
}

create_namespace <- function( cluster ) {
    library ( 'clusterProfiler' )
    library ( 'org.Hs.eg.db' )
    gsrename_ = bitr ( cluster , fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db  )
    namespace = gsrename_[,2]
    return( namespace )
}

namespace_from_file <- function ( filename = "../data/corr_de_analysis_bas.csv" ,
        value = "B.M.p" , cutoff=1. , name="name" ) {
    value = c(value)
    res_tab = read.table( filename )
    cluster = res_tab[ res_tab[ c(value) ]<cutoff ,c(name) ][-1]
    namespace = create_namespace ( cluster )
    return ( namespace )
}

translated_values_from_file <- function ( filename = "../data/corr_de_analysis_bas.csv" ,
        value = "B.M.p" , name="name" ) {

    library ( 'org.Hs.eg.db' )
    
    res_tab = read.table( filename )
    ld <- bitr ( as.vector( res_tab$name ) ,
            fromType = "SYMBOL"      ,
            toType   = c("ENTREZID") ,
            OrgDb    = org.Hs.eg.db  )

    val_df = NULL
    for (v in rownames(res_tab) ) {
        tmp_name = ld[ ld[,1]==res_tab[v,name],2 ][1]
        if (!is.na(tmp_name)) {
            tmp_df <- data.frame(res_tab[v,value],
                        row.names=ld[ ld[,1]==res_tab[v,name],2 ][1] )
            if (is.null(val_df)) {
                val_df <- tmp_df
            } else {
                val_df <- rbind( val_df, tmp_df )
            }
        }
    }
    return( val_df )
}

cluster_enrichment <- function( O , exclude = NULL ) {
    library ( 'clusterProfiler' )
    library (   'org.Hs.eg.db'  )

    if ( !is.null(exclude) ) {
        C <- list()
        for ( value in names( O ) ) {
            if ( !grepl(exclude, value) ) {
                C[[value]] <- O[[value]]
            }
        }
    } else {
        C <- O
    }
    ck = compareCluster(geneCluster = C, fun = "enrichGO", OrgDb='org.Hs.eg.db' )
    return(ck)
}

create_cnet_plot <- function( ego , filename = "../results/cnetPathways.pdf" , foldchange=NULL ) {
    library ( 'clusterProfiler' )
    library ( 'DOSE' )
    library ( 'org.Hs.eg.db' )
    library ( 'Cairo' )
    Cairo( file = filename , type="pdf" , width=6000 , height=6000 , units="px", dpi=200 )
    egox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID' )
    p <- cnetplot(
      egox ,
      showCategory = 5   ,
      foldChange = foldchange ,
      layout = "kk"      ,
      colorEdge = TRUE   ,
      circular = FALSE   ,
      node_label = "all"
    )
    print ( p )
}

draw_upset_plot <- function ( ego, filename = "../results/upsetPathways.pdf" ) {
    library ( 'enrichplot' )
    library ( 'Cairo' )
    Cairo( file = filename , type="pdf" , width=6000 , height=3000 , units="px", dpi=200 )
    u = upsetplot( ego, text.scale = 3. )
    print( u )
}

draw_dot_plot <- function( ck , filename="../results/cluster_Pathways.pdf" ) {
    library ( 'enrichplot' )
    library ( 'ggplot2' )
    library ( 'Cairo' )
    Cairo( file = filename , type="pdf" , width=5000 , height=3000 , units="px", dpi=200)
    d = dotplot( ck ) + theme(
        axis.text.x = element_text(angle=45, hjust=1) ,
        axis.text.y = element_text( angle=0, hjust=1) ,
        panel.grid = element_blank()
        )
    print ( d )
}

universal_cutoff <- 0.05

O  <- doClusterProfilingRehash(cutoff=universal_cutoff,use=".q")
print(O)
#
# RETRIEVE ALL FOLDCHANGES
foldchanges <- translated_values_from_file( value="B.T2D.NGT.fc" )
ck <- cluster_enrichment( O )

draw_dot_plot( ck )

namespace <- namespace_from_file()
#
# THIS IS THE TRADITIONAL ENRICHMENT 
# ANALYSIS QUANTIFICATION USING GO
# THESE RESULTS ARE ALMOST IDENTICALL TO
# WHAT CAN BE OBTAINED FROM THE PYTHON
# IMPETUOUS PACKAGE
#
value = "B.M.q"
ego <- enrichGO(gene          = O[[value]]   ,
                universe      = namespace    ,
                OrgDb         = org.Hs.eg.db ,
                ont           = "CC" ,
                pAdjustMethod = "BH" ,
                pvalueCutoff  = universal_cutoff ,
                qvalueCutoff  = universal_cutoff )
#
# GRAB THE PERTINENT FOLDCHANGES AND
# FEED THEM TO THE CNET
select_fc_values = foldchanges[c(O[[value]]),]
names(select_fc_values) <- c(O[[value]])
#
# COMPARTMENTS
create_cnet_plot( ego  , foldchange =  select_fc_values )
draw_upset_plot( ego )
#
# FUNCTION
ego_mf <- enrichGO(gene          = O[[value]]   ,
                universe      = namespace    ,
                OrgDb         = org.Hs.eg.db ,
                ont           = "MF" ,
                pAdjustMethod = "BH" ,
                pvalueCutoff  = universal_cutoff ,
                qvalueCutoff  = universal_cutoff )
create_cnet_plot( ego_mf  , foldchange =  select_fc_values ,  filename = "../results/cnetPathways_MF.pdf" )
draw_upset_plot( ego_mf  , filename = "../results/upsetPathways_MF.pdf" )
#
# Process
ego_bp <- enrichGO(gene          = O[[value]]   ,
                universe      = namespace    ,
                OrgDb         = org.Hs.eg.db ,
                ont           = "BP" ,
                pAdjustMethod = "BH" ,
                pvalueCutoff  = universal_cutoff ,
                qvalueCutoff  = universal_cutoff )
create_cnet_plot( ego_bp  , foldchange =  select_fc_values ,  filename = "../results/cnetPathways_BP.pdf" )
draw_upset_plot( ego_bp  , filename = "../results/upsetPathways_BP.pdf" )
