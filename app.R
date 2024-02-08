# REQUIRED LIBRARIES #
library(plotly)
library(RCurl)
library(RColorBrewer)

data = read.table(file = "interpro_domains.tsv", stringsAsFactors=F, header=T, sep="\t",quote="")
geneset = unique(data[,c("Transcript.stable.ID","GeneName")])
names(geneset) <- c("transcript","gene")
genenames=sort(unique(geneset$gene))

# INVOKING THE SHINY APP BY R #
shinyApp(

  # DESIGNING THE USER INTERFACE #
  ui = fluidPage(

    # TITLE OF THE APPLICATION #
    titlePanel("AlphaMissense pathogenicity predictions"),
    br(),

    sidebarLayout(
      sidebarPanel(
          selectizeInput("genename", "Gene to plot", choices=NULL, multiple = FALSE, selected = NULL),
          br(),
          uiOutput("radiobuttons"),
          uiOutput("links"),
          width = 4),
      mainPanel(fluidRow(br(),column(12, plotlyOutput("singlegene_plot",height = "500px"))))
    )
  ),

  server = function(input, output, session)
  {
    updateSelectizeInput(session, 'genename', choices = genenames, server = TRUE)

    output$radiobuttons <- renderUI({
      genesToPlot <- input$genename
      if(genesToPlot != "")
      {
        transcript <- geneset$transcript[which(geneset$gene == genesToPlot)]
        radioButtons("trans", "Available transcripts:", choices=transcript)
      }
    })

    output$links <- renderUI({
      gene <- input$genename
      transcript <- geneset$transcript[which(geneset$gene == gene)]
      if(gene != "")
      {
      tagList(a("Ensembl transcript info",
                href=paste0("http://useast.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=",transcript)),
              br(),
              a("Uniprot protein search",
                href=paste0("https://www.uniprot.org/uniprotkb?query=",gene,"&facets=model_organism%3A9606")),
              br(),
              a("Interpro domain search",
                href=paste0("https://www.ebi.ac.uk/interpro/search/text/",gene,"/?page=1#table"))
      )
    }
    })
    

    output$singlegene_plot <- renderPlotly({
      transcriptToPlot <- input$trans
      geneToPlot <- input$genename
      #transcriptToPlot <- "ENST00000321117"
      if(length(transcriptToPlot) > 0)
      {
        df = as.data.frame(read.table(file = paste0("https://storage.googleapis.com/tcga_shiny/AlphaMissense_transcripts/",transcriptToPlot,".tsv"), header = FALSE, check.names=FALSE))
        names(df) <- c("Aminoacid","Score")
        df$Category <- "Ambiguous"
        df$Category[which(df$Score < 0.34)] <- "Benign"
        df$Category[which(df$Score > 0.56)] <- "Pathogenic"
        df$AApos = as.numeric(gsub("\\D", "", df$Aminoacid))

        #get the domains for this transcript
        domains = data[data$Transcript.stable.ID==input$trans,]
        domains = domains[,c("Interpro.start","Interpro.end","Interpro.ID","Interpro.Short.Description","Interpro.Description")]
        names(domains) = c("start","end","id","shortdesc","desc")
        level=0;
        vert.pos=c()
        domains = domains[with(domains, order(start, end)),]
        #walk through, keeping track of intersections, and bump the rectangles 
        #vertically up or down as needed to avoid overlaps
        pos = c(domains$start,domains$end)
        pos.type = c(rep("start",nrow(domains)),rep("end",nrow(domains)))
        pos.type=pos.type[order(pos)]
        for(i in pos.type){
            if(i=="start"){
              level = level + 0.1
              vert.pos=c(vert.pos,level)
            } else { #pos.type=="end"
              level = level - 0.1
            }
        }
        domains$btm = 1.0 + vert.pos
        domains$top = domains$btm + 0.1
        
        domains$tooltip = paste0("Interpro ID:   ",domains$id,"<br />",
                                 "Description:  ",domains$desc,"<br />",
                                 "Short Desc:  ",domains$shortdesc,"<br />",
                                 "Position:       ",domains$start,"-",domains$end)

        domcol = rep(brewer.pal(8,"Pastel2"),ceiling(max(table(data$Transcript.stable.ID))/8)) 
        g <- ggplot(data=df,aes(x=AApos, y=Score,
                                    text = paste0("Amino acid: ",Aminoacid,"<br />Score: ",
                                                  Score, "<br />Category: ", Category))) +
          geom_point(aes(color = Category)) +
          scale_color_manual(values = c("#fd8d3c", "#fed976", "#bd0026")) +
          ylab("Pathogenicity score") + theme_bw() +
          ggtitle(geneToPlot) +
          geom_rect(data=domains,aes(xmin=start, xmax=end, ymin=btm, ymax=top, 
                                      text=domains$tooltip), 
                    fill=domcol[1:nrow(domains)], 
                    inherit.aes=FALSE) +
          theme(text=element_text(size=12, family="avenir", face="bold"),
                axis.text=element_text(size=12, family="avenir", face="bold"),
                axis.text.x = element_text(size=12, family="avenir", face="bold"),
                axis.title.x = element_blank()) +
          ylim(0,max(domains$top))  
            
        ggplotly(g, tooltip="text") %>% layout(hoverlabel = list(align = "left")) 

      }
    })
  }
)
