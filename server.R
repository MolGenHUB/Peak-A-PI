source("global.R")

shinyServer(function(input, output) {
  # dir <- reactive({ input$wd })
  
  output$uibam <- renderUI({
    filenames <- list.files(path = "dat", pattern = ".bam$")
    selectInput(inputId = "bam", label = "Choose Dataset", filenames)
  })
  
  sbf <- reactive({
    paste("dat", input$bam, sep = "/")
  })
  
  output$uichr <- renderUI({
    if (is.null(input$bam)) 
      return(NULL)
    seqnames <- getSeqnamesBam(sbf())
    selectInput(inputId = "chr", label = "Choose Chromosome", seqnames)
  })
  
  cvg <- reactive({
    chrcvg <- getCoverage(sbf(), input$chr)
    chrcvg
  })
  
  gr <- reactive({
    gr <- peakapi(dat = cvg(), chr = input$chr, bg = input$bg, min.se = input$se, 
                     min.ft = input$ft, Edges = input$end)
    ### Add codes to remove duplicate records based on sequences.
    names(gr) <- paste("FC", 1:length(gr), sep = "")
    gr
  })
  
  bedtable <- reactive({
    gr <- gr()
    table <- data.frame(chr = input$chr, start = start(gr), end = end(gr), name = names(gr), 
                        score = "*", strand = strand(gr))
    table
  })
  
  output$bed <- renderDataTable({
    bedtable()
  })
  
  output$dl <- downloadHandler(filename = function() {
    paste(format(Sys.time(), "%Y%m%d%I%M"), "bed", sep = ".")
  }, content = function(file) {
    # Howto add parameters as comments in the header of generated file?
    write.table(bedtable(), file, sep = "\t", quote = F, row.names = F, col.names = F)
  })
})
