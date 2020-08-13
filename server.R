#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
    
    
    
    # Make mqpar for data input -----------------------------------------------
    
    
    expDes <- reactive({
        req(input$mqpar_file)
        mqpar_xml <- xmlTreeParse(input$mqpar_file$datapath)
        topMqpar <- xmlRoot(mqpar_xml)
        topMqpar <- xmlApply(topMqpar, function(x) xmlSApply(x, xmlValue))
        mqpar_df <- data.frame(t(topMqpar), row.names = NULL)
        
        expDesign <- data.frame(mqpar_df$filePaths)
        expDesign <- data.frame(splitstackshape::cSplit(indt = expDesign,splitCols = "filePaths",sep = "\\"))
        expDesign <- data.frame(filenames=expDesign[,dim(expDesign)[2]])
        expDesign$experiments <- data.frame(experiments=mqpar_df$experiments,stringsAsFactors = F)[,1]
        expDesign$experiments <- as.character(expDesign$experiments)
        expDesign$Volume_H <- ""
        expDesign$Volume_L <- ""
        if (input$multiplicity == "triple"){
            expDesign$Volume_M <- ""   
        }
        return(expDesign)
    })
    
    
    # Hot-table changes -------------------------------------------------------
    
    MyChanges <- reactive({
        if(is.null(input$hotable1)){return(expDes())}
        else if(!identical(expDes(),input$hotable1)){
            # hot.to.df function will convert your updated table into the dataframe
            as.data.frame(hot.to.df(input$hotable1))
        }
    })
    output$hotable1 <- renderHotable({MyChanges()}, readOnly = F)
    output$tbl = DT::renderDataTable(MyChanges(), options = list(pageLength = 50))
    
    
    
    # Upload peptides.txt ---------------------------------------------------------
    
    output$pepfile <- renderDataTable({
        req(input$peptides_file)
        
        df <- read.delim(file = input$peptides_file$datapath,
                         header = T,
                         sep = "\t",
                         stringsAsFactors = F)
        peptidesTable <<- df
        return(df)
    })
    
    
    
    # Main calculations -------------------------------------------------------
    
    mixcalcTableDL <- reactive({
        req(input$calculate)
        InitialMix <- MyChanges()
        experiments <- InitialMix$experiments
        
        
        if (input$multiplicity == "double"){
            InitialMix$Proportion_H <- ""
            InitialMix$Proportion_L <- ""
            InitialMix$Corrected_volume_H <- ""
            InitialMix$Corrected_volume_L <- ""
            InitialMix$Volume_H <- as.numeric(InitialMix$Volume_H)
            InitialMix$Volume_L <- as.numeric(InitialMix$Volume_L)
            peptidesTable_intensities <- peptidesTable[,c(grep("Intensity.H.",colnames(peptidesTable)), 
                                                          grep("Intensity.L.",colnames(peptidesTable)))]
            peptidesTable_intensities_sums <- data.frame(t(colSums(peptidesTable_intensities)))
            
            for (i in experiments){
                temp <- peptidesTable_intensities_sums[,grep(i, colnames(peptidesTable_intensities_sums))]
                temp2 <- temp
                temp2[1,1] <- temp[1,1]/sum(temp[1,1],temp[1,2]) #heavy
                temp2[1,2] <- temp[1,2]/sum(temp[1,1],temp[1,2]) #light
                InitialMix[grep(i,InitialMix$experiments),"Proportion_H"] <- round((temp2[1,1]*100),2)
                InitialMix[grep(i,InitialMix$experiments),"Proportion_L"] <- round((temp2[1,2]*100),2)
                
                InitialMix[grep(i,InitialMix$experiments),"Corrected_volume_H"] <- round(InitialMix[grep(i,InitialMix$experiments),"Volume_H"]/(temp2[1,1]/0.50),2)
                InitialMix[grep(i,InitialMix$experiments),"Corrected_volume_L"] <- round(InitialMix[grep(i,InitialMix$experiments),"Volume_L"]/(temp2[1,2]/0.50),2)
            }
        }
        
        if (input$multiplicity == "triple"){
            InitialMix$Proportion_H <- ""
            InitialMix$Proportion_L <- ""
            InitialMix$Proportion_M <- ""
            InitialMix$Corrected_volume_H <- ""
            InitialMix$Corrected_volume_L <- ""
            InitialMix$Corrected_volume_M <- ""
            InitialMix$Volume_H <- as.numeric(InitialMix$Volume_H)
            InitialMix$Volume_L <- as.numeric(InitialMix$Volume_L)
            InitialMix$Volume_M <- as.numeric(InitialMix$Volume_M)
            
            peptidesTable_intensities <- peptidesTable[,c(grep("Intensity.H.",colnames(peptidesTable)), 
                                                          grep("Intensity.L.",colnames(peptidesTable)),
                                                          grep("Intensity.M[.]",colnames(peptidesTable)))]
            
            peptidesTable_intensities_sums <- data.frame(t(colSums(peptidesTable_intensities)))
            
            for (i in experiments){
                temp <- peptidesTable_intensities_sums[,grep(i, colnames(peptidesTable_intensities_sums))]
                temp2 <- temp
                temp2[1,1] <- temp[1,1]/sum(temp[1,1],temp[1,2],temp[1,3]) #heavy
                temp2[1,2] <- temp[1,2]/sum(temp[1,1],temp[1,2],temp[1,3]) #light
                temp2[1,3] <- temp[1,3]/sum(temp[1,1],temp[1,2],temp[1,3]) #med
                
                InitialMix[grep(i,InitialMix$experiments),"Proportion_H"] <- round((temp2[1,1]*100),2)
                InitialMix[grep(i,InitialMix$experiments),"Proportion_L"] <- round((temp2[1,2]*100),2)
                InitialMix[grep(i,InitialMix$experiments),"Proportion_M"] <- round((temp2[1,3]*100),2)
                
                InitialMix[grep(i,InitialMix$experiments),"Corrected_volume_H"] <- round(InitialMix[grep(i,InitialMix$experiments),"Volume_H"]/(temp2[1,1]/0.3333),2)
                InitialMix[grep(i,InitialMix$experiments),"Corrected_volume_L"] <- round(InitialMix[grep(i,InitialMix$experiments),"Volume_L"]/(temp2[1,2]/0.3333),2)
                InitialMix[grep(i,InitialMix$experiments),"Corrected_volume_M"] <- round(InitialMix[grep(i,InitialMix$experiments),"Volume_M"]/(temp2[1,3]/0.3333),2)
            }
        }
        return(InitialMix)
    })
    
    # Mixing Table Render -----------------------------------------------------------
    
    output$mixcalcTable <- DT::renderDataTable({
        mixTable <- mixcalcTableDL()
        # return(mixTable)
        # DT::datatable(mixTable, options = list(pageLength = 50)) %>% DT::formatStyle(c("Corrected_volume_H","Corrected_volume_L"), backgroundColor = '#E95420')
        DT::datatable(mixTable, options = list(pageLength = 50)) %>% 
            DT::formatStyle(grep("Corrected_volume",colnames(mixTable)), color = "#E95420", fontWeight = "bold")
    })
    
    # Mixing download -----------------------------------------------------------
    
    output$download <- downloadHandler(
        filename = function() {
            paste0("MixingCheckCalc.txt")
        },
        content = function(file) {
            write.table(mixcalcTableDL(), file, sep = "\t", row.names = FALSE)
        }
    )
    
    # Help Page --------------------------------------------------------------
    
    output$helppage <- renderUI({
        HTML("<font color='#E95420' size=8>",
             "Instructions for use",
             "</font>",
             "<br/>",
             "<font color='#000000' size=5>",
             "- Ensure that your experiment name has no spaces when setting up MaxQuant","<br/>",
             "- Select SILAC/Dimethyl Experiment type (i.e.: double, single","<br/>",
             "- Select <strong>peptides.txt</strong> file. This should be in the txt directory of your maxquant run,","<br/>",
             "- Select the <strong>mqpar.xml</strong> file. This should be in the same directory as your raw files","<br/>",
             "- The table on the left in the <strong>Input</strong> tab is editable. Input the volumes used for the Heavy, Light, and/or Medium constituent of each sample into the empty <strong>Volume</strong> columns. You may copy and paste values","<br/>",
             "- Press the <strong>Calculate Mixing</strong> button","<br/>",
             "- Select the <strong>Mixing Calculation</strong> tab to view the results. The new volumes for correct mixing are highlighted","<br/>",
             "- You may press the <strong>Download Mixing Calculations</strong> button at the bottom of the page to download this table in txt format (tab-separated)","<br/>",
             "</font>")
        
    })
    
    
    # About Page --------------------------------------------------------------
    
    output$aboutpage <- renderUI({
        HTML("<font color='#E95420' size=5>","Author:","</font>",
             "<font color='#000000' size=5>"," Brandon Dean Murugan","</font>",
             "<br/>",
             "<font color='#E95420' size=5>","Email:","</font>",
             "<font color='#000000' size=5>","brandonmurugan@gmail.com","</font>")
        
    })
    
})