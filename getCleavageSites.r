## This is a script to convert part of the refseq file into a .bed-like file centered on the poly-A cleavage site

source("/home/home/jeremy/Code/useful_R/useful_R.r")

get_command_args()
# inputFile, outputFile, Up, Down

input = read.delim(inputFile, stringsAsFactors=F)

## polyA site:
Ends = ifelse(input$strand=='+',input$txEnd,input$txStart)

Up = as.numeric(Up)
Down = as.numeric(Down)

# This is set up for 0-indexed data

# Determine positions of the upstream:
Up.left = Ends + ifelse(input$strand == '+', -(Up + 1), 0)
Up.right = Ends + ifelse(input$strand == '+', -1, Up)

# Determine positions of the downstream:
Down.left = Ends + ifelse(input$strand == '+', 0, -Down-1)
Down.right = Ends + ifelse(input$strand == '+', Down,-1)

nonNeg= (Down.left) > 0 & Down.right > 0 & Up.left > 0 & Up.right > 0

Output = rbind(
cbind(input$chr,Up.left,Up.right,paste(input$name,'Up',sep='_'),0,input$strand)[nonNeg,],
cbind(input$chr,Down.left,Down.right,paste(input$name,'Down',sep='_'),0,input$strand)[nonNeg,])

write.delim(Output,outputFile,col.names=F)
write.delim(Output[!duplicated(Output[,4]) ,],paste(outputFile,'noDups',sep='.'),col.names=F)

