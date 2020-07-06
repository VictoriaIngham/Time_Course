correct_gene_names3<- function (red_green_data, gene_names_file, positions)
{
    gene_names_file = as.matrix(gene_names_file)
  
  for(i in 1:length(positions))
    {
	if(length(positions)==0)
	{
		break
	}

	found<-grep(red_green_data$genes$ProbeName[positions[i]],gene_names_file[,1])
	for(j in 1:length(found))
	{
		red_green_data$genes$SystematicName[positions[i]]<- gene_names_file[found[j],2]
	}
    }
    return(red_green_data)
}