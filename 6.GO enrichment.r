#signature: your gene module file with two columns, Symbol and module
#genecol: the column index of gene symbols in your gene module file
#ontology: the pathway annotation file
#background: human_coding_gene_list.txt



mygoenrich<-function(signature,genecol=1,ontology,background){
	output=data.frame()
	mod=unique(signature$module)
	bggene=unique(c(background$symbol,background$prev_symbol));bggene=bggene[!is.na(bggene)]
	nbackground=dim(background)[1]
	for (item in mod){
		genelist=unique(signature[signature$module==item,genecol])
		genelist=genelist[genelist %in% bggene]
		sigsize=length(genelist)
		result=sapply(ontology$GeneSymbol,function(x){
			ontogene=unlist(strsplit(x,"; "))
			ontogene_size=length(ontogene)
			ontogene=ontogene[ontogene %in% bggene]
			ontogene_tsize=length(ontogene)
			overlap_gene=paste(intersect(genelist,ontogene),collapse=";")
			overlap_size=length(intersect(genelist,ontogene))
			fold=overlap_size/(ontogene_tsize*sigsize/nbackground)
			pv=phyper(overlap_size,ontogene_tsize,nbackground-ontogene_tsize,sigsize,lower.tail=F)
			return(c(fold,pv,sigsize,nbackground,ontogene_size,ontogene_tsize,overlap_size,overlap_gene))
		},USE.NAMES = F)
		result=t(result)
		result=data.frame(result)
		names(result)=c("Fold","Pvalue","Sig.Size","Background.Size","Ontology.Size","Ontology.TrueSize","Overlap","Overlap.Sig")
		result$MODULE=item;result$Ontology=ontology$GeneCategory
		result$Pvalue=as.numeric(as.character(result$Pvalue))
		result=result[order(result$Pvalue),]
		result$FDR=p.adjust(result$Pvalue,"BH")
		result=result[,c(9,10,1,2,11,3:8)]
		output=rbind(output,result)
	}
	return(output)
}

signature <- read.csv("modules_for_GO_mapped.csv", header = T,stringsAsFactors = F) #read the csv of genes mapped to human genes
#colnames(signature) <- c("Symbol", "module")
ontology <- read.delim('msigdb.pathway+gobp.txt', header = T, sep='\t', stringsAsFactors = F)
background <- read.delim('human_coding_gene_list.txt', header = T, sep = '\t', fill = T, stringsAsFactors = F)
go_modules <- mygoenrich(signature = signature, genecol = 1, ontology = ontology, background = background)
write.csv(go_modules, "GO_modules_ing.csv")
