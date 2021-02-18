#load in our files containing annotations

#read in annotations from different source

#bioMart
bioMart <- read.csv("/Users/pdubelko/Desktop/TBDM/Annotations/illumina_biomart.csv", header = TRUE, stringsAsFactors = FALSE)
head(bioMart,20)

#platform specific
company <- read.table("/Users/pdubelko/Desktop/Annotations.txt",
                      header = TRUE, sep = "\t")
company <- company[c(1,5)]

#geo
geo <-read.table("/Users/pdubelko/Desktop/TBDM/Annotations/geo_illumina_anno.txt", header = TRUE, sep = "\t")
head(geo)

#only for illumnia - used study annotation (I believe TBDM)
study <- read.table("/Users/pdubelko/Desktop/illumina_annotation_study.txt", header = TRUE, sep = "\t")


bioMart.df <- data.frame(probeID = bioMart$ProbeID, 
                         Symbol = bioMart$SYMBOL, stringsAsFactors = FALSE)

#clean up bioMart annotation
bioMart.df[1,1] == bioMart.df[2,1]
for (i in 2:nrow(bioMart.df)-1){
  if (bioMart.df[i,1] == bioMart.df[i+1,1]){
    if (bioMart.df[i,2] != bioMart.df[i+1,2]){
      bioMart.df[i,2] = paste(bioMart.df[i,2], "_" , bioMart.df[i+1,2])
    }
    else{
      bioMart.df <- bioMart.df[-i+1]
    }
  }
}
write.table(bioMart.df, file = "biomart_illumina", sep = "\t")

company.df <- data.frame(probeID = company$Transcript, 
                         Symbol = as.character(company$ILMN_Gene), stringsAsFactors = FALSE)

geo.df <- data.frame(probeID = geo$Transcript, Symbol = as.character(geo$Symbol), stringsAsFactors = FALSE)


study.df <- data.frame(probeID = study$ProbeID, Symbol = as.character(study$Study), stringsAsFactors = FALSE)

#MERGE all dataframes
#merge first two
bioCompany <- merge(bioMart.df,company.df, by = "probeID", all.x = TRUE, all.y = TRUE, stringsAsFactors = FALSE)
head(bioCompany)
names(bioCompany) <- c("probeID", "biomart", "company")
#add in third to data.frame already merged
complete <- merge(bioCompany, geo.df, by = "probeID", all.x = TRUE, all.y = TRUE, stringsAsFactors = FALSE)


final <- merge(complete, study.df, by = "probeID", all.x = TRUE, all.y = TRUE, stringsAsFactors = FALSE)
#clarify where symbols came from
names(final) <- c("probeID", "bioMart", "company", "GEO","study")

write.table(final, "all_annotations_illumina.txt", sep = "\t")
