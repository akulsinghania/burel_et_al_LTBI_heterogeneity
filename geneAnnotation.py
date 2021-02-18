import csv
# readIn function cleans collected annotation files and outputs dictionary to be passed into the condense function. This entire function takes into account missing genes from each data source and formatting issues to create a uniform dictionary. Order of file is BioMart (col 1), GEO (col2), Platform annotations (col3)

def readIn(file):
    text = open(file, "r")
    probes = {}
    for line in text:
        line = line.strip()
        info = line.split("\t")
	#account for probes that do not have associated genes from Biomart/GEO/Platform
	#do not want non genes included in final dictionary
        if "-" in info[1:]:
            info.remove("-")
        if "---" in info[1:]:
            info.remove("---")

	#info[0] is probe
        for item in info[1:]:
	    #account for multiple gene annotations from source, split by _
            if item.find("_") != -1:
                new = item.split("_")
                for val in new:
                    val = val.strip()
                    info.append(val)
                info.remove(item)
        for item in info[1:]:
	    #account for multiple gene annotations from source, split by ///
            if item.find("///") != -1:
                new = item.split("///")
                for val in new:
                    val = val.strip()
                    info.append(val)
                info.remove(item)
	#if probe already in dictionary, add new entries
        if info[0] in probes:
            for item in info[1:]:
                probes[info[0]].append(item)
	#first entry, append list of genes
        else:
            probes[info[0]] = info[1:]
    return probes

#Condense function is used to find our gene we have the most confidence in. We pass in a dictionary, each key is a probe and the values are lists of genes. Because of how they were read in (BioMart > GEO > Platform), we give higher confidence to genes with lower indexes in the case of same counts. This function loops through all probes and goes through the gene list, counts each item and see's if it appears more than half the time or if not how many times it appears and if it appears the most, or in the case of appearing the same amount as other genes, what source it came from. The output is our final annotation for a single platform.

def condense(probeList):
    #for each probe, loop through
    for key in probeList:
        values = probeList[key]
        count = 0
        mostFreq = None
	#look at each gene associated with current probe
        for item in values:
            if values.count(item) >= count:
                #is this gene the most frequent
                if mostFreq:
                    if values.count(item) == count:
			#higher confidence genes have more confidence, if count is =, default to lower index gene.
                        if values.index(item) < values.index(mostFreq):
                            mostFreq = item
                    else:
                        mostFreq = item
                else:
                    mostFreq = item
		#if statement if mostFreq > len(values)/2 you have found consensus gene. 
                if values.count(mostFreq) > len(values)/2:
                        probeList[key] = mostFreq
                        break
                count = values.count(mostFreq)
            probeList[key] = mostFreq
    return probeList


probes = readIn("all_annotations_illumina.txt")

condensed = condense(probes)

w = csv.writer(open("annotations_illumina.csv", "w"))
for key, val in condensed.items():
    w.writerow([key,val])
