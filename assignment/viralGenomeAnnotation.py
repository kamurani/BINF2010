# Import modules etc 
import sys
#import os
import re 
from Bio import Entrez
from Bio import SeqIO
from itertools import repeat
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

## Students may use some of these (mostly for the graphics)
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW

# Set up Entrez access (need to do this, need to be mindful of database/url access)
Entrez.email = "s.ballouz@unsw.edu.au"
Entrez.tool = "viral.py"


# Get arguments from commandline 
args = sys.argv
if(len(args) == 1 ):
    print("Missing arguments, please specify accession/taxid")
    exit()

# Extract "first" argument
input = args[1]

# Extract features and write to GTF file   
def extract_gtf(query_feature,query_seq, query):
    dict_codon={}
    output_handle_pro=open(query+"_protein.fa", "w")
    output_handle=open(query+".gtf", "w")
    for feature in query_feature:
        #print(feature)
        feat_int = feature["GBFeature_intervals"][0]
        if len(feat_int) == 2:
            fromi=int(feat_int["GBInterval_point"])
            toi=fromi
        else:
            fromi=int(feat_int["GBInterval_from"])
            toi=int(feat_int["GBInterval_to"])
        if fromi > toi:
                fromi2=toi
                toi=fromi
                fromi=fromi2
                strandi= "-"
        else:
            strandi= "+"    
        feat_key = feature["GBFeature_key"]
        ### Is there a better way to parse/check this? Breaks for CDS? 
        if feat_key == "gene" or feat_key == "CDS" or feat_key=="regulatory": 
            feat_quals = feature["GBFeature_quals"][0]['GBQualifier_value']
            feat_len = len(feature["GBFeature_quals"])
            if feat_key == "CDS":
                #print("?")
                gene_id = feature["GBFeature_quals"][feat_len-2]['GBQualifier_value']
                pro_seq = feature["GBFeature_quals"][feat_len-1]['GBQualifier_value']
                cds_seq=query_seq[fromi:toi]    
                GC_content(cds_seq)
                #print(dict_codon)
                dict_codon = merge_dict(dict_codon , codon_use(cds_seq,query))
                output_handle_pro.write(">" + feat_quals + "\n")
                output_handle_pro.write(pro_seq + "\n")
            else: 
                gene_id = feature["GBFeature_quals"][feat_len-1]['GBQualifier_value']
            #if feat_key == "regulatory":
            #    print(feature)
            record_line = query + "\t" + feat_key + "\t" + str(fromi) + "\t" + str(toi) + "\t.\t" + strandi + "\t.\tgene_id \"" + gene_id + "\"; gene_name \"" + feat_quals + "\"\n"    
            #print(record_line)            
            output_handle.write(record_line)
    #print(dict_codon)
    output_handle_pro.close()
    output_handle.close()
    ## print out         
    output_handle=open((query +"_codon_usage.txt"), "w")
    for k, v in dict_codon.items():
        output_handle.write(str(k) + '\t'+ str(v) + '\n')
    #output_handle.write(dict_codon)
    output_handle.close() 



# Extract features and draw plot with basic plotting
def plot_genome_gff_std(query_feature,query_seq,query_id):
    print("Plotting genome")
    query_seq=query_seq.lower()
    nmax=len(query_feature)
    viridis_colors= mpl.colormaps['viridis'].resampled(i)
    i=0
    features={}
    for feature in query_feature:
    #print(feature)
        feat_int = feature["GBFeature_intervals"][0]
        feat_key = feature["GBFeature_key"]
        ### Is there a better way to parse/check this? Breaks for CDS? 
        if feat_key == "gene":
            feat_quals = feature["GBFeature_quals"][0]['GBQualifier_value']
            x1=int(feat_int["GBInterval_from"])
            x2=int(feat_int["GBInterval_to"])
            point1=(x1,i-1)
            point2=(x1,i)
            point3=(x2,i)
            point4=(x2,i-1)
            points = (point1, point2, point3, point4,point1)
            features[feat_quals] = points
            i = i + 1
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(15,5), sharex=True, gridspec_kw={"height_ratios": [4, 1]}
    )
    x = np.arange(len(query_seq) - 50)
    y=list(repeat(i,len(x)))
    #ax1.plot(x,y)
    j=0
    ax1.set_ylim(bottom=-1, top=i)
    for feat in features.keys():
        coords=(features[feat])
        xs, ys = zip(*coords) 
        ax1.fill(xs,ys,color=viridis_colors(j), alpha=0.3)
        ax1.plot(xs,ys,color=viridis_colors(j), alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
        ax1.annotate(feat,(min(xs)+1,ys[0]))
        j +=1
    ax1.set_yticklabels([])
    ax1.set_xlabel('bp')
    gc = lambda s: 100.0 * len([c for c in s if c in "gc"]) / 50
    xx = np.arange(len(query_seq) - 50)
    yy = [gc(query_seq[x : x + 50]) for x in xx]
    ax2.fill_between(xx + 25, yy, alpha=0.3)
    ax2.set_ylim(bottom=0)
    ax2.set_ylabel("GC(%)")   
    fig.figure.savefig(query_id + ".png")     

# Extract features and draw plot with basic plotting v2
def plot_genome_gff_std2(query_feature,query_seq,query_id):
    print("Plotting genome")
    query_seq=query_seq.lower()
    bw=len(query_seq)*0.005
    i=0
    features={}
    for feature in query_feature:
    #print(feature)
        feat_int = feature["GBFeature_intervals"][0]
        feat_key = feature["GBFeature_key"]
        ### Is there a better way to parse/check this? Breaks for CDS? 
        if feat_key == "gene":
            feat_quals = feature["GBFeature_quals"][0]['GBQualifier_value']
            x1=int(feat_int["GBInterval_from"])
            x2=int(feat_int["GBInterval_to"])
            point1=(x1,i-1)
            point2=(x1,i)
            point3=(x2,i)
            point4=(x2,i-1)
            if x1 < x2: 
                point3a = (x2-bw,i)
                point4a = (x2-bw,i-1)
                point3b = (x2,i-0.5)
            if x1 >= x2: 
                point3a = (x2+bw,i)
                point4a = (x2+bw,i-1)
                point3b = (x2,i-0.5)
            
            # rectangle
            points = (point1, point2, point3, point4,point1)
            # arrow         
            points = (point1, point2, point3a, point3b, point4a,point1)
            
            features[feat_quals] = points
            i = i + 1
    viridis_colors= mpl.colormaps['viridis'].resampled(i)
    x = np.arange(len(query_seq) - 50)
    y = list(repeat(i,len(x)))
    fig, ax = plt.subplots(figsize =(15, 5))
    j=0
    for feat in features.keys():
        coords=(features[feat])
        xs, ys = zip(*coords) 
        ax.fill(xs,ys,color=viridis_colors(j), alpha=0.3)
        ax.plot(xs,ys,color="black", linewidth=1, solid_capstyle='round')
        ax.annotate(feat,(min(xs)+bw,ys[0]+0.25))
        j +=1
    plt.suptitle(query_id) 
    ax.set_yticklabels([])
    ax.set_xlabel('bp')
    fig.figure.savefig(query_id + "_v2.png")     

# Extract features and draw plot using dna_features_viewer
def plot_genome_gff(query_feature,query_seq,query_id):
    print("Plotting genome")
    query_seq=query_seq.lower()
    nmax=len(query_feature)
    i=0
    viridis_colors = mpl.colormaps['magma'].resampled(nmax)
    features=[]
    for feature in query_feature:
        #print(feature)
        feat_int = feature["GBFeature_intervals"][0]
        feat_key = feature["GBFeature_key"]
        ### Is there a better way to parse/check this? Breaks for CDS? 
        if feat_key == "gene":
            feat_quals = feature["GBFeature_quals"][0]['GBQualifier_value']
            record_line = GraphicFeature(start=int(feat_int["GBInterval_from"]), end=int(feat_int["GBInterval_to"]), strand=+1, color=viridis_colors(i), label=feat_quals)
            features.append(record_line)
        if feat_key == "regulatory":
            feat_quals = feature["GBFeature_quals"][0]['GBQualifier_value']
            record_line = GraphicFeature(start=int(feat_int["GBInterval_from"]), end=int(feat_int["GBInterval_to"]), strand=+1, color="green", label=feat_quals)
            features.append(record_line)    
        i = i + 1 
    graphic_record = GraphicRecord(sequence_length=len(query_seq), features=features)
    #graphic_record.plot(figure_width=15)  
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(15,5), sharex=True, gridspec_kw={"height_ratios": [4, 1]}
    )
    graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4)
    gc = lambda s: 100.0 * len([c for c in s if c in "gc"]) / 50
    xx = np.arange(len(query_seq) - 50)
    yy = [gc(query_seq[x : x + 50]) for x in xx]
    ax2.fill_between(xx + 25, yy, alpha=0.3)
    ax2.set_ylim(bottom=0)
    ax2.set_ylabel("GC(%)")
    fig.figure.savefig(query_id + ".png")        

# Draw plot using Bio.Graphics and a genbank record from SeqIO
def plot_genome(record):
    gd_diagram = GenomeDiagram.Diagram(record.id)
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()
    for feature in record.features:
        if feature.type != "gene":
            #Exclude this feature
            continue
        if len(gd_feature_set) % 2 == 0:
            color = colors.purple
        else:
            color = colors.lightcoral
        gd_feature_set.add_feature(feature, sigil="ARROW",
                                color=color, label=True,
                                label_size = 14, label_angle=0)        
    gd_diagram.draw(format="linear", pagesize='A4', fragments=4,
                    start=0, end=len(record))
    gd_diagram.write(record.id+"_linear_nice.pdf", "PDF")
    gd_diagram.write(record.id+"_linear_nice.eps", "EPS")
    gd_diagram.write(record.id+"_linear_nice.svg", "SVG")
    gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm),
                    start=0, end=len(record), circle_core = 0.5)
    gd_diagram.write(record.id+"_circular_nice.pdf", "PDF")
    gd_diagram.write(record.id+"_circular_nice.eps", "EPS")
    gd_diagram.write(record.id+"_circular_nice.svg", "SVG")

# Get codon table
def codon_use(query_seq,query) -> dict:
    dict_codon={}
    print("Count codons")
    for i in range(len(query_seq)):
        codon=str(query_seq[i:i+3])
        if len(codon) < 3: 
            break
        if codon not in dict_codon.keys():
            dict_codon[codon] = 1
        else:
            dict_codon[codon] += 1
    return(dict_codon)

# Merge dictionaries 
def merge_dict(dict1, dict2) -> dict:
    for key in dict1.keys():
        if key not in dict2.keys():
            dict2[key] = dict1[key]
        else:
            dict2[key] += dict1[key]
    return(dict2)




# Get GC content 
def GC_content(seq):
    #print(len(seq))
    gc = lambda s: 100.0 * len([c for c in s if c in "gc"])/len(s) 
    gc_total = round(gc(seq),2) 
    print("GC%:" + str(gc_total))

# Get first accession result from taxID search
def taxid2accession(txid: int, db="nuccore") -> str:
    print("Convert taxID to nucleic acid accession ID")
    handle = Entrez.esearch(db=db, term="txid"+str(txid)+"[Organism] AND refseq[filter] AND complete[title]")
    record = Entrez.read(handle)
    gi = record["IdList"][0]
    handle = Entrez.esummary(db=db, id=gi, retmode="xml")
    result = Entrez.read(handle)
    acc = result[0]['AccessionVersion']
    return str(acc)

# TODO - Get all accession result from taxID search - for multiple results, currently have not specified which to take, first is best, will adjust the spec
def taxid2accession_all(txid: int, db="nuccore") -> str:
    print("Convert taxID to nucleic acid accession ID")
    handle = Entrez.esearch(db=db, term="txid"+str(txid)+"[Organism] AND refseq[filter]")
    handle = Entrez.esearch(db=db, term="txid"+str(txid)+"[Organism] AND refseq[filter] AND complete[title]")
    record = Entrez.read(handle)
    print(len(record["IdList"]))
    #for gi in record["IdList"]:
    #    print(gi)
    #handle = Entrez.esummary(db=db, id=gi, retmode="xml")
    #result = Entrez.read(handle)
    #acc = result[0]['AccessionVersion']
    #print(acc)
#    return str(acc)



# Query database for GenBank entry, using nucleic acid accession ID 
def query_nc_id(query):
    print("Querying GenBank")
    handle = Entrez.esearch(db="nuccore", term=query)
    record = Entrez.read(handle)
    num_results = int(record['Count'])
    if(num_results > 0):
        query_ids = record['IdList']
    for query_id in query_ids:
        ### To print to FASTA format using the SeqIO module
        print("Writing FASTA file") 
        handle = Entrez.efetch(db="nuccore", id=query_id, rettype="gb", retmode="text")
        seq_record = SeqIO.read(handle, "gb")
        output_handle=open((query +".fasta"), "w")
        SeqIO.write(seq_record, output_handle, 'fasta')
        output_handle.close() 
        #plot_genome(seq_record)
        ## To extract GTF file   
        print("Getting features")        
        handle = Entrez.efetch(db="nuccore", id=query_id, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        query_seq = record[0]["GBSeq_sequence"]
        query_acc = record[0]["GBSeq_accession-version"]
        query_source = record[0]["GBSeq_source"]
        query_tx = record[0]["GBSeq_taxonomy"]  
        query_org = record[0]["GBSeq_organism"]  
        # Write GTF file  
        query_feature = record[0]["GBSeq_feature-table"] ## output as gff/gtf file 
        extract_gtf(query_feature,query_seq,query)
        plot_genome_gff(query_feature,query_seq,query)
        plot_genome_gff_std2(query_feature,query_seq,query)
        
 

# Query database for GenBank entry, using taxid (convert to nucc access first) 
def query_tax_id(query):
    acc=taxid2accession(query,"nuccore")
    query_nc_id(acc)

    
## Input checks 
## If matches accession ID 
if re.search("^\D", input):
    print("Nucleic acid accession ID: " + input)
    query_nc_id(input)

## If matches a TAXID 
elif re.search("^\d", input):
    print("Taxonomy ID: " + input)
    query_tax_id(input)

else:
    print("Missing")

# Example run 
# Option 1
## python viralGenomeAnnotation.py NC_045512.2  

# Option 2  
## python viralGenomeAnnotation.py 10376 
