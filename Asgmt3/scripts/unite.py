import csv

Up_file = open("./Up_genes.list")
Down_file = open("./Down_genes.list")

Up_dea_genes = open("./Up_dea_genes.csv","w",newline='')
Down_dea_genes = open("./Down_dea_genes.csv","w",newline='')

Pf = open("./Promoter.annotated")
Pf_list = Pf.readlines()

writer = csv.writer(Up_dea_genes)
writer.writerow(['Up-dea-gene'])
for upGene in Up_file.readlines():
    gene = upGene.replace('\n', r'\n')[:-2]
    #print(gene)
    for annotate_info in Pf_list:
        #print(gene)
        #print(annotate_info)
        if gene in annotate_info and gene != '':
            writer.writerow([gene])
            break

writer = csv.writer(Down_dea_genes)
writer.writerow(['Down-dea-gene'])
for downGene in Down_file.readlines():
    gene = downGene.replace('\n', r'\n')[:-2]
    #print(gene)
    for annotate_info in Pf_list:
        #print(gene)
        #print(annotate_info)
        if gene in annotate_info and gene != '':
            writer.writerow([gene])
            break


Up_file.close()
Down_file.close()
Up_dea_genes.close()
Down_dea_genes.close()
