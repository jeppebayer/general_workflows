template_file = "/faststorage/project/EcoGenetics/general_workflows/population_genetics/sfs/configuratons/EntNic/templates/EntNic_pooled_template.config.yaml"
pop_file = "/faststorage/project/EcoGenetics/general_workflows/population_genetics/sfs/configuratons/EntNic/templates/pop_id.txt"
pops = open(pop_file)
pop_list = [pop.strip("\n") for pop in pops]
sp = "EntNic"
align_path = "/faststorage/project/EcoGenetics/general_workflow_outputs/population_genetics/alignments/collembola/Entomobrya_nicoleti/grassland/"
for pop in pop_list:
    outfile = f"/faststorage/project/EcoGenetics/general_workflows/population_genetics/sfs/configuratons/EntNic/{sp}_pooled_{pop}.config.yaml"
    content = open(template_file)
    out = open(outfile,"w")
    for line in content:
        info = line.split(":")[0]
        if info == "population":
            out.write(f"population: {pop}\n")
        elif info == "population_bam_alignment":
            out.write(f"population_bam_alignment: {align_path}/{pop}/{pop}.filtered.bam\n")
        else:
            out.write(line)
    print(f"Done {pop}")



