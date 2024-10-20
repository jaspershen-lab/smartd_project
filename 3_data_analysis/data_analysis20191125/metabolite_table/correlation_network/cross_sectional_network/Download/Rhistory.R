mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec<-c("Erucamide","Phloroglucinocarboxaldehyde","D-Alanyl-D-alanine (D-Ala-D-Ala)","Imazamox","3-Acetoxypyridine","3'-O-methylcytidine","Pyruvate","1,4-Butynediol","Phytosphingosine","Isethionic acid sodium salt","Oxyquinoline","D-2-Phosphoglyceric acid","Suberic acid","O-Desmethylnaproxen","(Methylenecyclopropyl)acetyl-L-carnitine","Traumatic Acid","NICOTINATE","16-hydroxy hexadecanoic acid","N-Acetylmannosamine","3-Hydroxy-4-methoxybenzoic acid","D-Biotin","Alprazolam","Dodecanedioic acid","3-Hydroxyoctanoic acid","Brazilin","Senecioic acid","5-Fluorocytidine","Atipamezole","2-Hydroxyoctanoic acid","L-Methionine","Glycine","Methotrimeprazine","Dimethyl sulfoxide",".beta.-Carboline-1-propionic acid","DL-.alpha.-Aminopimelic acid","AMP (Adenosine 3'-Monophosphate Hydrate)","2,3-Dimethoxyphenol","Phenol","Confertifoline","5-Hydroxyindoleacetate","3-Hydroxypropionic acid (beta-lactic acid)","4-Hydroxy-L-proline","Glycodeoxycholic acid","(-)-Naringenin","4-Acetamidobutyric acid","3-Methylglutaric acid","PHIP","Crotonic acid","5-Hydroxymethyl-2-furancarboxylic acid","16.alpha.-Hydroxy-4-androstene-3,17-dione","Stearoyl-L-Carnitine","Beta-Alanine","3-Oxopentanoic acid","Estrone","Butylphosphonic acid","cis-4-Hydroxycyclohexanecarboxylic acid","4-Vinylphenol","Stearic acid","5-Hydroxyindoleacetic acid")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
mSet<-PlotPathSummary(mSet, "path_view_0_", "png", 72, width=NA)
mSet<-SaveTransformedData(mSet)
