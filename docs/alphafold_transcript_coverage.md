# Abstract

The AlphaFold/AlphaFold2 database from EMBL-EBI is comprised of a single protein structure per gene symbol. These are named based on UniProt accession numbers, with identical sequences to each AlphaFold structure; this allows for precise mapping of sequence-based information to the protein tertiary structures. While useful, there is enormous clinical value in building structures from alternative transcripts for each gene symbol, where we have knowledge that these are protein-coding and expressed. Often, alternate transcripts are of significant relevance in specific tissues, such as the brain or immune cells, and of specific interest to scRNA-Seq experiments. We explore the coverage of the protein-coding space in the human genome compared to the AlphaFold structures and show that 11.8% of actionable ClinVar variants are currently not represented and contained within alternate transcripts.  


# Steps


## Files downloaded

MANE Select transcripts  

-   MANE.GRCh38.v0.95.summary.txt.gz
-   <https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v0.95.summary.txt.gz>
-   Accessed 10/11/2021

AlphaFold2 Homo sapiens Structures  

-   UP000005640\_9606\_HUMAN.tar
-   <https://ftp.ebi.ac.uk/pub/databases/alphafold/UP000005640_9606_HUMAN.tar>
-   Accessed 9/11/2021

GENCODE protein-coding transcripts  

-   gencode.v39.pc\_transcripts.fa
-   <http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v39.pc_transcripts.fa.gz>
-   Accessed 20/1/22

GENCODE protein-coding translations  

-   gencode.v39.pc\_translations.fa
-   <http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v39.pc_translations.fa.gz>
-   Accessed 20/1/22

GENCODE comprehensive gene annotation (gtf format)  

-   gencode.v39.annotation.gtf
-   <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz>
-   Accessed 20/1/22

ClinVar variant summary table  

-   variant\_summary.txt.gz
-   <https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz>
-   Accessed 25/1/22


## Unpack the AlphaFold2 transcripts

Unzip the AlphaFold2 structures into a folder. In my case, `~/data/af2/`  

    tar -xvf UP000005640_9606_HUMAN.tar


## Get the UniProt ID&rsquo;s corresponding to AF2 structures

Quick shell command to parse the names from the AF2 structure filenames.  

    ls | awk -F"-" '{print $2}' | sort | awk '!a[$0]++' > uniprots.txt

20,504 unique UniProt IDs to use.  


## Get UniProt metadata and matched ENSGs

To match the UniProt accessions to GENCODE, I&rsquo;m using the ENSG&rsquo;s paired with each UniProt structure in the UniProtKB idmapping table. For whatever reasons, matching GENCODE ENST&rsquo;s to UniProtKB&rsquo;s listed ENST&rsquo;s directly is very lossy (even disregarding transcript version numbers). Here we&rsquo;re just getting the idmapping table, filtering to UniProts where we have an AlphaFold2 structure and splitting out the ENSG from `ENSG000001234; ENSG000001235` into separate rows.  

Total unique ENSGs matched to the available AlphaFold structures is 21,649.  

    # get_uniprot_metadata.R
    # M. Silk
    
    library(data.table)
    library(magrittr)
    library(purrr)
    
    
    uniprots = fread("~/data/af2/uniprots.txt", header = F)
    names(uniprots) = "u"
    uniprots %>% nrow
    
    
    idmapping = fread("~/data/HUMAN_9606_idmapping_selected.tab")
    idmapping_names = c("UniProtKB-AC",
                        "UniProtKB-ID",
                        "GeneID (EntrezGene)",
                        "RefSeq",
                        "GI",
                        "PDB",
                        "GO",
                        "UniRef100",
                        "UniRef90",
                        "UniRef50",
                        "UniParc",
                        "PIR",
                        "NCBI-taxon",
                        "MIM",
                        "UniGene",
                        "PubMed",
                        "EMBL",
                        "EMBL-CDS",
                        "Ensembl",
                        "Ensembl_TRS",
                        "Ensembl_PRO",
                        "Additional PubMed")
    names(idmapping) = idmapping_names
    
    # Select only UniProt entries contained in the AlphaFold structure set
    idm = idmapping[`UniProtKB-AC` %in% uniprots$u]
    
    # Identify which are not in the current UniProt metadata file
    uniprots[!u %in% idmapping$`UniProtKB-AC`]
    
    
    # Get table of UniProt ID's and corresponding Ensembl ID's.
    idm %<>% .[, .(u = `UniProtKB-AC`,
                   g = Ensembl)]
    idm %<>% .[!is.na(g)]
    idm %<>% .[, .(u,
                   g = strsplit(g, "; ", fixed = T)[[1]]),
               seq_len(nrow(idm))]
    idm[, seq_len := NULL]
    idm %<>% .[!is.na(g)]
    
    idm$e %>% unique %>% length
    idm$u %>% unique %>% length
    
    fwrite(idm, "uniprot_ensg_matches.txt", sep = "\t")


## Convert GENCODE FASTA to data.table

Convert the GENCODE transcripts into a data.table for easy string matching.  

    # convert_gencode_fasta_to_datatable.R
    # M. Silk
    
    library(data.table)
    library(magrittr)
    library(purrr)
    library(seqinr)
    
    ## gencode_transcripts = read.fasta("~/data/gencode/gencode.v39.pc_transcripts.fa", as.string = T)
    gencode_translations = read.fasta("~/data/gencode/gencode.v39.pc_translations.fa", as.string = T)
    
    
    COLNAMES_TRANSLATIONS = c("ensp", "enst", "ensg",
                       "havana_g", "havana_t", "transcript",
                       "gene", "aalength", "sequence")
    
    
    convert_entry = function(a) {
      header = attr(a, "name") %>% strsplit("|", fixed = T) %>% .[[1]]
      s = a %>% as.character %>% toupper
      vals = c(header, s)
      vals %<>% as.list %>% as.data.table
      names(vals) = COLNAMES_TRANSLATIONS
      return(vals)
    }
    
    
    gencode_translations_dt = lapply(gencode_translations, convert_entry) %>%
      rbindlist(use.names = T, fill = T)
    
    fwrite(gencode_translations_dt, "gencode_translations.txt", sep = "\t")


## Match GENCODE transcripts to UniProt&rsquo;s ENSGs

Next confirm that the GENCODE protein-coding transcripts can be matched to the UniProt-associated ENSGs to see how many GENCODE protein-coding sequences have a valid ENSG.  
ENSGs with a structure: 21,649  
ENSGs with a structure and a GENCODE sequence: 19,434  

GENCODE peptide sequences: 20,365  
GENCODE peptide sequences with an ENSG matched to a structure: 19,434  

Of these, we need to check which of the MANE select transcripts also have a match.  
Loading up the MANE select transcripts, the ENSGs listed are not unique. Some duplicates exist for 56 MANE Plus Clinical transcripts, where an additional transcript to the MANE select has been assigned.  

ENSGs in MANE select: 18,584  
ENSGs in MANE select with an ENSG matched to a structure: 18,462.  

Of the 122 ENSG&rsquo;s with MANE transcripts that don&rsquo;t have a structure in AlphaFold&rsquo;s set, these seem to have a valid UniProt identifier that doesn&rsquo;t match to the GENCODE set. Interesting, and maybe a set of structures to look at after this. For example, `Q8N743` corresponds to the gene KIR3DL3. It&rsquo;s the UniProt identifier used for the AlphaFold structure. But, UniProt lists two ENSG&rsquo;s tied to this; ENSG00000274639 and ENSG00000276196, both of which differ to GENCODE&rsquo;s ENSG of ENSG00000242019. I did a quick look at omitting the gene symbols corresponding to these 122 ENSGs, which only affects 39 of the ClinVar variants in the later steps (negligible).  

    # match_gencode_to_ensgs.R
    # M. Silk
    
    library(data.table)
    library(magrittr)
    library(purrr)
    
    gencode = fread("gencode_translations.txt")
    uniprot = fread("uniprot_ensg_matches.txt")
    
    gencode[, ensg := strsplit(ensg, ".", fixed = T)[[1]][1],
            seq_len(nrow(gencode))]
    
    # which observed ensgs from the structures are in the gencode set?
    uniprot$g %>% unique %>% length
    uniprot[g %in% gencode$ensg]$g %>% unique %>% length
    
    gencode$ensg %>% unique %>% length
    gencode[ensg %in% uniprot$g]$ensg %>% unique %>% length
    
    mane = fread("~/data/gencode/MANE.GRCh38.v0.95.summary.txt")
    mane[, Ensembl_Gene := strsplit(Ensembl_Gene, ".", fixed = T)[[1]][1],
         seq_len(nrow(mane))]
    
    mane$Ensembl_Gene %>% unique %>% length
    mane[duplicated(Ensembl_Gene) | duplicated(Ensembl_Gene, fromLast = T)]
    
    mane$Ensembl_Gene %>% unique %>% length
    mane[Ensembl_Gene %in% uniprot$g]$Ensembl_Gene %>% unique %>% length
    
    mane[!Ensembl_Gene %in% uniprot$g]$Ensembl_Gene %>% unique


## Read AlphaFold2 peptide chains

First get filelist of AF2 structures to read, use this to read in each and collect the peptide sequences into a table.  

    cd ~/data/af2/
    ls *.pdb.gz > filelist.txt

    # read_alphafold_peptide_chains.R
    
    library(data.table)
    library(magrittr)
    library(purrr)
    library(bio3d)
    
    
    # Custom amino acid conversion
    .aa1 = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    
    .aa3 = c("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
             "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
             "PRO", "SER", "THR", "TRP", "TYR", "VAL")
    names(.aa1) = .aa3
    names(.aa3) = as.character(.aa1)
    
    
    aa1to3 = function(a) {
      return(as.character(.aa3[toupper(a)]))
    }
    
    aa3to1 = function(a) {
      return(as.character(.aa1[toupper(a)]))
    }
    
    
    filelist = fread("filelist.txt", header = F)
    names(filelist) = "f"
    
    
    load_pdb = function(p) {
      read.pdb(paste0("~/data/af2/", p))[["atom"]] %>%
        as.data.table %>%
        return
    }
    
    
    extract_chain = function(x) {
      ch = x[, .(resno, resid)] %>% .[!duplicated(.)]
      ch[, resid := aa3to1(resid)]
      return(ch$resid %>% paste(collapse = ""))
    }
    
    
    filelist[, ch := load_pdb(f) %>% extract_chain, 1:nrow(filelist)]
    fwrite(filelist, "alphafold_peptide.txt", sep = "\t")


## Combine AlphaFold2 peptide chains

AlphaFold2 peptide chains that were too large to be calculated in a single run (2700aa < ) were calculated as multiple staggered segments of the overall sequence, overlapping. For example a structure of 3,500aa is calculated as 1 - 2000, 200 - 2200, etc (numbers are just an example).  
Merging the longer sequences into a single entity for each.  

Probably an easier way but I wasn&rsquo;t sure if the fragment sizes would be consistent between each piece generated. So, this uses a recursive function to iteratively align the first two fragments then align to the next repeatedly until a single sequence is returned.  

    # merge_alphafold_fragments.R
    # M. Silk
    
    library(data.table)
    library(magrittr)
    library(purrr)
    library(msa)
    
    alpha = fread("alphafold_peptide.txt")
    alpha[, c("j1", "uniprot", "fragment", "j2") := tstrsplit(f, "-", fixed = T)]
    alpha[, j1 := NULL][, j2 := NULL]
    
    alpha[, fragment := parse_number(fragment)]
    setkey(alpha, uniprot, fragment)
    
    
    # Recursive function to merge each fragment with the following and align
    merge_fragments = function(s) {
      if (length(s) == 1) {
        return(s)
      }
    
      aastrings = AAStringSet(s[1:2])
      names(aastrings) = paste0("F", seq_len(2))
      aastrings.align = msa(aastrings, method = "ClustalW", order = "input", verbose = F) %>%
        as.character %>% strsplit("") %>% as.data.table
      aastrings.align[, M := sort(c(F1, F2))[2], seq_len(nrow(aastrings.align))]
      aastrings.merged = aastrings.align$M %>% paste(collapse = "")
    
      if (length(s) == 2) {
        return(aastrings.merged)
      }
      else {
        merge_fragments(c(aastrings.merged, s[3:length(s)]))
      }
    }
    
    
    alpha %<>% .[, .(f,
                     fragment,
                     fullch = merge_fragments(ch)),
                 uniprot]
    alpha %<>% .[fragment == 1]
    alpha[, fragment := NULL]
    
    fwrite(alpha, "alphafold_peptide_mergedfrags.txt", sep = "\t")


## Match GENCODE to AlphaFold2 chains

Check exact matches between the GENCODE transcripts and AlphaFold peptide sequences. For each GENCODE ENST, match it to an AlphaFold structure based on ENSG (if available), check if strings match, and if so return the UniProt ID corresponding to the exact-matching structure.  

Note that there are 4 ENSGs that each have 2 UniProt structures matchiong based on ENSG which is interesting.  
These are  

1: ENSG00000288649  2  
2: ENSG00000104880  2  
3: ENSG00000287363  2  
4: ENSG00000164112  2  

All MANE select transcripts confirmed to be represented in the GENCODE transcript set.  

Across the GENCODE transcript set, we have 27,181 unique ENSTs with a exactly-matching structure in the AlphaFold database, 1,856 transcripts have ENSGs that aren&rsquo;t in the UniProt idmapping set, and 78,304 GENCODE ENSTs that don&rsquo;t have a match. It is unlikely that we need the full 100,000+ GENCODE transcripts, so need some ways to select out the transcripts with likely clinical relevance (GTEx / pext?)  

Specific to the MANE Select transcripts, 122 do not have a matching ENSG, 16,611 are an exact match and 1,907 are non-matching. Of the non-matching transcripts, 1,340 have an alternate transcript with an exact match, 689 do not.  

Looking specifically at the MANE clinical plus set, 16 have exact matches and 40 do not.  

    # match_gencode_to_alphafold.R
    # M. Silk
    
    library(data.table)
    library(magrittr)
    library(purrr)
    
    uniprot = fread("uniprot_ensg_matches.txt")
    gencode = fread("gencode_translations.txt")
    alpha = fread("alphafold_peptide_mergedfrags.txt")
    mane = fread("~/data/gencode/MANE.GRCh38.v0.95.summary.txt")
    mane[, Ensembl_Gene := strsplit(Ensembl_Gene, ".", fixed = T)[[1]][1],
         seq_len(nrow(mane))]
    mane[, Ensembl_nuc := strsplit(Ensembl_nuc, ".", fixed = T)[[1]][1],
         seq_len(nrow(mane))]
    
    
    # remove transcript version numbers
    gencode[, ensg := strsplit(ensg, ".", fixed = T)[[1]][1],
            seq_len(nrow(gencode))]
    gencode[, enst := strsplit(enst, ".", fixed = T)[[1]][1],
            seq_len(nrow(gencode))]
    
    
    # for a transcript, check if exact matches an alphafold transcript
    # note that there are 4 ENSGs with multiple structures
    
    # for a given GENCODE transcript, check if it matches to any AlphaFold structures
    check_transcript = function(e) {
      print(paste0(":: ", e), quote = F)
      ensgs = gencode[enst == e]$ensg %>% unique
      peptide = gencode[enst == e]$sequence %>% unique
    
    
      if (length(peptide) > 1) {
        return("Error: two different peptide sequences?")
      }
    
      uniprots = uniprot[g %in% ensgs]$u %>% unique
    
      if (length(uniprots) == 0) {
        return("ENSG not in UniProt set")
      }
    
      structures = alpha[uniprot %in% uniprots]$fullch %>% unique
    
      if (length(structures) == 0) {
        return("UniProt ID not in AF2 structures")
      }
    
      if (peptide %in% structures) {
        return("Match")
      }
    
      else {
        return("No match")
      }
    }
    
    
    # ENST have an exact match to an AlphaFold structure?
    gencode[, matchcheck := check_transcript(enst),
            seq_len(nrow(gencode))]
    
    # ENSG have at least 1 matching ENST to AlphaFold structure?
    gencode[, anymatches := any(matchcheck == "Match"),
            ensg]
    
    # Is it in the MANE select transcripts?
    gencode[, mane_select := ifelse(enst %in% mane[MANE_status == "MANE Select"]$Ensembl_nuc, TRUE, FALSE)]
    gencode[, mane_clin := ifelse(enst %in% mane[MANE_status == "MANE Plus Clinical"]$Ensembl_nuc, TRUE, FALSE)]
    gencode[, mane_any := ifelse(enst %in% mane$Ensembl_nuc, TRUE, FALSE)]
    
    # Check counts - all
    gencode[!duplicated(enst)]$matchcheck %>% table
    
    # Check counts - specific to the MANE set (all) and then specifically plus clinical
    gencode[!duplicated(enst)][mane_any == TRUE]$matchcheck %>% table
    gencode[!duplicated(enst)][mane_clin == TRUE]$matchcheck %>% table
    
    # check counts - MANE select with no match but an alternate tx that does
    gencode[!duplicated(enst)][mane_select == TRUE][matchcheck != "Match"]$anymatches %>% table
    
    fwrite(gencode, "gencode_translations_matched.txt", sep = "\t")


## Confirm MANE mismatches are valid

1,907 MANE transcripts don&rsquo;t appear to have an AlphaFold2 structure. It&rsquo;s likely that some of these were missed due to mismatches in ENSG&rsquo;s used between the two datasets, however would be good to confirm that the majority of these have absolutely no matches across any of the AlphaFold2 structures disregarding any UniProt naming (ie. simply matching each MANE select transcript against all AlphaFold2 sequences).  

Only 3 missed MANE select transcripts had a valid AlphaFold2 structure, mismatched by ENSG (or otherwise). So, the remaining 1,904 MANE select transcripts do not have a exactly matching structure.  

    # confirm_mane_mismatches.R
    # M. Silk
    
    library(data.table)
    library(magrittr)
    library(purrr)
    
    alpha = fread("alphafold_peptide_mergedfrags.txt")
    gencode = fread("gencode_translations_matched.txt")
    
    gencode[mane_select == T & matchcheck != "Match"]
    
    # 1990 mismatched transcripts
    
    aa = gencode[mane_select == T & matchcheck != "Match"][sequence %in% alpha$fullch]
    aa
    
    ##                  matchcheck anymatches mane_select mane_clin mane_any
    ##  1: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ##  2: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ##  3: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ##  4: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ##  5: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ##  6: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ##  7: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ##  8:                No match       TRUE        TRUE     FALSE     TRUE
    ##  9: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 10: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 11: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 12: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 13: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 14: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 15: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 16: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 17: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 18: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 19: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 20: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 21:                No match       TRUE        TRUE     FALSE     TRUE
    ## 22:                No match       TRUE        TRUE     FALSE     TRUE
    ## 23: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ## 24: ENSG not in UniProt set      FALSE        TRUE     FALSE     TRUE
    ##                  matchcheck anymatches mane_select mane_clin mane_any
    
    # Only a handful with no ENSG matches, and 3 that do have a valid AlphaFold2 structure.


## Filter and format ClinVar variant summary

Filter the ClinVar variants as per what was suggested to 2\* reviewer status variants, and P/LP only. For P/LP, ClinicalSignificance was parsed for any non-conflicting associations. In this case, I&rsquo;ve included any where there is no conflict, so &ldquo;Pathogenic, other&rdquo; as a category is fine. ClinVar was then filtered to GRCh38 only. Lastly, 9 variants were duplicated. Most of these are PAR X/Y, all 9 were removed for simplicity.  

    # edit_clinvar.R
    # M. Silk
    
    library(data.table)
    library(magrittr)
    library(purrr)
    library(rtracklayer)
    
    clinvar = fread("~/data/clinvar/variant_summary.txt.gz")
    
    
    TWO_STAR_PLUS = c("practice guideline",
                      "reviewed by expert panel",
                      "criteria provided, multiple submitters, no conflicts")
    
    CLINVAR_LP = c("Likely pathogenic",
                   "Likely pathogenic, Affects",
                   "Likely pathogenic, drug response",
                   "Likely pathogenic, other",
                   "Likely pathogenic, risk factor",
                   "Pathogenic",
                   "Pathogenic, Affects",
                   "Pathogenic, association",
                   "Pathogenic, association, protective",
                   "Pathogenic, confers sensitivity",
                   "Pathogenic, drug response",
                   "Pathogenic, drug response, other",
                   "Pathogenic, other",
                   "Pathogenic, protective",
                   "Pathogenic, risk factor",
                   "Pathogenic/Likely pathogenic",
                   "Pathogenic/Likely pathogenic, drug response",
                   "Pathogenic/Likely pathogenic, other",
                   "Pathogenic/Likely pathogenic, risk factor")
    
    clinvar %<>% .[ReviewStatus %in% TWO_STAR_PLUS]
    clinvar %<>% .[ClinicalSignificance %in% CLINVAR_LP]
    clinvar %<>% .[Assembly == "GRCh38"]
    
    dupnames = clinvar[duplicated(Name) | duplicated(Name, fromLast = T)]$Name
    clinvar %<>% .[!Name %in% dupnames]
    
    fwrite(clinvar, "clinvar_twostar_plp.txt", sep = "\t")


## Check ClinVar variants in GENCODE

Since GENCODE provides a GTF which can be easily read into R with GenomicRanges (GRanges), we can format ClinVar to be the same, and match ClinVar against the GENCODE set fairly instantly to ascertain if their coordinates are present within GENCODE, filtered to regions covered by AlphaFold2.  

What we&rsquo;re looking for here is the proportion of ClinVar that isn&rsquo;t represented by the transcripts matching to AlphaFold. So, we&rsquo;ve selected transcripts based on exact matches between their peptide sequences, but for ClinVar matches, we want to be able to also include the variants that are within the exons of transcripts even if these are not part of the resulting CDS sequence. These variants are still relevant, and can also be considered covered by these transcripts. Then for comparison, we want to know which of these variants would otherwise be contained within the exons of the full GENCODE transcript set.  

Within exons:  
&ldquo;af2\_in&rdquo;: within current AF2 structures  
&ldquo;af2\_out&rdquo;: within exons of alternate GENCODE transcripts  

           af2_out
    af2_in  FALSE  TRUE
      FALSE  3321  4641
      TRUE   1881 23298

So, 4,641 are contained within exons of alternate transcripts.  

We also want to know this comparison within the CDS regions, whereby these are the variants that could be visualised on a structure. Reran after changing TARGET\_SEQUENCES at the top of the script.  

Within CDS sequences:  
&ldquo;af2\_in&rdquo;: within current AF2 structures  
&ldquo;af2\_out&rdquo;: within CDS sequencs of alternate GENCODE transcripts  

           af2_out
    af2_in  FALSE  TRUE
      FALSE  4131  3903
      TRUE   3762 21345

Just to summarise this, we have 31,141 unique ClinVar variants (unique by `Name` column) that are P/LP and 2\* or greater reviewer status. 25,107 (75.8%) of these currently overlap with a CDS sequence represented by the AlphaFold transcripts. 3,903 (11.8%) of the ClinVar variants are not represented across any of the coding sequences covered by current AlphaFold transcripts, but are within alternate transcripts within GENCODE. 4,131 (12.5%) of the ClinVar variants are not present within either and likely exist within intronic and non-translated regions. Checking this set of ClinVar, all variants have &rsquo;unusual&rsquo; c. names (all start with g. or n. or have a + or - in the c. name).  

    # check_clinvar_in_gencode.R
    # M. Silk
    
    library(data.table)
    library(magrittr)
    library(purrr)
    library(rtracklayer)
    
    # Edit this to select what sequence types to match GENCODE against ClinVar variants.
    # eg. CDS will select only ClinVar variants within the CDS regions (only coding sequences)
    # exon will select both coding and non-coding exons and will skip introns.
    TARGET_SEQUENCES = "CDS"
    
    gencode = fread("gencode_translations_matched.txt")
    gtf = rtracklayer::import("~/data/gencode/gencode.v39.annotation.gtf")
    clinvar = fread("clinvar_twostar_plp.txt")
    
    # remove transcript versioning from enst
    transcript_id_no_version = gtf$transcript_id %>% lapply(., function(a) (strsplit(a, ".", fixed = T)[[1]][1])) %>% unlist
    mcols(gtf, level = "within")$transcript_id = transcript_id_no_version
    
    # Filter to exons only for coordinates
    gtf %<>% .[mcols(gtf)$type == TARGET_SEQUENCES]
    
    # filter to transcript ID's with a matching AlphaFold structure vs those without
    gtf_in = gtf[mcols(gtf)$transcript_id %in% gencode[matchcheck == "Match"]$enst]
    gtf_out = gtf[!mcols(gtf)$transcript_id %in% gencode[matchcheck == "Match"]$enst]
    
    
    # make GRanges object from clinvar
    clin_gr = makeGRangesFromDataFrame(clinvar,
                                       keep.extra.columns = T,
                                       seqnames.field = "Chromosome",
                                       start.field = "Start",
                                       end.field = "Stop")
    seqlevelsStyle(clin_gr) = "UCSC"
    
    
    # overlap and convert to a column
    oo_in = findOverlaps(clin_gr, gtf_in)
    oo_in_hits = queryHits(oo_in) %>% unique
    oo_in_col = seq_len(nrow(mcols(clin_gr))) %in% oo_in_hits
    
    oo_out = findOverlaps(clin_gr, gtf_out)
    oo_out_hits = queryHits(oo_out) %>% unique
    oo_out_col = seq_len(nrow(mcols(clin_gr))) %in% oo_out_hits
    
    mcols(clin_gr)$af2_in = oo_in_col
    mcols(clin_gr)$af2_out = oo_out_col
    
    # print comparison table
    mcols(clin_gr)[, c("af2_in", "af2_out")] %>% table

