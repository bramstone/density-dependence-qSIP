#------------
#   Stone et al. 
#   2022
#   Nutrients strengthen density dependence of per-capita growth and mortality rates in the soil bacterial community
#   Oecologia
# -----------
# collection and processing of NCBI prokaryote genome data for 16S copy number

# load packages
library(Biostrings)
library(phangorn)
library(castor)
library(data.table)


# load qSIP sequences
qsip <- fread('../data/qsip_growth_mortality.csv')
qsip_seqs <- readDNAStringSet('../data/dim_asv.fasta')
qsip_seqs <- qsip_seqs[qsip[, unique(taxon_id)]]
qsip_tax <- fread('../data/taxonomy_table.csv')[taxonID %in% unique(qsip$taxon_id)]

# download the Ribosomal rRNA database (rrnDB version 5.7 - last modified 2021-01-18)
# URLs:
# https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.7_16S_rRNA.fasta.zip
# https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.7.tsv.zip

rrn_seqs <- readDNAStringSet('rrnDB-5.7_16S_rRNA.fasta')  # 104,233 16S sequences
rrn <- fread('rrnDB-5.7.tsv', fill = T, sep = '\t')  # 20,642 prokaryotes

# download list of NCBI RefSeq prokaryote genomes and limit the rrn to those whose status is complete
# for this publication, download was on 2022-02-25
# URL:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

refseq <- fread('prokaryotes.txt', quote = '')
refseq <- refseq[Status == 'Complete Genome'
                 ][, record_id_num := sub('GCA_|GCF_', '', `Assembly Accession`)
                   ][, record_id_num := sub('\\.\\d+', '', record_id_num)]
# 27,635 complete genomes in refseq as of the download from 2022-02-25


# Organize data-----------------------------------------------------------------

# keep track of Archaea and Bacteria for rooting phylogenetic trees
rrn[grepl('Bacteria|Domain', `RDP taxonomic lineage`) | grepl('Bacteria|Domain', `Data source lineage`), domain := 'Bacteria'
    ][is.na(domain), domain := 'Archaea']

rrn_dat <- rrn[, .(`Data source record id`, `NCBI scientific name`, `16S gene count`, domain)]
setnames(rrn_dat, new = c('record_id', 'ncbi_name', 'num_16s', 'domain'))
rrn_dat[, record_id_num := sub('GCA_|GCF_', '', record_id)
        ][, record_id_num := sub('\\.\\d+', '', record_id_num)]

# match to complete genomes according to RefSeq
rrn_dat <- merge(rrn_dat, refseq[, .(record_id_num, Status)], by = 'record_id_num', all.x = T)

# link GCA / GCF numbers to sequences
seq_names <- fread(text = names(rrn_seqs), sep = '|', fill = T, header = F)
seq_names <- seq_names[, full := names(rrn_seqs)
                       ][, .(V2, V1, full)
                         ][, length_16s := width(rrn_seqs)]
setnames(seq_names, old = c('V2', 'V1', 'full'), new = c('record_id', 'ncbi_name', 'full_seq_name'))
#
seq_names <- merge(seq_names, rrn_dat[, .(record_id, domain, Status)], by = 'record_id', all.x = T)


# quality filtering
seq_names <- seq_names[Status == 'Complete Genome'
                       ][length_16s <= 2000]
rrn_seqs <- rrn_seqs[seq_names$full_seq_name]
rrn_dat <- rrn_dat[record_id %in% unique(seq_names$record_id)]
# 20,048 complete genomes (identified from refseq) in rrnDB


# Align qSIP sequences and rrnDB sequences--------------------------------------

# keep only one 16S sequences per "species" - the longest one as per Louca et al. 2018
rep_seqs <- rrn_seqs[seq_names[domain == 'Bacteria'
                               ][order(-length_16s)
                                 ][!duplicated(ncbi_name), full_seq_name]]
# 16S sequences from 8,518 taxa

# combine with qSIP taxa 16S sequences (V4 region)
seqs <- c(qsip_seqs, rep_seqs)

# align sequences using predefined guide tree to prevent AlignSeqs making one and running out of memory
gT <- lapply(order(width(seqs)),
             function(x) {
               attr(x, "height") <- 0
               attr(x, "label") <- names(seqs)[x]
               attr(x, "members") <- 1L
               attr(x, "leaf") <- T
               x
             })
#
attr(gT, "height") <- 0.5
attr(gT, "members") <- length(seqs)
class(gT) <- "dendrogram"

set.seed(22822)
align_seqs <- DECIPHER::AlignSeqs(seqs,
                                  guideTree = gT,  # recommended for > 10,000 sequences - limits memory use
                                  iterations = 0,   # can be 0 when using a guide tree
                                  refinements = 0,  # can be 0 when using a guide tree
                                  restrict = c(-500, 2, 10),  # limits run-time by avoiding unnecessary comparisons
                                  anchor = 0.7,  # higher values limits the memory use, default: 0.7
                                  processors = 2)

# convert the number of mismatched base bairs into "evolutionary distance" using the rules of the Jukes Cantor model
# the output value represents an estimate of the proportion of nucleotide substitutions in the sequence
set.seed(3122)
seq_dist <- dist.ml(phyDat(as(align_seqs, 'matrix'), type = 'DNA'), model = 'JC69', gamma = T)


#  Phylogenetic distance--------------------------------------------------------

# Get distances directly from matrix.............
# this works because the tree takes two taxa with distance of e.g. 0.2 and
# creates a shared node between them with 2 edges each of length 0.1.
# This should actually get to the same result because Louca et al. used FastTree which
# estimates distances between aligned sequences also using the Jukes-Cantor model, as was applied above

# may take a minute or two - there are 58 million distances...
nstd <- data.table(t(combn(names(seqs), 2)), 
                   dist = c(seq_dist))

# limit to distances between qSIP taxa and reference genomes, then get minimum distance
nstd <- nstd[!V2 %in% names(qsip_seqs)
             ][!V1 %in% names(rep_seqs)
               ][order(dist)
                 ][, .(best_match = V2[1], nstd = dist[1]), by = V1]
setnames(nstd, 'V1', 'taxonID')

# match with SILVA taxonomy
nstd <- merge(nstd, qsip_tax, all.x = T, by = 'taxonID')


# combine with 16S gene copy number
copy_num <- merge(rrn_dat[, .(record_id, num_16s)], 
                  seq_names[full_seq_name %in% names(rep_seqs), .(record_id, full_seq_name)], 
                  all.y = T, by = 'record_id')

nstd <- merge(nstd, copy_num, all.x = T, by.x = 'best_match', by.y = 'full_seq_name')
setcolorder(nstd, c('taxonID', 'record_id', 'best_match', 'num_16s', 'nstd'))

# # save output
# fwrite(nstd, file = '../data/qsip_taxa_16s_gcn_nstd.csv', quote = F)


# plot how well we can trust our estimates of 16S rRNA gene copy number assignment
ggplot(nstd[taxonID %in% unique(qsip$taxon_id)], aes(num_16s, nstd)) +
  geom_jitter(width = 0.3, shape = 21, fill = gray(.85)) +
  geom_hline(yintercept = 0.15, color = hsv(0, .8, .6)) +
  ylab('NSTD (Prop. substitutions per site)') +
  xlab('16S rRNA GCN per genome estimate') +
  labs(tag = 'A') +
  annotate('text', x = Inf, y = Inf, 
           size = 3.25,
           label = paste('Perc. ASVs above 0.15 NSTD: ', 
                         round(nstd[taxonID %in% unique(qsip$taxon_id), sum(nstd > 0.15) / .N], 3) * 100, '%'), 
           hjust = 1.1, vjust = 2)
# 37.4% (851) of ASVs have untrustworthy 16S gene copy numbers

# what proportion of the community do imprecise copy number estimates make up?
qsip <- merge(qsip, nstd[, .(taxonID, num_16s, nstd)], all.x = T, by.x = 'taxon_id', by.y = 'taxonID')

qsip[, 1 - (sum(abund_16s[nstd > 0.15]) / sum(abund_16s)), by = .(ecosystem, treatment, rep, sample_id)]
