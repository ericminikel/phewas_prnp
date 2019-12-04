options(stringsAsFactors=F)
if(interactive()) {
  setwd('~/d/sci/src/phewas_prnp')
}
library(sqldf)

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}

path = 'data/otg/'

assoc = read.table(paste0(path,'ENSG00000171867-associated-studies.tsv'),sep='\t',header=T,quote='',comment.char='')
colnames(assoc) = gsub('[^a-z0-9_*]','_',tolower(colnames(assoc)))
coloc = suppressWarnings(read.table(paste0(path,'ENSG00000171867-colocalising-studies.tsv'),sep='\t',header=T,quote='',comment.char=''))
colnames(coloc) = gsub('[^a-z0-9_*]','_',tolower(colnames(coloc)))


# coloc data by snp+study for the 2 studies with coloc reported above
qtlcoloc_files = list.files(path = path, pattern='qtl-coloc-*')
qtlcoloc = read.table(paste0(path,qtlcoloc_files[1]),sep='\t',header=T,quote='',comment.char='')[0,]
qtlcoloc$variant = character(0)
for (qfile in qtlcoloc_files) {
  temp = read.table(paste0(path,qfile),sep='\t',header=T,quote='',comment.char='')
  temp$variant = gsub('.+-','',gsub('\\.tsv','',qfile))
  qtlcoloc = rbind(qtlcoloc, temp)
}
colnames(qtlcoloc) = gsub('[^a-z0-9_*]','_',tolower(colnames(qtlcoloc)))

# how many lead SNPs are there in the assoc table
assoc_snps = sqldf("
                   select   index_variant_id, count(*) n
                   from     assoc
                   group by 1
                   order by 1
                   ;")
assoc_snps
# cat(paste0('https://genetics.opentargets.org/variant/',assoc_snps$index_variant_id,'\n'))
# used the above code to manually grab the TSV files of v2g data for all 32 unique lead SNPs

# v2g data for the 32 lead SNPs for the 59 GWAS hits listed for PRNP
v2g_files = list.files(path = path, pattern='*-assigned-genes-summary.tsv')
v2g = read.table(paste0(path,v2g_files[1]),sep='\t',header=T,quote='',comment.char='')[0,]
v2g$variant = character(0)
for (qfile in v2g_files) {
  temp = read.table(paste0(path,qfile),sep='\t',header=T,quote='',comment.char='')
  temp$variant = gsub('-assigned-genes-summary\\.tsv','',qfile)
  v2g = rbind(v2g, temp)
}
colnames(v2g) = gsub('[^a-z0-9_*]','_',tolower(colnames(v2g)))

# some descriptive stats to understand the datasets
sqldf("
select   variant, sum(overall_v2g) overall_v2g_total
from     v2g
group by 1
order by 1
;")
# they don't add to 1, although they are distributed reasonably tightly around 1
# i wonder if it is a sum of not-quite-mutually-exclusive probabilities

# how strong is PRNP?
sqldf("
select   gene, variant, overall_v2g
from     v2g
where    gene = 'PRNP'
order by 1, 2
;")
# .26 and .44 for the two SNPs associated to CJD; all < 0.1 for the SNPs associated to other phenotypes

# how strong is PRNP as a proportion of total?
sqldf("
select  variant, sum(case when gene = 'PRNP' then overall_v2g else 0 end) / sum(overall_v2g) prnp_proportion
from     v2g
group by 1
order by 1
;")
# 0.28 and 0.40 for the CJD-associated SNPs, all < 0.1 for SNPs associated to other phenotypes
# and 0.00 (i.e. not even a row for PRNP) for many of them

# how about in qtlcoloc for all sources of info
sqldf("
select   variant, sum(case when gene = 'PRNP' then h4/h3 else 0 end) / sum(h4/h3) prnp_proportion
from     qtlcoloc
group by 1
order by 1
;")

# and just for GTEx?
sqldf("
select   variant, sum(case when gene = 'PRNP' then h4/h3 else 0 end) / sum(h4/h3) prnp_proportion
from     qtlcoloc
where    source = 'GTEX_v7'
group by 1
order by 1
;")

# ------------

# summary table
assoc$pos = as.integer(gsub('20_','',gsub('_[ACGT]+_[ACGT]+','',assoc$index_variant_id)))
prnp_tss = 4686350
assoc$prnp_tss_dist = assoc$pos - prnp_tss

prnp_v2g = sqldf("
                 select   variant, sum(case when gene = 'PRNP' then overall_v2g else 0 end) / sum(overall_v2g) prnp_proportion
                 from     v2g
                 group by 1
                 order by 1
                 ;")

prnp_summary = sqldf("
                   select   a.index_variant_id, a.prnp_tss_dist, count(*) n, group_concat(reported_trait,', ') phenotypes, v.prnp_proportion, a.pos
                   from     assoc a, prnp_v2g v
                   where    a.index_variant_id = v.variant
                   group by 1
                   order by 1
                   ;")

prnp_summary$phenotypes = paste0(substr(prnp_summary$phenotypes,1,60),ifelse(nchar(prnp_summary$phenotypes)>60,'...',''))

write.table(prnp_summary,'output/prnp_otg_summary.tsv', sep='\t', col.names=T, row.names=F, quote=F)


gstraint = read.table('data_large/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz',sep='\t',header=T) # read.table natively handles bgzipped files : )

png('figures/prnp_otg_v2g_view.png',width=800,height=500,res=100)

xlims = range(prnp_summary$pos) + c(-1.25e5,+1e5)

grch38_37_offset = 4699605 - 4680251 # based on rs1799990
gns = subset(gstraint, chromosome==20 & end_position > min(xlims) & end_position < max(xlims))
gns$midpoint = (gns$start_position + gns$end_position) / 2
gns$start_position = gns$start_position + grch38_37_offset
gns$end_position = gns$end_position + grch38_37_offset
gns$midpoint = gns$midpoint + grch38_37_offset

gns = gns[with(gns, order(start_position)),]
gns$level = as.integer(NA)
gns$level[1] = 0
for (i in 2:nrow(gns)) {
  if (gns$start_position[i] < gns$end_position[i-1] + 2e5) {
    if (gns$start_position[i] > max(gns$end_position[gns$level==0] + 2e5, na.rm=T)) {
      gns$level[i] = 0
    } else {
      gns$level[i] = gns$level[i-1] + 1
    }
  } else {
    gns$level[i] = 0
  }
}

prnp_summary$color = '#777777'
prnp_summary$color[grepl('Creutz',prnp_summary$phenotypes)] = '#FF2019'


layout_matrix = matrix(c(1,2),nrow=2,byrow=T)
layout(layout_matrix, heights=c(1,1))

par(mar=c(0,4,2,1))
plot(NA, NA, xlim=xlims, ylim=c(0,0.5), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=(0:5)/10, labels=percent((0:5)/10), lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=3.0, text='PRNP share of overall v2g')
points(x=prnp_summary$pos, y=prnp_summary$prnp_proportion, type='h', lwd=3, lend=1, col=prnp_summary$color)
legend('topright', c('prion disease','all other phenotypes'), col=c('#FF2019','#777777'), text.col=c('#FF2019','#777777'), lwd=3)

par(mar=c(4,4,0,1))
plot(NA, NA, xlim=xlims, ylim=c(-0.1,1.1), axes=F, ann=F, xaxs='i', yaxs='i')
segments(x0=gns$start_position, x1=gns$end_position, y0=1-gns$level/max(gns$level), y1=1-gns$level/max(gns$level), lwd=1)
text(x=gns$midpoint, y=1-gns$level/max(gns$level)+0.05, pos=1, labels=gns$gene, cex=0.5)
axis(side=1, at=seq(3e6,8e6,by=5e5), labels=formatC(seq(3e6,8e6,by=5e5), format='fg', big.mark=','))
mtext(side=1, line=2.5, text='chr20 position (GRCh38)')

dev.off()




# plot to illustrate colocalization concept
ccf_col = '#00CDCD'
eso_col = '#55220D'
ret_col = '#660000'

# GTEx v8 PRNP esophagus mucosa eQTL
gtex = read.table('data_large/gtex/Esophagus_Mucosa.v8.signif_variant_gene_pairs.txt.gz',sep='\t',header=T)
gtex = gtex[gtex$gene_id %in% c('ENSG00000171867.16','ENSG00000088826.17'),]
gtex$gene[gtex$gene_id == 'ENSG00000171867.16'] = 'PRNP'
gtex$gene[gtex$gene_id == 'ENSG00000088826.17'] = 'SMOX'
gtex$chr = gsub('_.*','',gsub('chr','',gtex$variant_id))
gtex$pos = as.integer(gsub('chr.+_','',gsub('_[ACGT]+_[ACGT]+_b38','',gtex$variant_id)))


apparent_gtex_threshold = max(gtex$pval_nominal)
# Astle 2016 summary stats: https://www.phpc.cam.ac.uk/ceu/haematological-traits/
# specifically reticulocyte count: ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/human/2017-12-12/hematological_traits/ret/ret_N170641_narrow_form.tsv.gz
# which was linked from ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/human/2017-12-12/hematological_traits/ret/
# then: 
# gzcat data/astle2016/ret_N170641_narrow_form.tsv.gz | head -1 > data/astle2016/20.tsv
# gzcat data/astle2016/ret_N170641_narrow_form.tsv.gz | grep ^20 | head -100000 >> data/astle2016/20.tsv
ret = read.table('data/astle2016/ret_20.tsv',sep='\t',header=T)
colnames(ret) = gsub('[^a-z0-9_]','_',tolower(colnames(ret)))


png('figures/prnp_colocalization_demo.png',width=800,height=500,res=100)

layout_matrix = matrix(c(1,2),nrow=2,byrow=T)
layout(layout_matrix, heights=c(1,1))

par(mar=c(0,4,2,8))
plot(NA, NA, xlim=xlims, ylim=c(0,30), axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=0:6*5, las=2)
mtext(side=2, line=2.5, text='-log10(P)')
abline(h=-log10(5e-8), lwd=1, lty=3)
mtext(side=4, at=-log10(5e-8), las=2, line=0, text='GWsig threshold')
abline(h=-log10(apparent_gtex_threshold), lwd=1, lty=3)
mtext(side=4, at=-log10(apparent_gtex_threshold), las=2, line=0, text='GTEx FDR threshold')
points(ret$bp, -log10(ret$p), pch=20, col=ret_col)
points(gtex$pos, -log10(gtex$pval_nominal), pch=20, col=ccf_col)
#abline(v=4182655)
mtext(side=3,at=gtex$pos[gtex$pval_nominal==min(gtex$pval_nominal[gtex$gene=='SMOX'])]-150000, line=0, text=expression(italic('SMOX')~'eQTL  '), col=ccf_col, cex=0.9)
mtext(side=3,at=gtex$pos[gtex$pval_nominal==min(gtex$pval_nominal[gtex$gene=='PRNP'])]+150000, line=0, text=expression(italic('  PRNP')~'eQTL'), col=ccf_col, cex=0.9)
segments(x0=gtex$pos[gtex$pval_nominal==min(gtex$pval_nominal[gtex$gene=='SMOX'])]-150000, x1=gtex$pos[gtex$pval_nominal==min(gtex$pval_nominal[gtex$gene=='SMOX'])]-50000, y0=28, y1=15, col=ccf_col)
segments(x0=gtex$pos[gtex$pval_nominal==min(gtex$pval_nominal[gtex$gene=='PRNP'])]+150000, x1=gtex$pos[gtex$pval_nominal==min(gtex$pval_nominal[gtex$gene=='PRNP'])]+50000, y0=28, y1=9, col=ccf_col)

legend('topright',c('esophagus mucosa eQTLs','reticulocyte count'),col=c(ccf_col,ret_col),text.col=c(ccf_col,ret_col),pch=20,bty='n')

par(mar=c(4,4,0,8))
plot(NA, NA, xlim=xlims, ylim=c(-0.1,1.1), axes=F, ann=F, xaxs='i', yaxs='i')
segments(x0=gns$start_position, x1=gns$end_position, y0=1-gns$level/max(gns$level), y1=1-gns$level/max(gns$level), lwd=1)
text(x=gns$midpoint, y=1-gns$level/max(gns$level)+0.05, pos=1, labels=gns$gene, cex=0.5)
axis(side=1, at=seq(3e6,8e6,by=5e5), labels=formatC(seq(3e6,8e6,by=5e5), format='fg', big.mark=','))
mtext(side=1, line=2.5, text='chr20 position (GRCh38)')

dev.off()


