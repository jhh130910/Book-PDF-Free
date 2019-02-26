
# --*-- coding: utf-8 --*--

import argparse
import HTSeq
import random
from Bio import SeqIO
import os.path
import os
import rpy2.robjects as robjects

#guoyang 20130723

parser = argparse.ArgumentParser(description='coverage density distribution -- circos guoyang@novogene.cn')
parser.add_argument('--bam',required=True,help="the bam file from mapping")
parser.add_argument('--fa',required=True,help="genome fasta file")
parser.add_argument('--n',type=int,default=10,help="n chrom or contig will be drew default:10")
parser.add_argument('--iv',type=int,default=1000,help="the interval for coverage default:1000")
parser.add_argument('--r',type=int,default=10000,help="r reads for reads distribution track default:10000")
parser.add_argument('--output',default=None,help="Change the output filename.")
args=parser.parse_args()

bam=args.bam
assert os.path.isfile(bam)

fa=args.fa
assert os.path.isfile(fa)

n =args.n
iv=args.iv
r =args.r
output=args.output

mean=robjects.r['mean']
sd=robjects.r['sd']

#genome file length
len_genome={}
chrom_count=0
for seq_record in SeqIO.parse(fa,'fasta'):
        chrom_count+=1
        if len(seq_record.seq) not in len_genome:
                len_genome[len(seq_record.seq)]=[]
        len_genome[len(seq_record.seq)].append(seq_record.id)
n=min(n,chrom_count)

b2s_len=reversed(sorted(len_genome.keys()))
genome_len={}
flag=0
for i in b2s_len:
        for j in len_genome[i]:
                genome_len[j]=i
                if len(genome_len) == n:
                        flag=1
                        break
        if flag == 1:
                break

total_len=sum(genome_len.values())
if total_len > 3116677:
        units=int(round(float(total_len)/3116677))*10000
else:
        units=10000

genome=[]
flag=0
for each in genome_len:
        if flag == 0:
                genome.append('chr - '+each+' '+each+' 0 '+str(genome_len[each]-1)+' vlgrey\n')
                flag=1
                continue
        if flag == 1:
                genome.append('chr - '+each+' '+each+' 0 '+str(genome_len[each]-1)+' grey\n')
                flag=0
open('genome.txt','w').writelines(genome)

#bam file
bam_file=HTSeq.BAM_Reader(bam)
coverage = HTSeq.GenomicArray( "auto", stranded=True, typecode="i" )
reads=[]
read_count=0
for almnt in bam_file:
        if almnt.aligned and almnt.iv.chrom in genome_len:
                coverage[almnt.iv]+=1
                read_count+=1
                if almnt.iv.strand == '+':
                        reads.append(almnt.iv.chrom+'\t'+str(almnt.iv.start)+'\t'+str(almnt.iv.end-1)+'\tstroke_color=vdred\n')
                if almnt.iv.strand == '-':
                        reads.append(almnt.iv.chrom+'\t'+str(almnt.iv.start)+'\t'+str(almnt.iv.end-1)+'\tstroke_color=vdblue\n')

#read tile distribution
r=min(r,read_count)
read_sample=random.sample(reads,r)
open('reads.txt','w').writelines(read_sample)

#line density
line_plus=[]
line_minus=[]
coverage_plus_total=[]
coverage_minus_total=[]
for each in genome_len:
        flag=1
        start=-iv
        while True:
                start+=iv
                if start+iv <= genome_len[each]:
                        window_plus=HTSeq.GenomicInterval(each,start,start+iv,'+')
                        coverage_plus=list(coverage[window_plus])
                        coverage_plus_mean=mean(robjects.IntVector(coverage_plus))[0]
                        coverage_plus_total.append(coverage_plus_mean)
                        line_plus.append(each+'\t'+str(start)+'\t'+str(start+iv-1)+'\t'+str(coverage_plus_mean)+'\n')
                        window_minus=HTSeq.GenomicInterval(each,start,start+iv,'-')
                        coverage_minus=list(coverage[window_minus])
                        coverage_minus_mean=mean(robjects.IntVector(coverage_minus))[0]
                        coverage_minus_total.append(coverage_minus_mean)
                        line_minus.append(each+'\t'+str(start)+'\t'+str(start+iv-1)+'\t'+str(coverage_minus_mean)+'\n')
                else:
                        window_plus=HTSeq.GenomicInterval(each,start,genome_len[each],'+')
                        coverage_plus=list(coverage[window_plus])
                        coverage_plus_mean=mean(robjects.IntVector(coverage_plus))[0]
                        coverage_plus_total.append(coverage_plus_mean)
                        line_plus.append(each+'\t'+str(start)+'\t'+str(genome_len[each]-1)+'\t'+str(coverage_plus_mean)+'\n')
                        window_minus=HTSeq.GenomicInterval(each,start,genome_len[each],'-')
                        coverage_minus=list(coverage[window_minus])
                        coverage_minus_mean=mean(robjects.IntVector(coverage_minus))[0]
                        coverage_minus_total.append(coverage_minus_mean)
                        line_minus.append(each+'\t'+str(start)+'\t'+str(genome_len[each]-1)+'\t'+str(coverage_minus_mean)+'\n')
                        flag=0
                if flag == 0:
                        break

#filter
plus_cutoff=mean(robjects.FloatVector(coverage_plus_total))[0]+3*sd(robjects.FloatVector(coverage_plus_total))[0]
minus_cutoff=mean(robjects.FloatVector(coverage_minus_total))[0]+3*sd(robjects.FloatVector(coverage_minus_total))[0]
line_plus_filtered=[]
line_minus_filtered=[]
for eachLine in line_plus:
        if float(eachLine.split()[-1].strip()) < plus_cutoff:
                line_plus_filtered.append(eachLine)
for eachLine in line_minus:
        if float(eachLine.split()[-1].strip()) < minus_cutoff:
                line_minus_filtered.append(eachLine)

open('line_plus.txt','w').writelines(line_plus_filtered)
open('line_minus.txt','w').writelines(line_minus_filtered)

ideogram_position='''
radius           = 0.775r
thickness        = 30p
fill             = yes
fill_color       = black
stroke_thickness = 2
stroke_color     = black
'''
open('ideogram.position.conf','w').write(ideogram_position)

ticks='''

show_ticks          = yes
show_tick_labels    = yes

<ticks>

radius           = dims(ideogram,radius_outer)
orientation      = out
label_multiplier = 1e-4
color            = black

thickness        = 2p
label_offset     = 5p
format           = %d

<tick>
spacing        = 1u
show_label     = no
size           = 10p
</tick>

<tick>
spacing        = 5u
show_label     = yes
label_size     = 15p
size           = 15p
</tick>

<tick>
spacing        = 10u
show_label     = yes
label_size     = 20p
size           = 20p
</tick>

</ticks>

'''
open('ticks.conf','w').write(ticks)

ideogram='''

<ideogram>

<spacing>
default = 0.01r
#break   = 0.5r
</spacing>

<<include ideogram.position.conf>>
<<include ideogram.label.conf>>

radius*       = 0.75r

</ideogram>

'''
open('ideogram.conf','w').write(ideogram)

ideogram_label='''
show_label       = yes
label_font       = default
label_radius     = dims(ideogram,radius) + 0.1r
label_size       = 24
label_parallel   = yes
label_case       = lower

'''
open('ideogram.label.conf','w').write(ideogram_label)

circos='''

<<include etc/colors_fonts_patterns.conf>>

<<include ideogram.conf>>
<<include ticks.conf>>

karyotype   = genome.txt

<image>
<<include etc/image.conf>>
</image>

chromosomes_units = %s

chromosomes_display_default = yes

<plots>

<plot>
type      = tile
file      = reads.txt
r1        = 0.95r
r0        = 0.75r
layers    = 4
margin    = 0.2u
thickness = 4
padding   = 2
layers_overflow  = grow
orientation      = in
stroke_thickness = 1p
#stroke_color     = grey
#color            = orange

<backgrounds>
<background>
color     = vvlgrey
#y0        = 0.006
</background>

</backgrounds>

</plot>

<plot>
type        = histogram
file        = line_plus.txt
r1          = 0.70r
r0          = 0.55r
#min         = 0
#max         = 1
extend_bin  = yes
fill_color  = vdorange
color       = vdorange
thickness   = 0
orientation = out

<axes>
<axis>
color     = lgrey_a2
thickness = 1
spacing   = 0.1r
</axis>
</axes>

</plot>

<plot>
type        = histogram
file        = line_minus.txt
r1          = 0.55r
r0          = 0.40r
#max         = 1
#min         = 0
extend_bin  = yes
fill_color  = vdgreen
color       = vdgreen
thickness   = 0
orientation = in

<axes>
<axis>
color     = lgrey_a2
thickness = 1
spacing   = 0.1r
</axis>
</axes>
</plot>

</plots>

<<include etc/housekeeping.conf>>

''' % (units)
open('circos.conf','w').write(circos)

if output:
        os.system('circos -conf circos.conf -outputfile %s' % output)
else:
        os.system('circos -conf circos.conf')
