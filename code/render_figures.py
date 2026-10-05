"""Render three main figures from frozen display inputs; no model fitting."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import xml.etree.ElementTree as ET
import cairosvg

# 1. 输入与统一样式：所有数据来自冻结图源，只有版式变化。
NS = 'http://www.w3.org/2000/svg'
ET.register_namespace('', NS)
W = 510
BLUE, TEAL, GRAY, ORANGE = '#286782', '#187e89', '#edf0f2', '#AD633C'
INK, MUTED = '#222222', '#485C64'

def node(parent, tag, **attrs):
    return ET.SubElement(parent, '{'+NS+'}'+tag, {k.replace('_','-'):str(v) for k,v in attrs.items()})

def text(parent, x, y, value, size=7, color=INK, weight='normal', anchor='start', **attrs):
    e = node(parent,'text',x=x,y=y,font_family='Arial',font_size=size,fill=color,
             font_weight=weight,text_anchor=anchor,**attrs)
    e.text = str(value)
    return e

def canvas(height):
    root = ET.Element('{'+NS+'}svg',width=f'{W}pt',height=f'{height}pt',viewBox=f'0 0 {W} {height}')
    node(root,'rect',x=0,y=0,width=W,height=height,fill='white')
    return root

def heading(root, letter, y, title):
    text(root,7,y,letter,10,weight='bold')
    text(root,26,y,title,8.5,weight='bold')

def box(root,x,y,width,height,lines,fill='#EDF4F7',bold_lines=()):
    node(root,'rect',x=x,y=y,width=width,height=height,rx=4,fill=fill,stroke='#70929D',stroke_width=.65)
    first = y+height/2 - (len(lines)-1)*6 + 2.5
    for i,line in enumerate(lines):
        text(root,x+width/2,first+12*i,line,7.5 if i in bold_lines else 7,
             weight='bold' if i in bold_lines else 'normal',anchor='middle')

def arrow(root,x1,y1,x2,y2):
    node(root,'path',d=f'M{x1} {y1} L{x2} {y2}',stroke=BLUE,stroke_width=1,fill='none')
    node(root,'path',d=f'M{x2-4} {y2-2.5} L{x2} {y2} L{x2-4} {y2+2.5}',stroke=BLUE,stroke_width=1,fill='none')

def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def save(root,number,out):
    svg=ET.tostring(root,encoding='utf-8',xml_declaration=True)
    (out/f'Figure_{number}.svg').write_bytes(svg)
    cairosvg.svg2pdf(bytestring=svg,write_to=str(out/f'Figure_{number}.pdf'))
    cairosvg.svg2png(bytestring=svg,write_to=str(out/f'Figure_{number}.png'),scale=600/96)

def render(source,out):
    out.mkdir(parents=True,exist_ok=True)
    inputs=['Fig1_analysis_populations.tsv','Fig2_M2_membership_sensitivity.tsv',
            'Fig3_a_pair_composition.tsv','Fig3_b_chr19_incidence_matrix.tsv',
            'Fig3_c_chr6_variant_trait_counts.tsv']
    hashes={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in inputs}
    populations={r['population']:int(r['n']) for r in read(source/inputs[0])}
    assert populations=={'original_model':1652,'masked_model':1259,'selected_input_audit':283}
    pair_rows=read(source/inputs[2]);counts=[int(r['n_pairs']) for r in pair_rows]
    assert counts==[464,80,2680] and {int(r['denominator']) for r in pair_rows}=={3224}
    events=read(source/'Fig3_b_chr19_event_edges.tsv')
    hashes['Fig3_b_chr19_event_edges.tsv']=hashlib.sha256((source/'Fig3_b_chr19_event_edges.tsv').read_bytes()).hexdigest()
    assert len(events)==154

    # 2. 图 1：两个独立总体；两种触发条件按 OR 组合，不暗示因果分解。
    root=canvas(346)
    heading(root,'a',18,'Source-mapped gene support')
    box(root,16,34,140,44,['11 non-HCC liver traits','Fixed source cis records'])
    box(root,185,34,140,44,['Aggregated gene evidence','Top 5% defines support'])
    box(root,354,34,140,44,['56 seed gene identities','Excluded from comparisons'])
    arrow(root,157,56,183,56);arrow(root,326,56,352,56)
    text(root,255,94,'Mapped input ≠ significant eQTL link or effector-gene assignment',color=MUTED,anchor='middle')
    heading(root,'b',119,'Distinct analysis populations')
    box(root,16,132,228,69,['Conditional breadth comparison',
        f"{populations['original_model']} genes → {populations['masked_model']} after masking",
        '≥4 versus 1–3 supporting traits'],bold_lines=(0,1))
    box(root,266,132,228,69,['Selected input audit',
        '283 genes · 1465 support events','3224 within-gene trait pairs'],bold_lines=(0,1))
    text(root,255,216,'HCC association assessed separately in European and East Asian data',color=MUTED,anchor='middle')
    heading(root,'c',243,'Seed-gene-body membership mask (Rule A)')
    box(root,16,259,146,40,['Annotated TSS overlaps','a seed gene body'],fill='#F2F5EF')
    text(root,174,282,'OR',7,weight='bold',anchor='middle')
    box(root,186,259,153,40,['Any frozen mapped input','overlaps a seed gene body'],fill='#F2F5EF')
    # The OR separator joins the conditions; the arrow targets the whole-gene decision.
    arrow(root,341,279,364,279)
    box(root,367,259,127,40,['Exclude the entire gene'],fill='#E6EFF3',bold_lines=(0,))
    text(root,255,320,'Support flags, scores and global HCC ranks remain fixed.',color=MUTED,anchor='middle')
    text(root,255,334,'Population-composition sensitivity; not a causal decomposition.',color=MUTED,anchor='middle')
    save(root,1,out)

    # 3. 图 2：保留四组点估计、区间、祖源色与原始/屏蔽形状。
    models={(r['ancestry'],r['part']):r for r in read(source/inputs[1])}
    order=[('EUR','original'),('EUR','retained'),('EAS','original'),('EAS','retained')]
    assert len(models)==4
    root=canvas(262)
    text(root,16,18,'M2 conditional breadth before and after membership masking',8.5,weight='bold')
    text(root,351,39,'Estimate [95% CI]',7,weight='bold')
    xp=lambda value:101+(value+.9)/2.3*225
    node(root,'path',d=f'M{xp(0)} 46 V221',stroke='#777777',stroke_width=.8,stroke_dasharray='3 3')
    model_records=[]
    for y,(ancestry,part) in zip([65,105,163,203],order):
        row=models[(ancestry,part)];estimate,low,high=[float(row[k]) for k in ['estimate','ci_low','ci_high']]
        assert row['analysis_role']=='current_manuscript' and low<=estimate<=high
        color=BLUE if ancestry=='EUR' else ORANGE
        label='original' if part=='original' else 'masked'
        text(root,87,y+2.5,f'{ancestry} · {label}',anchor='end')
        node(root,'path',d=f'M{xp(low)} {y} H{xp(high)} M{xp(low)} {y-3} V{y+3} M{xp(high)} {y-3} V{y+3}',stroke=color,stroke_width=1.3,fill='none',data_role='confidence_interval',data_ancestry=ancestry,data_population=part,data_low=low,data_high=high)
        if part=='original':
            node(root,'circle',cx=xp(estimate),cy=y,r=2.6,fill=color,data_role='estimate',data_estimate=estimate)
        else:
            node(root,'rect',x=xp(estimate)-2.6,y=y-2.6,width=5.2,height=5.2,fill=color,data_role='estimate',data_estimate=estimate)
        label=f'{estimate:.3f} [{low:.3f}, {high:.3f}]'
        text(root,351,y+2.5,label)
        model_records.append(dict(ancestry=ancestry,population=part,estimate=estimate,ci_low=low,ci_high=high,label=label))
    node(root,'path',d=f'M{xp(-.8)} 223 H{xp(1.4)}',stroke=INK,stroke_width=.8,fill='none')
    for v in [-.5,0,.5,1]:
        node(root,'path',d=f'M{xp(v)} 223 V227',stroke=INK,stroke_width=.8)
        text(root,xp(v),237,f'{v:.1f}'.replace('-','−'),anchor='middle')
    text(root,218,249,'Conditional breadth contrast (rank-evidence units)',anchor='middle')
    save(root,2,out)

    # 4. 图 3：压缩说明区，不删减基因/变异；颜色表示事件有无。
    root=canvas(620)
    heading(root,'a',16,'Minimum-P input composition')
    text(root,26,30,'283 selected candidates | 1465 support events | 3224 within-gene pairs')
    defs=node(root,'defs')
    x=26;bar_width=468
    for label,count,fill in zip(['shared_without_floor','shared_floor','disjoint'],counts,[BLUE,BLUE,'#D6DADD']):
        width=bar_width*count/3224
        node(root,'rect',x=x,y=43,width=width,height=15,fill=fill,stroke='white',stroke_width=.6,data_category=label,data_count=count)
        if label=='shared_floor':
            clip=node(defs,'clipPath',id='floor_clip')
            node(clip,'rect',x=x,y=43,width=width,height=15)
            hatch=node(root,'g',clip_path='url(#floor_clip)')
            for offset in range(-15,int(width)+16,4):
                node(hatch,'path',d=f'M{x+offset} 58 L{x+offset+15} 43',stroke='white',stroke_width=.7)
        x+=width
    node(root,'path',d='M26 64 H494',stroke=INK,stroke_width=.65)
    for v in range(0,101,20):
        x=26+bar_width*v/100
        node(root,'path',d=f'M{x} 64 V67',stroke=INK,stroke_width=.65)
        text(root,x,76,str(v),anchor='middle')
    text(root,260,88,'Within-gene cross-trait pairs (%)',anchor='middle')
    text(root,26,101,'Shared minimum input: 544 / 3224 (16.9%)',color=BLUE,weight='bold')
    text(root,258,101,'464 without floor ties; 80 with floor-tied minima')
    text(root,26,113,'Disjoint minimum coordinates: 2680 / 3224 (83.1%); not a count of independent signals.')

    heading(root,'b',133,'chr19: reuse of one minimum-P input')
    text(root,26,147,'TM6SF2 E167K (rs58542926) | GRCh37 chr19:19379549 C/T')
    text(root,26,160,'50 candidates | 154 observed support events across six traits')
    matrix=read(source/inputs[3]);genes=[k for k in matrix[0] if k!='trait']
    traits=['CIRR','LIVER_FAT','NAFLD','ALT','ALP','BILI']
    labels=['Cirrhosis','Liver fat','NAFLD*','ALT','ALP','Bilirubin']
    values={r['trait']:r for r in matrix}
    assert len(genes)==50 and len(matrix)==6
    for offset,subset in [(0,genes[:25]),(244,genes[25:])]:
        x,y,cw,rh=87+offset,196,27,8.4
        for j,label in enumerate(labels):
            tx=x+j*cw+16;ty=y-7
            text(root,tx,ty,label,transform=f'rotate(-40 {tx} {ty})')
        for i,gene in enumerate(subset):
            text(root,x-4,y+i*rh+6.5,gene,anchor='end')
            for j,trait in enumerate(traits):
                value=int(values[trait][gene]);assert value in (0,1)
                node(root,'rect',x=x+j*cw,y=y+i*rh,width=cw,height=rh,fill=BLUE if value else GRAY,stroke='white',stroke_width=1,
                     data_gene=gene,data_trait=trait,data_value=value)
    node(root,'rect',x=26,y=414,width=8,height=7,fill=BLUE)
    text(root,38,420,'Verified minimum input in this gene–trait event')
    node(root,'rect',x=294,y=414,width=8,height=7,fill=GRAY)
    text(root,306,420,'No such event for this variant')

    heading(root,'c',443,'chr6: trait-specific and shared strongest inputs')
    text(root,26,457,'GRCh37 chr6:25684587–26791233 | all 10 verified minimum inputs')
    matrix6=read(source/inputs[4]);assert len(matrix6)==10
    x,y,cw,rh=191,482,74,10
    for j,label in enumerate(['Cirrhosis','Bilirubin','AIH','PSC']):text(root,x+(j+.5)*cw,475,label,7.5,anchor='middle')
    for i,row in enumerate(matrix6):
        key=row['variant_key'];parts=key.split(':')
        text(root,x-6,y+i*rh+7.4,parts[2]+' '+parts[3]+'/'+parts[4],7.5,anchor='end')
        for j,trait in enumerate(['CIRR','BILI','AIH','PSC']):
            n=int(row[trait])
            node(root,'rect',x=x+j*cw,y=y+i*rh,width=cw,height=rh,fill=TEAL if n else GRAY,stroke='white',stroke_width=1,
                 data_variant=key,data_trait=trait,data_count=n)
            if n:text(root,x+(j+.5)*cw,y+i*rh+7.5,n,7.5,color='white',anchor='middle')
    text(root,26,598,'Numbers: distinct candidate genes using each variant as a minimum input for that trait.')
    text(root,26,611,'Blank cells: no observed event. Different inputs may still tag the same genetic signal.')
    save(root,3,out)
    assert sum(int(values[t][g]) for t in traits for g in genes)==154
    assert sum(int(row[t]) for row in matrix6 for t in ['CIRR','BILI','AIH','PSC'])==194
    (out/'render_record.json').write_text(json.dumps({'source_sha256':hashes,'model_rows':model_records,
        'figure_sizes_pt':[[510,346],[510,262],[510,620]],'models_run':0,
        'chr19_genes':50,'chr19_events':154,'chr6_variants':10,'chr6_events':194},indent=2)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parents[1]/'figure_source_data')
    parser.add_argument('--output-dir',type=Path,default=Path('figures'))
    args=parser.parse_args();render(args.source_dir,args.output_dir)
