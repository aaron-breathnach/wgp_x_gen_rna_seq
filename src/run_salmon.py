import os
from pathlib import Path

SALMON = '''
    salmon quant 
    -i salmon/mouse_index 
    -l A 
    -1 reads/processed/{prefix}_1.qcd.fq.gz 
    -2 reads/processed/{prefix}_2.qcd.fq.gz 
    -p {threads} 
    --validateMappings 
    -o quants/{prefix}_quant
'''.replace('\n', '').replace('    ', '')

salmon = {'cmd': SALMON, 'target': 'quants/{prefix}_quant/quant.sf'}

with open('run_salmon.sh', 'w') as f:
    for x in Path('reads/raw').rglob('*_1.fq.gz'):
        prefix = x.as_posix().split('/')[-1:][0].split('_1.fq.gz')[0]
        target = salmon['target'].format(prefix=prefix)
        if not os.path.exists(target):
            cmd = salmon['cmd'].format(prefix=prefix, threads=os.cpu_count())
            f.write(cmd + '\n')
