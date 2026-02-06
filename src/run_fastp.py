import os
from pathlib import Path

FASTP = '''
    docker run -v $PWD:/data staphb/fastp 
    fastp 
    -i reads/raw/{prefix}_1.fq.gz 
    -I reads/raw/{prefix}_2.fq.gz 
    -o reads/processed/{prefix}_1.qcd.fq.gz 
    -O reads/processed/{prefix}_2.qcd.fq.gz 
    -j reads/processed/{prefix}.json 
    -h reads/processed/{prefix}.html
'''.replace('\n', '').replace('    ', '')

fastp = {'cmd': FASTP, 'target': 'reads/processed/{prefix}_1.qcd.fq.gz'}

with open('run_fastp.sh', 'w') as f:
    if not os.path.exists('reads/processed'):
        os.makedirs('reads/processed')
    for x in Path('reads/raw').rglob('*_1.fq.gz'):
        prefix = x.as_posix().split('/')[-1:][0].split('_1.fq.gz')[0]
        target = fastp['target'].format(prefix=prefix)
        if not os.path.exists(target):
            cmd = fastp['cmd'].format(prefix=prefix)
            f.write(cmd + '\n')

