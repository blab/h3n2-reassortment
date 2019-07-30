"""
This script runs the entire suite of seasonal flu builds via AWS batch
One batch run is created per HA/NA combination, ie
h3n2_ha_2y + h3n2_na_2y is one build and
h3n2_ha_6y + h3n2_na_6y is another
"""

import subprocess
import argparse
import os

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run h3n2-reassortment flu builds on aws', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--system', type = str, default = 'local', help='where to run, local or batch')
    parser.add_argument('-s', '--segments', nargs='+', type = str, help ="flu segments to include", default=['ha', 'na', 'mp', 'pb1', 'pb2', 'pa', 'ns', 'np'])
    params = parser.parse_args()

    if not os.path.exists("logs"):
        os.makedirs("logs")
    else:
        rmc = 'rm -rf logs/*'
        subprocess.call(rmc, shell=True)

    lineage = 'h3n2'
    resolution = '2y'

    cpus = len(params.segments)
    memory = 1800 * cpus
    if params.system == 'local':
        call = ['nextstrain', 'build', '--native', '.', '--jobs', '1']
    elif params.system == 'batch':
        call = ['nextstrain', 'build', '--aws-batch', '--aws-batch-cpus', str(cpus), '--aws-batch-memory', str(memory), '.', '--jobs', str(cpus)]
    targets = []
    for segment in params.segments:
        targets.append('targets/flu_seasonal_%s_%s_%s'%(lineage, segment, resolution))
    call.extend(targets)
    print(' '.join(call))
    log = open('logs/live_%s.txt'%(lineage), 'w')
    if params.system == 'local':
        pro = subprocess.call(call)
    if params.system == 'batch':
        pro = subprocess.Popen(call, stdout=log, stderr=log)
