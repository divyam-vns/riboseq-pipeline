#!/usr/bin/env python3
"""
get_srr.py

Uses Biopython Entrez to retrieve SRR run accessions for given GEO GSE IDs.
Outputs a simple list (one SRR per line).
"""
import argparse
from Bio import Entrez
import time


def gse_to_srr(gse, email):
    Entrez.email = email
    # search GEO GDS db for the GSE
    handle = Entrez.esearch(db="gds", term=gse)
    rec = Entrez.read(handle)
    handle.close()
    if not rec['IdList']:
        return []
    gds_id = rec['IdList'][0]
    # link to sra
    handle = Entrez.elink(dbfrom='gds', db='sra', id=gds_id)
    linkrec = Entrez.read(handle)
    handle.close()
    srrs = []
    if linkrec and linkrec[0].get('LinkSetDb'):
        links = linkrec[0]['LinkSetDb'][0]['Link']
        for link in links:
            sra_id = link['Id']
            # fetch runinfo
            ef = Entrez.efetch(db='sra', id=sra_id, rettype='runinfo', retmode='text')
            runinfo = ef.read()
            ef.close()
            # parse Run column lines
            for line in runinfo.splitlines():
                if line.startswith('Run,'):
                    continue
                parts = line.split(',')
                if parts:
                    run = parts[0]
                    if run.startswith('SRR'):
                        srrs.append(run)
            time.sleep(0.4)
    return sorted(list(set(srrs)))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--geo', nargs='+', required=True, help='GEO GSE ids')
    p.add_argument('--email', required=True, help='Entrez email')
    p.add_argument('--out', required=True, help='output file for SRR list')
    args = p.parse_args()

    all_srr = []
    for g in args.geo:
        print(f'Fetching SRR for {g}...')
        s = gse_to_srr(g, args.email)
        print(f'  found {len(s)} runs')
        all_srr.extend(s)

    all_srr = sorted(list(set(all_srr)))
    with open(args.out, 'w') as fh:
        for r in all_srr:
            fh.write(r + '\n')
    print(f'Wrote {len(all_srr)} SRR ids to {args.out}')
