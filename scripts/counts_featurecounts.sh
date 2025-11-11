#!/usr/bin/env bash
set -euo pipefail

GTF=${1}
shift
BAMS="$@"

featureCounts -T 8 -t CDS -g gene_id -a $GTF -o counts.txt $BAMS

# output: counts.txt
