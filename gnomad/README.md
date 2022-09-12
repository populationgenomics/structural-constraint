# gnomad V4 power analysis

* generate a site-only hail Table from gnomADv4_0
* repartition
* VEP 101 annotate
* extract SFS of missense variation to 10 selected transcripts

## Workflow:

repartition step:

```sh
analysis-runner --dataset constraint --output-dir "."\
    --access-level test 
    --description "repartition 10 transcripts site-onlt Table" \
    python3 10transcripts_repartition.py
```
