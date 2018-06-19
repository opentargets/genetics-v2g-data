Get cis-regulatory datasets
===========================

Process cis-regulatory datasets and upload them to 'genetics-portal-data
' bucket.

#### Datasets

### Interval datasets
- PCHiC from Javierre et al 2016
  - [Publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123897/)
  - [Dataset](ftp://ftp.ebi.ac.uk/pub/contrib/pchic/CHiCAGO/))
- Fantom5 enhancer-TSS associations from Andersson et al 2014
  - [Publication](https://www.nature.com/articles/nature12787#linking-enhancer-usage-with-tss-expression)
  - [Dataset](http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed)
- DHS correlation with gene promoters from Thurman et al 2012
  - [Publication](https://www.nature.com/articles/nature11232#a-map-of-distal-dhstopromoter-connections)
  - [Dataset](http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz)

Output file naming convention:
  `{data type}/{experiment type}/{source}/{cell|tissue type}/{chrom}.{processed|raw}.tsv.gz`

Interval output columns:
  - `chrom`, `start`, `end`: 1-based genomic co-ordinates. Start and end are inclusive.
  - `ensembl_id`: Ensembl gene ID
  - `score`: feature score from original dataset (unmodified)
  - `cell_type`: cell or tissue type of this experiment

### QTL datasets
- GTEx v7
  - TODO
- pQTL from Sun et al 2018
  - TODO
  - [Publication](https://www.nature.com/articles/s41586-018-0175-2)
  - [Data available from](http://www.phpc.cam.ac.uk/ceu/proteins/)
  - Data downloaded to `gs://genetics-portal-data/pqtl_sun2018`
  - Notes:
    - Do we only keep cis-associations?
    - Do we need to exclude "pQTLs without evidence against binding effects"
      - Main body of paper: "Of 1,927 pQTLs, 549 (28%) were cis-acting (Supplementary Table 4). Genetic variants that change protein structure may result in apparent cis pQTLs owing to altered aptamer binding rather than true quantitative differences in protein levels. We found evidence against such artefactual associations for 371 (68%) cis pQTLs (Methods, Supplementary Tables 4, 7, 8). The results were materially unchanged when we repeated downstream analyses but excluded pQTLs without evidence against binding effects."
      - Methods section title "Evidence against aptamer-binding effects at cis pQTLs".
    - What significance threshold should we use?
    - What to do when there is cross-reactivity?
    - Some genes have multiple SOMA IDs (assayed more than once)

Output file naming convention:
  `{data_type}/{exp_type}/{source}/{cell_type}/{ensembl_id}/{chrom}.pval{pval}.{cis_reg|trans_reg|all}.tsv.gz`

Output files:
  - `{chrom}.pval{pval}.cis_reg.tsv.gz`: cis-regulatory associations for `{ensembl_id}` within 1Mb (`config['sun2018_cis_window']`) of gene TSS. Filtered by pval < 2.5e-5 (`config['sun2018_cis_pval']`).
  - `{chrom}.pval{pval}.trans_reg.tsv.gz`: trans-regulatory associations for `{ensembl_id}` outside 1Mb (`config['sun2018_cis_window']`) of gene TSS. Filtered by pval < 5e-8 (`config['sun2018_trans_pval']`).
  - `{chrom}.pval1.all.tsv.gz`: all associations for `{ensembl_id}` in common format.

QTL output columns:
  - `chrom`, `pos`: genomic coordinates on GRCh37
  - `other_allele`: non-effect allele
  - `effect_allele`: allele on which the effect (`beta`) is measured
  - `beta`: the effect size and direction
  - `se`: standard error of `beta`
  - `pval`: p-value of association

#### Usage

```
# Install dependencies into isolated environment
conda env create -n get_cisreg_data --file environment.yaml

# Activate environment
source activate get_cisreg_data

# Alter configuration file
nano config.yaml

# Authenticate google cloud storage
gcloud auth application-default login

# Download and annotate manifest for Sun pQTL dataset
snakemake -s scripts/sun2018_pqtl.make_manifest.Snakefile

# Execute workflow (locally)
snakemake
```

#### Notes

I am getting the following error when trying to access google storage in parallel:

```
Traceback (most recent call last):
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/snakemake/executors.py", line 344, in _callback
    callback(job)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/snakemake/scheduler.py", line 323, in _proceed
    self.get_executor(job).handle_job_success(job)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/snakemake/executors.py", line 357, in handle_job_success
    super().handle_job_success(job)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/snakemake/executors.py", line 151, in handle_job_success
    assume_shared_fs=self.assume_shared_fs)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/snakemake/jobs.py", line 843, in postprocess
    self.dag.handle_remote(self, upload=upload_remote)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/snakemake/dag.py", line 489, in handle_remote
    f.upload_to_remote()
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/snakemake/io.py", line 322, in upload_to_remote
    self.remote_object.upload()
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/snakemake/remote/GS.py", line 88, in upload
    if not self.bucket.exists():
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/google/cloud/storage/bucket.py", line 171, in exists
    query_params=query_params, _target_object=None)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/google/cloud/_http.py", line 299, in api_request
    headers=headers, target_object=_target_object)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/google/cloud/_http.py", line 193, in _make_request
    return self._do_request(method, url, headers, data, target_object)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/google/cloud/_http.py", line 223, in _do_request
    body=data)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/google_auth_httplib2.py", line 198, in request
    uri, method, body=body, headers=request_headers, **kwargs)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/httplib2/__init__.py", line 1509, in request
    (response, content) = self._request(conn, authority, uri, request_uri, method, body, headers, redirections, cachekey)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/httplib2/__init__.py", line 1259, in _request
    (response, content) = self._conn_request(conn, request_uri, method, body, headers)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/site-packages/httplib2/__init__.py", line 1212, in _conn_request
    response = conn.getresponse()
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/http/client.py", line 1331, in getresponse
    response.begin()
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/http/client.py", line 297, in begin
    version, status, reason = self._read_status()
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/http/client.py", line 258, in _read_status
    line = str(self.fp.readline(_MAXLINE + 1), "iso-8859-1")
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/socket.py", line 586, in readinto
    return self._sock.recv_into(b)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/ssl.py", line 1009, in recv_into
    return self.read(nbytes, buffer)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/ssl.py", line 871, in read
    return self._sslobj.read(len, buffer)
  File "/Users/em21/miniconda3/envs/get_cisreg_data/lib/python3.6/ssl.py", line 631, in read
    v = self._sslobj.read(len, buffer)
ssl.SSLError: [SSL: WRONG_VERSION_NUMBER] wrong version number (_ssl.c:2273)
```
