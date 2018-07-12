Testing CWL remote files
========================

Simple `cp` CWL tool to test using remote files as input.

In the [CWL specification](https://www.commonwl.org/v1.0/Workflow.html#File) the `location` property of `class` `File` allows for different URI schemes to be used.

This works for local (file://) and http:// inputs, but not for other remote inputs (e.g. ftp:// or gs://). I have also tried using `toil-cwl-runner` in place of the reference cwl runner but this doesn't work either.

There is a [CWL runner script from Google Genomics](https://github.com/googlegenomics/pipelines-api-examples/tree/master/cwl_runner) which runs a workflow on Google Cloud Platform. However, this works by creating a disk, starting a compute engine, then copying all inputs from a google bucket before running the workflow.

```
# Run using a local input (works)
cwl-runner cp.cwl cp_from_local-job.yml

# Run using http file as input (works)
cwl-runner cp.cwl cp_from_http-job.yml

# Run using ftp file as input (doesn't work)
cwl-runner cp.cwl cp_from_ftp-job.yml

# Run using GCS as input (doesn't work)
cwl-runner cp.cwl cp_from_gcs-job.yml
```
