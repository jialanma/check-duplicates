# check_duplicates

The pipeline generates read counts and QC reports and marks the duplicated reads in RNA-seq data.

The pipeline runs on Hail Batch on GCP.

How to run:

`python3 main.py --billing-project [your GCP billing project] --requester-pays-project [your GCP project to pay for access] --file-dir [a temporary GCP directory] 
--tissue [tissue type] --out-dir [GCP output patb] [IDs of samples to process]`
