

def parse_anderson(cfg, gene_index, lift, proximity_limit):
    """
    Parse the anderson file and return a dataframe with the intervals.

    :param anderson_file: Path to the anderson file (.bed).
    :return: Spark Dataframe

    **Summary of the logic**

    - Reading .bed file (input)
    - Parsing the names column -> chr, start, end, gene, score
    - Mapping the coordinates to the new build -> liftover
    - Joining with target index by gene symbol (some loss as input uses obsoleted terms)
    - Dropping rows where the gene is on other chromosomes
    - Dropping rows where the gene TSS is too far from the midpoint of the intervals
    - Adding constant columns for this dataset
    - Return spark dataframe.
    """

    # Read the anderson file:
    anderson_df = (
        SparkSession.getActiveSession().createDataFrame(pd.read_csv(cfg.data_file, sep='\t', header=0, low_memory=False, skiprows=1))
        # Parsing score column and casting as float:
        .withColumn('score', f.col('score').cast('float') / f.lit(1000))

        # Parsing the 'name' column:
        .withColumn('parsedName', f.split(f.col('name'), ';'))
        .withColumn('gene_symbol', f.col('parsedName')[2])
        .withColumn('location', f.col('parsedName')[0])
        .withColumn('chrom', f.regexp_replace(f.split(f.col('location'), ':|-')[0], 'chr', ''))
        .withColumn('start', f.split(f.col('location'), ':|-')[1].cast(t.IntegerType()))
        .withColumn('end', f.split(f.col('location'), ':|-')[2].cast(t.IntegerType()))

        # Select relevant columns:
        .select('chrom', 'start', 'end', 'gene_symbol', 'score')

        # Drop rows with non-canonical chromosomes:
        .filter(f.col('chrom').isin([str(x) for x in range(1, 23)] + ['X', 'Y', 'MT']))

        # For each region/gene, keep only one row with the highest score:
        .groupBy('chrom', 'start', 'end', 'gene_symbol')
        .agg(f.max('score').alias('score'))

        .orderBy('chrom', 'start')
        .persist()
    )

    # Prepare gene set:
    genes = gene_index.withColumnRenamed('gene_name', 'gene_symbol').select('gene_symbol', 'chr', 'gene_id', 'TSS')

    return (
        # Lift over the intervals:
        lift.convert_intervals(anderson_df, 'chrom', 'start', 'end')
        .drop('start', 'end')
        .withColumnRenamed('mapped_start', 'start')
        .withColumnRenamed('mapped_end', 'end')
        .distinct()

        # Joining with the gene index (unfortunately we are losing a bunch of genes here due to old symbols):
        .join(genes, on='gene_symbol', how='left')
        .filter(
            # Drop rows where the gene is not on the same chromosome
            (f.col('chrom') == f.regexp_replace(f.col('chr'), 'chr', ''))
            # Drop rows where the TSS is far from the start of the region
            & (
                f.abs(
                    (f.col('start') + f.col('end')) / 2 - f.col('TSS')
                ) <= proximity_limit
            )
        )

        # Adding constant values:
        .withColumn('dataset_name', f.lit(cfg.dataset_name))
        .withColumn('data_type', f.lit(cfg.data_type))
        .withColumn('experiment_type', f.lit(cfg.experiment_type))
        .withColumn('pmid', f.lit(cfg.pubmed_id))
        .withColumn('bio_feature', f.lit(cfg.bio_feature))

        # Select relevant columns:
        .select(
            'chrom', 'start', 'end', 'gene_id', 'score', 'dataset_name', 'data_type', 'experiment_type', 'pmid', 'bio_feature'
        )
        .persist()
    )
