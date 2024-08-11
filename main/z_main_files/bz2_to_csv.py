import sys
from pyspark.sql import SparkSession
from pyspark.sql.functions import row_number
from pyspark.sql.window import Window

input_file = sys.argv[1]
output_file = sys.argv[2]

spark = SparkSession.builder.appName("FormatConverterSingleFile") \
                            .master("local[5]") \
                            .getOrCreate()

df = spark.read.option("header", "true") \
               .option("sep", "\t") \
               .option("compression", "bzip2") \
               .csv(input_file) \
               .select('id', 'smiles')

windowSpec = Window.orderBy('id')

df = df.withColumn('db_index', row_number().over(windowSpec) - 1) \
       .select('db_index', 'id', 'smiles')

df.coalesce(1).write.mode("overwrite") \
                    .option("header", "true") \
                    .option("sep", ",") \
                    .option("compression", "gzip") \
                    .csv(output_file)

spark.stop()
