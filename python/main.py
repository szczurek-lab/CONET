import conet.src.data_converter.data_converter as dc
import pandas

x = dc.DataConverter("notebooks/data/SA501X3F_filtered_corrected_counts.csv", ',', 150000, 3095677412, True)
y = x.corrected_counts

df = pandas.read_csv('indeksy', header=None)[0].tolist()
x.create_CoNET_input_files(df, "./")
print("dupa")