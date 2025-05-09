from pandas import read_csv
from matplotlib import pyplot

# Load spectra-cn.hist file
file_path = "./results/hifiasm/merqury[spectra-cn].hist"
dataframe = read_csv(file_path, sep="\t", header=None, names=["Copy Number", "Count"])

# Plot the k-mer spectrum
pyplot.figure(figsize=(10, 6), tight_layout=True)
pyplot.scatter(dataframe["Copy Number"], dataframe["Count"], marker='o')
pyplot.title("K-mer Copy Number Spectrum")
pyplot.xlabel("Copy Number")
pyplot.ylabel("K-mer Count")
pyplot.xscale("log")
pyplot.yscale("log")
pyplot.savefig("./results/images/distribution.png")
pyplot.close()
