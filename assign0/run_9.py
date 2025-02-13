from pickle import load


if __name__ == "__main__":
    with open("./results/GRCh38.naive.0.pkl", "rb") as file:
        locations, count = load(file=file)
        print(locations, count)

    with open("./results/GRCh38.kmp.0.pkl", "rb") as file:
        locations, count = load(file=file)
        print(locations, count)

    with open("./results/CHM13v2.naive.0.pkl", "rb") as file:
        locations, count = load(file=file)
        print(locations, count)

    with open("./results/CHM13v2.kmp.0.pkl", "rb") as file:
        locations, count = load(file=file)
        print(locations, count)
