from pickle import dump

from assign0.code import load_task, naive, fuzzy_location_search


if __name__ == "__main__":
    _, path_2, sequence = load_task()
    locations, count = fuzzy_location_search(file_path=path_2, pattern=sequence, method=naive,
                                             length=len(sequence), stride=1, tolerance=1)
    with open("./results/CHM13v2.naive.1.pkl", "wb") as file:
        dump(obj=(locations, count), file=file)

