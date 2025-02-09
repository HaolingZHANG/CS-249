from pickle import dump

from assign0.code import load_task, q_gram, fuzzy_location_search


if __name__ == "__main__":
    path_1, _, sequence = load_task()
    locations, count = fuzzy_location_search(file_path=path_1, pattern=sequence, method=q_gram,
                                             length=len(sequence), stride=1, tolerance=1)
    with open("./results/GRCh38.q-gram.1.pkl", "wb") as file:
        dump(obj=(locations, count), file=file)

