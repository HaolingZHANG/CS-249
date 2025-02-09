from pickle import dump

from assign0.code import load_task, naive, exact_location_search


if __name__ == "__main__":
    _, path_2, sequence = load_task()
    locations, count = exact_location_search(file_path=path_2, pattern=sequence, method=naive,
                                             length=len(sequence), stride=1)
    with open("/hy-tmp/results/CHM13v2.naive.0.pkl", "wb") as file:
        dump(obj=(locations, count), file=file)
