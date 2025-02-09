from pickle import dump

from assign0.code import load_task, kmp, exact_location_search


if __name__ == "__main__":
    _, path_2, sequence = load_task()
    locations, count = exact_location_search(file_path=path_2, pattern=sequence, method=kmp,
                                             length=len(sequence) * 10, stride=len(sequence) * 8)
    with open("/hy-tmp/results/CHM13v2.kmp.0.pkl", "wb") as file:
        dump(obj=(locations, count), file=file)
