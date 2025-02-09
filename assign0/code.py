from os import path, makedirs, remove
from numpy import ndarray, array, zeros, min, ceil
from requests import get
from typing import Generator, Tuple
from zipfile import ZipFile


def download_genome(url: str,
                    folder: str,
                    file_name: str) \
        -> str:
    """
    Download genome data from given url and save it in given folder.

    :param url: genome data url using NCBI API.
    :type url: str

    :param folder: folder to save fasta file.
    :type folder: str

    :param file_name: file name of FASTA file.
    :type file_name: str

    :return: file path.
    :rtype: str
    """
    makedirs(folder, exist_ok=True)

    file_path = path.join(folder, file_name)
    if not path.exists(file_path):
        response = get(url, stream=True)
        if response.status_code == 200:
            with open(path.join(folder, "temp.zip"), "wb") as file:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    file.write(chunk)

            print("Downloaded: \"" + file_name + "\".")

            with ZipFile(path.join(folder, "temp.zip")) as zip_file:
                for target_file in zip_file.namelist():
                    if target_file.endswith(".fna"):
                        with zip_file.open(target_file) as source, open(file_path, "wb") as target:
                            target.write(source.read())

            print("Extracted: \"" + file_name + "\".")
            remove(path.join(folder, "temp.zip"))

        else:
            print("Failed to download \"" + file_name + "\" with status code " + str(response.status_code) + ".")
    else:
        print("Downloaded and extracted: \"" + file_name + "\".")

    return file_path


def load_genome(file_path: str,
                batch_length: int,
                batch_stride: int) \
        -> Generator[Tuple[int, int, str], None, None]:
    """
    Yield k-mer genome data from given file path.

    :param file_path: FASTA genome data file path.
    :type file_path: str

    :param batch_length: length of k-mer in the genome.
    :type batch_length: int

    :param batch_stride: stride length of k-mer in the genome.
    :type batch_stride: int

    :return: k-mer genome data generator.
    :rtype: Generator
    """
    saved_segment, sequence_index, sequence_location = "", -1, 0
    with open(file_path) as handle:
        for line in handle.readlines():
            line = line.strip().upper()

            if line.startswith(">"):
                sequence_index += 1

                # yield remaining k-mers from the previous sequence.
                while len(saved_segment) >= batch_length:
                    yield sequence_index, sequence_location, saved_segment[:batch_length]
                    if len(saved_segment) >= batch_stride + batch_length:
                        sequence_location += batch_stride
                        saved_segment = saved_segment[batch_stride:]
                    else:
                        sequence_location += len(saved_segment) - batch_length
                        saved_segment = saved_segment[len(saved_segment) - batch_length:]

                saved_segment = ""

            else:
                saved_segment += line.upper()
                while len(saved_segment) >= batch_length:
                    yield sequence_index, sequence_location, saved_segment[:batch_length]
                    if len(saved_segment) >= batch_stride + batch_length:
                        sequence_location += batch_stride
                        saved_segment = saved_segment[batch_stride:]
                    else:
                        break

    # yield any remaining k-mers at the end of the file.
    while len(saved_segment) >= batch_length:
        yield sequence_index, sequence_location, saved_segment[:batch_length]
        if len(saved_segment) >= batch_stride + batch_length:
            sequence_location += batch_stride
            saved_segment = saved_segment[batch_stride:]
        else:
            sequence_location += len(saved_segment) - batch_length
            saved_segment = saved_segment[len(saved_segment) - batch_length:]


def exact_location_search(file_path: str,
                          pattern: str,
                          method: callable,
                          length: int,
                          stride: int) \
        -> Tuple[ndarray, int]:
    """
    Exact search for the location of the template that can match the pattern.

    :param file_path: FASTA file path of template.
    :type file_path: str

    :param pattern: a sequence pattern used for detection.
    :type pattern: str

    :param method: used matching method.
    :type method: callable.

    :param length: length of k-mer in the genome.
    :type length: int

    :param stride: stride length of k-mer in the genome.
    :type stride: int

    :return: matched locations and the matching count.
    :rtype: numpy.ndarray, int
    """
    matched_locations = set()
    for index, location, segment in load_genome(file_path=file_path, batch_length=length, batch_stride=stride):
        matched_local_positions = method(template=segment, detected=pattern)
        for local_position in matched_local_positions:
            matched_locations.add((index, location + local_position, location + local_position + len(pattern)))

    return array(list(matched_locations)), len(matched_locations)


def fuzzy_location_search(file_path: str,
                          pattern: str,
                          method: callable,
                          length: int,
                          stride: int,
                          tolerance: int) \
        -> Tuple[ndarray, int]:
    """
    Fuzzy search for the location of the template that can match the pattern.

    :param file_path: FASTA file path of template.
    :type file_path: str

    :param pattern: a sequence pattern used for detection.
    :type pattern: str

    :param method: used matching method.
    :type method: callable.

    :param length: length of k-mer in the genome.
    :type length: int

    :param stride: stride length of k-mer in the genome.
    :type stride: int

    :param tolerance: maximum tolerable number of edit errors.
    :type tolerance: int

    :return: matched locations and the matching count.
    :rtype: numpy.ndarray, int
    """
    assert tolerance > 0

    matched_locations = set()
    for index, location, segment in load_genome(file_path=file_path, batch_length=length, batch_stride=stride):
        matched_local_positions = method(template=segment, detected=pattern, tolerance=tolerance)
        for local_position in matched_local_positions:
            matched_locations.add((index, location + local_position, location + local_position + len(pattern)))

    return array(list(matched_locations)), len(matched_locations)


def naive(template: str,
          detected: str) \
        -> list:
    """
    Calculate the mapping between template and detected sequences by the naive algorithm.

    :param template: template sequence.
    :type template: str

    :param detected: detected sequence.
    :type detected: str

    :return: matched locations.
    :rtype: list
    """
    if len(template) == len(detected):
        for info_1, info_2 in zip(template, detected):
            if info_1 != info_2:
                return []
        return [0]

    elif len(template) > len(detected):
        matches = []
        for location in range(len(template) - len(detected) + 1):
            flag = True
            for index in range(len(detected)):
                if template[location + index] != detected[index]:
                    flag = False
                    break

            if flag:
                matches.append(location)

        return matches

    else:
        raise ValueError("The length of the template sequence should not be less than that of the detected sequence!")


def kmp(template: str,
        detected: str) \
        -> list:
    """
    Calculate the mapping between template and detected sequences by the Knuth–Morris–Pratt algorithm.

    :param template: template sequence.
    :type template: str

    :param detected: detected sequence.
    :type detected: str

    :return: matched locations.
    :rtype: list
    """
    if len(template) == len(detected):
        for info_1, info_2 in zip(template, detected):
            if info_1 != info_2:
                return []
        return [0]

    elif len(template) > len(detected):
        longest_prefix_suffix, length, location = zeros(shape=(len(detected),), dtype=int), 0, 1
        while location < len(detected):
            if detected[location] == detected[length]:
                length += 1
                longest_prefix_suffix[location] = length
                location += 1
            else:
                if length != 0:
                    length = longest_prefix_suffix[length - 1]
                else:
                    longest_prefix_suffix[location] = 0
                    location += 1

        index_1, index_2, matches = 0, 0, []

        while index_1 < len(template):
            if template[index_1] == detected[index_2]:
                index_1, index_2 = index_1 + 1, index_2 + 1

                if index_2 == len(detected):
                    matches.append(index_1 - index_2)
                    index_2 = longest_prefix_suffix[index_2 - 1]
            else:
                if index_2 != 0:
                    index_2 = longest_prefix_suffix[index_2 - 1]
                else:
                    index_1 += 1

        return matches

    else:
        raise ValueError("The length of the template sequence should not be less than that of the detected sequence!")


def q_gram(template: str,
           detected: str,
           tolerance: int) \
        -> list:
    """
    Calculate the mapping between template and detected sequences by the q-gram lemma algorithm.

    :param template: template sequence.
    :type template: str

    :param detected: detected sequence.
    :type detected: str

    :param tolerance: maximum tolerable number of edit errors.
    :type tolerance: int

    :return: matched locations.
    :rtype: list
    """
    if tolerance == 0:
        raise ValueError("Error tolerance needs to be greater than 0, otherwise please use the exact location search!")

    best_q = ceil(len(detected) / (tolerance + 1)).astype(int)

    if len(template) >= len(detected):
        detected_q_gram_table = {}
        for index in range(len(detected) - best_q + 1):
            q_mer = detected[index: index + best_q]
            if q_mer in detected_q_gram_table:
                detected_q_gram_table[q_mer] += 1
            else:
                detected_q_gram_table[q_mer] = 1

        matches = []
        for location in range(len(template) - len(detected) + 1):
            template_q_gram_table = {}
            for index in range(len(detected) - best_q + 1):
                q_mer = template[location + index: location + index + best_q]
                if q_mer in template_q_gram_table:
                    template_q_gram_table[q_mer] += 1
                else:
                    template_q_gram_table[q_mer] = 1

            overlap = 0
            for q_mer, count in detected_q_gram_table.items():
                if q_mer in template_q_gram_table:
                    overlap += min([count, template_q_gram_table[q_mer]])

            if overlap >= len(detected) - best_q + 1 - tolerance * best_q:
                matches.append(location)

        return matches

    else:
        raise ValueError("The length of the template sequence should not be less than that of the detected sequence!")


def load_task() \
        -> Tuple[str, str, str]:
    address = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
    params = "/download?include_annotation_type=GENOME_FASTA"

    file_path_1 = download_genome(url=address + "GCF_000001405.39" + params, folder="./genomes/", file_name="GRCh38.fa")
    file_path_2 = download_genome(url=address + "GCF_009914755.1" + params, folder="./genomes/", file_name="CHM13v2.fa")

    aluy_sequence = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGA" + \
                    "TCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAA" + \
                    "AAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGG" + \
                    "CTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGC" + \
                    "CACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAA" + \
                    "AAAAAAAAAAA"  # given from the assignment 0.

    return file_path_1, file_path_2, aluy_sequence
