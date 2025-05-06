from argparse import ArgumentParser


def evaluate_assembly(number):
    data = []
    for index in range(1, number + 1):
        with open("./tmp/report_" + str(index) + ".txt") as file:
            values = [[], []]
            for line_index, line in enumerate(file.readlines()):
                if line_index > 1:
                    infos = line.strip().split("  ")
                    values[0].append(infos[0])
                    values[1].append(infos[-1])
        print(len(values[1]))
        data.append(values)

    for line_index in range(30):
        string = "| " + data[0][0][line_index] + " | " + data[0][1][line_index]
        for index in range(1, number):
            string += " | " + data[index][1][line_index]
        string += " |"
        print(string)

        if line_index == 0:
            print(("|--" * number)[:-2])


if __name__ == "__main__":
    parser = ArgumentParser(description="QUAST merge")
    parser.add_argument("-n", "--number", required=True, help="number of QUAST report")
