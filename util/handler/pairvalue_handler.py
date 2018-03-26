"""
Pair-value handler

@auth: Yu-Hsiang Fu
@date: 2014/09/27
@update: 2018/03/26
"""


def read_pairvalue(file_path, is_int=False):
    pair_value = {}

    try:
        with open(file_path, mode="r") as f:
            for line in f:
                pair = line.strip().split()
                key = int(pair[0])

                if len(pair) == 2:
                    pair_value[key] = int(pair[1]) if is_int else float(int(pair[1]))
                else:
                    pair_value[key] = tuple([int(p) if is_int else float(p) for p in pair[1:]])

            f.close()
    except:
        print('[Error] The file can not be read ...')
        print('[Error] Please check this:  ' + str(file_path))

    return pair_value


def write_pairvalue(pair_value, file_path):
    try:
        with open(file_path, mode="w") as f:
            for key in pair_value.keys():
                value = " ".join([str(p) for p in pair_value[key]])
                f.write("{0} {1}\n".format(key, value))

            f.flush()
            f.close()
    except:
        print('[Error] The file can not be written ...')
        print('[Error] Please check this: ' + str(file_path))
