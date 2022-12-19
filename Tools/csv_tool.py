def read_csv(path: str, skip_first: bool = True):
    data = []
    with open(path) as csv:
        raw_data = csv.readlines()
        for row in raw_data[skip_first:]:
            temp = row.replace('\n', '').split(',')
            data.append(temp)
    return data
