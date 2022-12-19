from Tools import console_appearance as ca

length = 30


def loading_bar(x: float, x_max: float):
    pr = int(x * 100 // x_max)
    bar = '█' * int(length * x // x_max) + '-' * int((length - (length * x // x_max) - 1))
    if pr == 100:
        color = ca.GREEN
        end = f'Completed {ca.END}\n'
    else:
        color = ca.RED
        end = f'{pr}% - {x}/{x_max}'
    print(f'\r{color}{bar} {end}', end='')
