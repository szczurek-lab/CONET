

import synthetic_data.difference_distribution
import CoNET.conet_parameters as p

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
    #x = synthetic_data.difference_distribution.DifferenceDistribution([12, 3], [-2, 3], [1, 2], 8)
    y = p.CoNETParameters()
    print(y.to_string())
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
