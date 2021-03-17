
import conet.src.synthetic_data.ratios_distribution as rd
import conet.src.synthetic_data.data_generator as dg

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
    data_generator = dg.DataGenerator()
    ratios_distribution = rd.RatiosDistribution([0.2, 0.5, 0.3], [1.0, 3, 2.5], [0.1, 0.6, 0.3], 0.3)
    tree = data_generator.generate_random_tree(40, 20)
    input_data = data_generator.generate_random_ratios(40, tree, ratios_distribution, 80)
    input_data.save("./")
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
