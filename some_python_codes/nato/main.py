"""
high level support for doing this and that.
"""
# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import pandas
#import pylint


# Myfunc
def print_hi(name):
    """A dummy docstring."""
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


def generate():
    """A dummy docstring."""
    name = input("Adon ne?").upper()
    try:

        mylist = [nato_dic[c] for c in name]

    except KeyError:
        print("Invalid name")
    else:
        print(mylist)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

    # See PyCharm help at https://www.jetbrains.com/help/pycharm/
    df = pandas.read_csv('nato_phonetic_alphabet.csv')
    # for (index, row) in df.iterrows():
    #     print(row.code)

    nato_dic = {row.letter: row.code for (index, row) in df.iterrows()}
    print(nato_dic)

ENTER = True
while ENTER:
    generate()
    res = input('ginada').lower()
    if res == 'exit':
        ENTER = False
