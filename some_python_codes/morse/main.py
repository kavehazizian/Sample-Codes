# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import pandas as pd
import pathlib
print('current path')
print(pathlib.Path().absolute())
df=pd.read_csv('morse\morse.csv', encoding='utf-8',error_bad_lines=False)
def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.

def convert_toMorse(msg):
    print(type(msg))
    umsg=msg.upper()

    str=[]
    for char in umsg:
        try:
           ms=df.loc[df['word'] == char,'morse'].values[0]
           print(ms)

        except IndexError:
            print(f'sishdin {char} not found!')
        else:
            str.append(ms)
    return tuple(str)
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
       str=convert_toMorse('Kaveh, Azizian')

       print(str)
#print(df[df['word']=='D'])