# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import pandas as pd
import matplotlib.pyplot as plt

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/

df=pd.read_csv('QueryResults.csv', names=['DATE', 'TAG', 'POSTS'], header=0)
#print(df.shape)
#print(df.count())
gb=df.groupby('TAG')
#print(gb.count())
print(df['DATE'][1])
df['DATE']=pd.to_datetime(df['DATE'])
#print(df.head(5))
reshaped_df=df.pivot(index='DATE',columns='TAG',values='POSTS')
#print(reshaped_df.head(6))
plt.plot(reshaped_df.index, reshaped_df['c++'])