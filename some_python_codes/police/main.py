# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import pathlib
import os
import datetime
import pandas as pd
def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
try:
    df=pd.read_csv('../../Desktop/python_chini/pol/police/Deaths_by_Police_US.csv')
except FileNotFoundError:
    path=os.path.abspath("../../Desktop/python_chini/pol/police/Deaths_by_Police_US.csv")
    print(f'File not found in {path}')
else:

     kj=pathlib.Path(__file__).parent.absolute()

     print(kj)
finally:
    manner=df['manner_of_death'].unique()
    #print(manner)
    race=manner=df['race'].value_counts()
    #print((race))
    #print(df['city'].isnull)
   # print(df.shape)
    null_data=df.isnull
    #print(null_data)
    cldata=df['race'].dropna()
    #print(cldata.shape)
    mental = df['city'].value_counts()
   # print(mental)
    age= df[(df['age'] <=30) & (df['age'] >=20)]
    #print(age['age'].value_counts())
    print(df['age'].mean())
    print(age['gender'].value_counts())
    df['date']=pd.to_datetime(df['date'])
    date=df.sort_values('date')

    x = datetime.datetime(2016, 1, 1)
    #print(date[date['date']< x]['date'])
    state_dat=df.groupby('state')
    print(state_dat['date'])