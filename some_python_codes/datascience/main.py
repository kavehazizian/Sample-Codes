# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import pandas as pd


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    df = pd.read_csv('salaries_by_college_major.csv')
    #print(df.shape)
    df.isna()
    clean_df = df.dropna()
    print(clean_df['Starting Median Salary'].max())
    print(clean_df['Starting Median Salary'].idxmax())
    #print(clean_df.iloc[43])
    print(clean_df['Mid-Career Median Salary'].max())
    print(clean_df['Mid-Career Median Salary'].idxmax())
    print(clean_df['Mid-Career Median Salary'][8])

    print(clean_df['Mid-Career Median Salary'].min())
    print(clean_df['Mid-Career Median Salary'].idxmin())
    print(clean_df.iloc[18])

    new_df = clean_df['Mid-Career 90th Percentile Salary'] - clean_df['Mid-Career 10th Percentile Salary']
    clean_df.insert(1, 'spread', new_df)
    #print(clean_df.head(5))
    clean_df.to_csv('new.csv',index=False)
    dfs=clean_df.sort_values(by=['Starting Median Salary'],ascending=False)
    #print(dfs.head(5))
    dfs.to_csv('sorted.csv',index=False)
    gh = clean_df.sort_values(by=['spread'], ascending=False)
    print(gh[['Starting Median Salary','Undergraduate Major', 'spread']].head())

    print(clean_df.groupby('Group').mean())
    #print(df)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
