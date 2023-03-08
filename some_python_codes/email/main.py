# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import smtplib
import random
import datetime as dt
import pandas as pd
import fileinput

email = 'kaveh.azizian20@gmail.com'
password = 'Kian+laval559@'
text_list = ['letter_1.txt', 'letter_2.txt', 'letter_3.txt']


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


def send_birthday(name):
    mytext = random.choice(text_list)
    try:
        with open(f'letter_templates/{mytext}') as txt:
            lines = txt.readlines()
            file = open('letterto_send.txt', 'w')
            for line in lines:
                linetowrite = line.replace('[NAME]', name)

                file.writelines(linetowrite)
            file.close()
    except FileNotFoundError:
        print('Sample text letter not found!')

    with open('letterto_send.txt', 'r') as to_send:
        all_lines = to_send.read()
        with smtplib.SMTP("smtp.gmail.com") as sm:
            sm.starttls()
            sm.login(email, password)
            sm.sendmail(from_addr=email,to_addrs=email,msg=f'Subject: HappyBirthDay! \n\n{all_lines}')

def send_email():
    if weekday == 1:
        with open('quotes.txt') as quote_file:
            all_quotes = quote_file.readlines()
            quote = random.choice(all_quotes)
            print(quote)
        with smtplib.SMTP("smtp.gmail.com") as sm:
            sm.starttls()
            sm.login(email, password)
            sm.sendmail(from_addr=email, to_addrs=email, msg=f"Subject:Hi There \n\n{quote}")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    now = dt.datetime.now()
    weekday = now.weekday()
    today = now.day
    this_month = now.month

try:
    df = pd.read_csv('birthdays.csv')
except FileNotFoundError:
    print('No File found')
else:
    day = df['day']

    #print(df)
    for i, j in day.items():
        print(i)
        print(j)
        if (j == today) and (df['month'][i] == this_month):

            name = df['name'][i]
            print(name)

            send_birthday(name)

# send_birthday('lk')
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
