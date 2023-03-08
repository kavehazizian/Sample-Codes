# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
from tkinter import *
import pandas as pd
import random

BACKGROUND_COLOR = "#B1DDC6"

try:
    df = pd.read_csv('data/Reaming.csv')
except FileNotFoundError:
    df = pd.read_csv('data/french_words.csv')

to_learn = df.to_dict(orient='records')
rand_num = 0
learnt_list = set()
rand_range = 99


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


def isKnown():
    global learnt_list, rand_range
    learnt_list.add(rand_num)
    pick_word()
    listto_learn = df.drop(learnt_list)
    rand_range= len(listto_learn )

    listto_learn.to_csv('data/Reaming.csv', index=False)


def pick_word():
    global rand_num

    rand_num = random.randint(0, rand_range)
    print(rand_num)
    french_word = df['French'][rand_num]
    print(french_word)
    canvas.itemconfig(card_title, text='French', fill='black')
    canvas.itemconfig(card_word, text=french_word, fill='black')
    canvas.itemconfig(card_background, image=cardfront_image)
    flip_timer = window.after(3000, flip_card)


def chose_word():
    current_card = random.choice(to_learn)
    print(current_card)


def flip_card():
    window.after_cancel(flip_timer)
    english_word = df['English'][rand_num]
    canvas.itemconfig(card_title, text='English', fill='white')
    canvas.itemconfig(card_word, text=english_word, fill='white')
    canvas.itemconfig(card_background, image=cardback_image)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    window = Tk()
    window.title('Flashy')
    window.config(padx=50, pady=50, bg=BACKGROUND_COLOR)

    flip_timer = window.after(3000, flip_card)
    canvas = Canvas()
    canvas.config(width=800, height=526)
    cardfront_image = PhotoImage(file="images/card_front.png")
    cardback_image = PhotoImage(file='images/card_back.png')
    card_background = canvas.create_image(400, 263, image=cardfront_image)
    # canvas.create_image(400, 263, image=cardback_image)
    card_title = canvas.create_text(400, 150, text='Title', font=('Ariel', 40, 'italic'))
    card_word = canvas.create_text(400, 263, text='Word', font=('Ariel', 60, 'bold'))
    canvas.config(bg=BACKGROUND_COLOR, highlightthickness=0)
    canvas.grid(row=0, column=0, columnspan=2)

    cross_image = PhotoImage(file='images/wrong.png')
    unknown_button = Button(image=cross_image, highlightthickness=0, command=pick_word)

    unknown_button.grid(row=1, column=0)
    check_image = PhotoImage(file='images/right.png')
    correct_button = Button(image=check_image, highlightthickness=0, command=isKnown)
    correct_button.grid(row=1, column=1)

    # pick_word()

    window.mainloop()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
