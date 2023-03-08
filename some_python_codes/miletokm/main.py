# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
from tkinter import *

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
def button_clicked():
    print("I got clicked")

    converted = float(input.get())*1.8
    my_label2.config(text=str(round(converted,3)))

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    window = Tk()
    window.title("Mile to Km")
    window.minsize(width=100, height=50)
    window.config(padx=10, pady=10)

    # Label
    my_label = Label(text="Mile", font=("Arial", 24, "bold"))
    my_label.config(text="Mile")
    my_label.grid(column=2, row=0)
    #my_label.config(padx=50, pady=50)

    my_label2 = Label(text="  ", font=("Arial", 24))
    # my_label.config(text="Mile")
    my_label2.grid(column=1, row=1)
    #my_label2.config(padx=50, pady=50)

    my_label3 = Label(text="is equal to ", font=("Arial", 18))
    # my_label.config(text="Mile")
    my_label3.grid(column=0, row=1)
    my_label3.config(padx=50, pady=50)

    my_label1 = Label(text="KM", font=("Arial", 24, "bold"))
    #my_label.config(text="Mile")
    my_label1.grid(column=2, row=1)
    my_label1.config(padx=50, pady=50)

    # Button
    button = Button(text="Calculate", command=button_clicked)
    button.grid(column=1, row=2)

    input = Entry(width=10)
    print(input.get())
    input.grid(column=1, row=0)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
    window.mainloop()