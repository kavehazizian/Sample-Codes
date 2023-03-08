# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


def change_name(name ):
    path = "Mail Merge Project Start/Input/Letters/starting_letter.txt"
    new_Path = 'Mail Merge Project Start\Output\ReadyToSend\email_'+name+'.txt'
    orig_mail = open(path, 'r')
    new_mail = open(new_Path, "w+")
    for lines in orig_mail:
        x=lines.replace("[name]", name)
        new_mail.write(x)

    new_mail.close()

    orig_mail.close()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    # myfile = open("Mail Merge Project Start/Input/Letters/starting_letter.txt", "a+")
    # print(myfile.readline()[5::])
    # myfile.close()
    name = open("Mail Merge Project Start/Input/Names/invited_names.txt", "r")
    for line in name:
        change_name(line.rstrip('\n'))

    name.close()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
