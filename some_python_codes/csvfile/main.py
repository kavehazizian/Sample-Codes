# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import csv
import pandas
import turtle


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


def get_mouse_click_coor(x, y):
    print(x, y)




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #
    data = pandas.read_csv('2018_Central_Park_Squirrel_Census_-_Squirrel_Data.csv')
    print(data[data['Primary Fur Color'] == 'Cinnamon'])
    gary_sqirls_count = len(data[data['Primary Fur Color'] == 'Gray'])
    cinnamon_sqirls_count = len(data[data['Primary Fur Color'] == 'Cinnamon'])
    black_sqirls_count = len(data[data['Primary Fur Color'] == 'Black'])
    data_dic = {"Fur color": ["Gray", "Cinnamon", "Black"], "count": [gary_sqirls_count, cinnamon_sqirls_count,
                                                                      black_sqirls_count]}
    df = pandas.DataFrame(data_dic)
    df.to_csv('SquirlsCount.csv')
    # See PyCharm help at https://www.jetbrains.com/help/pycharm/
    screen = turtle.Screen()
    screen.title('U.S. States')
    image = "us-states-game-start/blank_states_img.gif"
    # image = "eglite.gif"
    # image = "us-states-game-start/sample.jpg"

    screen.addshape(image)
    turtle.shape(image)
    # turtle.onscreenclick(get_mouse_click_coor)
    # turtle.mainloop()
    state_list=[]
    while len(state_list) < 50:
         answer = screen.textinput(f'{len(state_list)}/50 correct', 'Enter a state name').title()
    # screen.exitonclick()
         map = pandas.read_csv('us-states-game-start/50_states.csv')
         all_state = map.state.to_list()
         if answer=="Exit":
            missing_states=[]
            for stat in all_state:
                if stat not in state_list:
                    missing_states.append(stat)
            df=pandas.DataFrame(missing_states)
            df.to_csv('states_to_learn.csv')
            break

         if answer in all_state:
             state_list.append(answer)
             t = turtle.Turtle()
             t.hideturtle()
             t.penup()
             state_data = map[map.state == answer]
             t.goto(int(state_data.x), int(state_data.y))
             t.write(answer)
