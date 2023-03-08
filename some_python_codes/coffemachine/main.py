# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import time

MENU = {
    "espresso": {
        "ingredients": {
            "water": 50,
            "coffee": 18,
        },
        "cost": 1.5,
    },
    "latte": {
        "ingredients": {
            "water": 200,
            "milk": 150,
            "coffee": 24,
        },
        "cost": 2.5,
    },
    "cappuccino": {
        "ingredients": {
            "water": 250,
            "milk": 100,
            "coffee": 24,
        },
        "cost": 3.0,
    }
}
resources = {
    "water": 300,
    "milk": 200,
    "coffee": 100,
}


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


def order():
    return input('What would you like? \n "L" for Latte \n "C" for Cappuccino \n "E" for Espresso \n "R" for Report \n '
                 '"O" for turning off \n')


def report():
    print('Here is the current stock ')
    for res, val in resources.items():
        print(res + ':' + str(val))


def check_coin(name='', inserted_coin=0.0):
    result = True
    refund = 0.0
    for adi, x in MENU.items():
        if adi == name:
            for z, y in x.items():
                if z == 'cost' and y > inserted_coin:
                    print('Not enough coins!')
                    result = False
                elif z == 'cost' and y <= inserted_coin:
                    refund = inserted_coin - y

    return [result, refund]


def insert_coin():
    quarters = input("How many quarters")
    dimes = input("How many dimes")
    nickels = input("How many nickels")
    pennies = input("How many pennies")
    total = int(quarters) * 0.25 + int(dimes) * 0.1 + int(nickels) * 0.05 + int(pennies) * 0.01
    return total


def check_resource(name=''):
    result = True
    for adi, x in MENU.items():
        if adi == name:
            for z, y in x.items():
                if not z == 'cost':
                    for item, val in y.items():
                        if item == 'water' and val > resources['water']:
                            print("No sufficient water")
                            result = False
                            break
                        elif item == 'milk' and val > resources['milk']:
                            print("No sufficient milk")
                            result = False
                            break
                        elif item == 'coffee' and val > resources['coffee']:
                            print("No sufficient coffee")
                            result = False
    return result


def operate(name=''):
    tot = insert_coin()
    res = check_coin(name, tot)
    if res[0]:
        print('here is your {} coffee and your refund of {}'.format(name, res[1]))


def print_dic():
    for x in MENU.values():
        for z, y in x.items():
            if not z == 'cost':
                print(y)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # order()
    # print_dic()
    machine_off = False
    # check_resource('cappuccino')
    while not machine_off:
        time.sleep(2.5)
        res = order()

        if res.lower() == 'o':
            machine_off = True
            print("Machine is off")
        elif res.lower() == 'r':
            report()
        elif res.lower() == 'l':
            if check_resource('latte'):
                # tot = insert_coin()
                # res = check_coin('latte', tot)
                # if res[0]:
                #     print('here is your coffee and your refund of {}'.format(res[1]))
                operate('latte')
        elif res.lower() == 'c':
            if check_resource('cappuccino'):
                operate('cappuccino')
        elif res.lower() == 'e':
            if check_resource('espresso'):
                operate('espresso')
        else:
            print('Not a correct option')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
