'''
writer : yizimi - yuanxin
Instructor : Mr. Mao, Place of Tang Dynasty and CITers
'''

import tkinter as tk
import numpy as np
from tkinter.messagebox import showinfo

N = 9
layout = [[0 for j in range(N)] for k in range(N)]

root = tk.Tk()
root.title("数独游戏")

frametop = tk.Frame(root)
gridvar = [[tk.StringVar() for column in range(N)] for row in range(N)]
frame = [tk.Frame(frametop) for row in range(N)]
grid = [[tk.Button(frame[row], width = 3, textvariable = gridvar[row][column],\
     relief = tk.GROOVE, command = lambda row = row,\
          column = column:gridclick(row, column),\
          font = ('Helvetica', '12')) for column in range(N)] for row in range(N)]
for row in range(N):
    for column in range(N):
        grid[row][column].pack(side = tk.LEFT)
    frame[row].pack(side = tk.TOP)
frametop.pack(side=tk.TOP, pady = 10)

framemiddle = tk.Frame(root)
selections = [tk.Button(framemiddle, width = 3, text = '%d' % number, relief = tk.RAISED,\
     command = lambda number=number:numberclick(selections[number - 1]),\
          font = ('Helvetica', '12')) for number in range(1, 10)]
for each in selections:
    each.pack(side = tk.LEFT)
framemiddle.pack(side=tk.TOP, pady = 15)


framebottom = tk.Frame(root)
erase = tk.Button(framebottom, text = '删除', relief = tk.RAISED,\
     font = ('HeiTi', '14', 'bold'), width = 7, height = 1, bg = 'darkgreen', fg = 'white')
erase.pack(side = tk.LEFT, padx = 15)
check = tk.Button(framebottom, text = '核查', relief = tk.RAISED,\
     font = ('HeiTi', '14', 'bold'), width = 7, height = 1, bg = 'darkblue', fg = 'white')
check.pack(side = tk.LEFT, padx = 15)
ok = tk.Button(framebottom, text = '退出', relief = tk.RAISED, command = exit,\
     font = ('HeiTi', '14', 'bold'), width = 7, height = 1, bg = 'darkred', fg = 'white')
ok.pack(side = tk.LEFT, padx = 15)
framebottom.pack(side = tk.TOP, pady = 5)

def gridclick(row, column):
    number = ''
    for i in range(N):
        if selections[i]["relief"] == tk.SOLID:
            number = '%d' % (i + 1)
            break
    gridvar[row][column].set(number)
    if number == '':
        layout[row][column] = 0
    else :
        layout[row][column] = int(number)


def numberclick(selectionbutton):
    for i in range(N):
        selections[i]['relief'] = tk.RAISED
    erase["relief"] = tk.RAISED
    selectionbutton["relief"] = tk.SOLID

def eraseclick(event):
    for i in range(N):
        selections[i]['relief'] = tk.RAISED


def checkclick(event):
    correct = verify() 
    if correct:
        showinfo("核查结果", "答案正确")
        print("correct")
    else:
        showinfo("核查结果", "答案不正确")
        print("wrong")
check.bind("<Button-1>", checkclick)

def readlayout(filename):
    layoutfile = open(filename, 'r')
    lines = layoutfile.readlines()
    for row in range(N):
        line = lines[row].strip('\n')
        if line != '':
            for i in range(len(line)):
                if line[i] != '' and line[i] != ' ':
                    layout[row][i] = int(line[i])
    return layout

def showlayout(layout, gridvar):
    for row in range(N):
        for column in range(N):
            gridvar[row][column].set(str(layout[row][column]) if (layout[row][column] != 0) else '')
            if layout[row][column] != 0:
                grid[row][column]["state"] = tk.DISABLED


def verifyrow():
    correct = True
    for row in range(N):
        line = layout[row].copy()
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct

def verifycolumn():
    correct = True
    for column in range(N):
        line = list((np.array(layout))[:, column])
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct

def verifyblock():
    correct = True
    for blockindex in range(N):
        block = getblock(blockindex)
        line = list((np.array(layout))[block[0] : block[1] + 1, block[2] : block[3] + 1].reshape(N))
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct


def getblock(index):
    rowstart = index // 3 * 3
    rowend = rowstart + 2
    columnstart = index % 3 * 3
    columnend = columnstart + 2
    return rowstart, rowend, columnstart, columnend

def verify():
    return verifyrow() & verifycolumn() & verifyblock()

erase.bind("<Button-1>", eraseclick)

layout = readlayout('sudoku.txt')
showlayout(layout, gridvar)

tk.mainloop()