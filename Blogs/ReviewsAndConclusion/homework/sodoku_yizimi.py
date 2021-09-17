'''
writer : yizimi - yuanxin
Instructor : Mr. Mao, Place of Tang Dynasty and CITers
'''

import tkinter as tk
import numpy as np
from tkinter.messagebox import showinfo
import copy

N = 9
layout = [[0 for j in range(N)] for k in range(N)]
# ======================================================
fix = [[False for j in range(N)] for k in range(N)]
# ======================================================

root = tk.Tk()
root.title("数独游戏")

frametop = tk.Frame(root)
gridvar = [[tk.StringVar() for column in range(N)] for row in range(N)]
frame = [tk.Frame(frametop) for row in range(N)]
grid = [[tk.Button(frame[row], width=3, textvariable=gridvar[row][column],\
     relief=tk.GROOVE, command=lambda row=row,\
          column=column:gridclick(row, column),\
          font=('Helvetica', '12')) for column in range(N)] for row in range(N)]
for row in range(N):
    for column in range(N):
        grid[row][column].pack(side=tk.LEFT)
    frame[row].pack(side=tk.TOP)
frametop.pack(side=tk.TOP, pady=10)

framemiddle = tk.Frame(root)
selections = [tk.Button(framemiddle, width=3, text='%d' % number, relief=tk.RAISED,\
     command=lambda number=number:numberclick(selections[number - 1]),\
          font=('Helvetica', '12')) for number in range(1, 10)]
for each in selections:
    each.pack(side=tk.LEFT)
framemiddle.pack(side=tk.TOP, pady=15)

framebottom = tk.Frame(root)
erase = tk.Button(framebottom, text='删除', relief=tk.RAISED,\
     font=('HeiTi', '14', 'bold'), width=6, height=1, bg='darkgreen', fg='white')
erase.pack(side=tk.LEFT, padx=5)
check = tk.Button(framebottom, text='核查', relief=tk.RAISED,\
     font=('HeiTi', '14', 'bold'), width=6, height=1, bg='darkblue', fg='white')
check.pack(side=tk.LEFT, padx=5)
ans = tk.Button(framebottom, text="答案", relief=tk.RAISED,\
     font=('HeiTi', '14', "bold"), width=6, height=1, bg='black', fg='white')
ans.pack(side=tk.LEFT, padx=5)
ok = tk.Button(framebottom, text='退出', relief=tk.RAISED, command=exit,\
     font=('HeiTi', '14', 'bold'), width=6, height=1, bg='darkred', fg='white')
ok.pack(side=tk.LEFT, padx=5)
framebottom.pack(side=tk.TOP, pady=5)


def gridclick(row, column):
    number = ''
    for i in range(N):
        if selections[i]["relief"] == tk.SOLID:
            number = '%d' % (i + 1)
            break
    gridvar[row][column].set(number)
    if number == '':
        layout[row][column] = 0
    else:
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
    correct = verify(layout)
    if correct:
        showinfo("核查结果", "答案正确")
        print("correct")
    else:
        showinfo("核查结果", "答案不正确")
        print("wrong")
    # check['relief'] = tk.RAISED


check.bind("<Button-1>", checkclick)


# ========================================================
def ansclick(event):
    showlayout(solve(layout, fix), fix, gridvar)


ans.bind("<Button-1>", ansclick)
# ========================================================


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


def showlayout(layout, fix, gridvar):
    for row in range(N):
        for column in range(N):
            gridvar[row][column].set(
                str(layout[row][column]) if (layout[row][column] != 0) else '')
            if fix[row][column]:
                grid[row][column]["state"] = tk.DISABLED


def verifyrow(layout):
    correct = True
    for row in range(N):
        line = layout[row].copy()
        # print(line)
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct


def verifycolumn(layout):
    correct = True
    for column in range(N):
        line = list((np.array(layout))[:, column])
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct


def verifyblock(layout):
    correct = True
    for blockindex in range(N):
        block = getblock(blockindex)
        line = list((np.array(layout))[block[0]:block[1] + 1,
                                       block[2]:block[3] + 1].reshape(N))
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


# =====================必要性添加========================


# def verify():
#     return verifyrow() & verifycolumn() & verifyblock()
def verify(layout):
    if not verifyrow(layout):
        return False
    elif not verifycolumn(layout):
        return False
    elif not verifyblock(layout):
        return False
    else:
        return True


def findfix(layout):
    fix = [[False for j in range(N)] for k in range(N)]
    for row in range(N):
        for column in range(N):
            if layout[row][column] > 0:
                fix[row][column] = True
    return fix


# =======================================================

# =======================暴力枚举法=======================

# def getnext(layout, fix):
#     nextchanged = False
#     for row in range(N - 1, -1, -1):
#         for column in range(N - 1, -1, -1):
#             if (not fix[row][column]):
#                 layout[row][column] += 1
#                 if layout[row][column] > 9:
#                     layout[row][column] = 1
#                 else :
#                     nextchanged = True
#                     break
#         if nextchanged:
#             break
#     return layout

# def solve(layout, fix):
#     for row in range(N):
#         for column in range(N):
#             if (not fix[row][column]):
#                 layout[row][column] = 1
#     startgrid = [0, 0]
#     for index in range(N * N):
#         if not fix[index // N][index % N]:
#             startgrid = [index // N, index % N]
#             break

#     solved = False
#     while (not solved) & (layout[startgrid[0]][startgrid[1]] <= 9):
#         if verify():
#             solved = True
#             break
#         else:
#             layout = getnext(layout, fix)
#     return layout

# ==========================================================

# ========================枚举剪枝法=========================


def findpossiblevalues(layout, fix):
    values = copy.deepcopy(layout)
    for row in range(N):
        for column in range(N):
            if (not fix[row][column]):
                values[row][column] = [i + 1 for i in range(N)]
            else:
                values[row][column] = [layout[row][column]]
    for row in range(N):
        for column in range(N):
            if fix[row][column]:
                # print(values[row][column])
                number = values[row][column][0]
                for i in range(N):
                    if (i != column) & (number in values[row][i]):
                        values[row][i].remove(number)
                for i in range(N):
                    if (i != row) & (number in values[i][column]):
                        values[i][column].remove(number)
                block = getblock_ij(row, column)
                for i in range(block[0], block[1] + 1):
                    for j in range(block[2], block[3] + 1):
                        if ((i != row) |
                            (j != column)) & (number in values[i][j]):
                            values[i][j].remove(number)
    return values


def getblock_ij(row, column):
    rowstart = row // 3 * 3
    rowend = rowstart + 2
    columnstart = column // 3 * 3
    columnend = columnstart + 2
    return rowstart, rowend, columnstart, columnend


def getnext(layout, fix, values):
    nextchanged = False
    for row in range(N - 1, -1, -1):
        for column in range(N - 1, -1, -1):
            if (not fix[row][column]):
                currentvalues = values[row][column]
                i = currentvalues.index(layout[row][column])
                if i < (len(currentvalues) - 1):
                    layout[row][column] = currentvalues[i + 1]
                    nextchanged = True
                    break
                else:
                    layout[row][column] = currentvalues[0]
        if nextchanged:
            break
    return nextchanged, layout


def solve(layout, fix):
    values = findpossiblevalues(layout, fix)

    for row in range(N):
        for column in range(N):
            if not fix[row][column]:
                layout[row][column] = values[row][column][0]

    solved = False
    nextchanged = True
    while (not solved) & (nextchanged):
        if verify(layout):
            solved = True
            break
        else:
            nextchanged, layout = getnext(layout, fix, values)

    return layout


# ==========================================================

erase.bind("<Button-1>", eraseclick)

layout = readlayout('solution_revised.txt')
fix = findfix(layout)

showlayout(layout, fix, gridvar)

tk.mainloop()
